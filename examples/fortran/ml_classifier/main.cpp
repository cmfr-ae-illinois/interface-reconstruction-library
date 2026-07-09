/*
 * Reads VTK data of gas-liquid flow simulation and classifies the interface
 * type using the torch-free Fortran implementation of the ML-classifier.
 */
#include <vtkCellCenters.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDataArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "examples/fortran/ml_classifier/ml_classifier_c_api.h"
#include "irl/generic_cutting/generic_cutting.h"

#include <iostream>
#include <vector>

static const int N = 5;

inline std::string classid_to_string(const int id) {
  if (id == 1) {
    return "ligament";
  } else if (id == 2) {
    return "droplet";
  } else if (id == 3) {
    return "sheet/film";
  } else if (id == 4) {
    return "ligament end";
  } else if (id == 5) {
    return "sheet end";
  } else if (id == 0) {
    return "well-resolved";
  }
  return "unknown";
}

int main(int argc, char* argv[]) {
  std::array<double, N * N * N> vfrac;
  std::array<double, N * N * N * 3> liq_bary;

  std::string datafile =
      "../examples/ml_minimal_classifier_notorch/test_data/data.vtr";
  std::string plicfile =
      "../examples/ml_minimal_classifier_notorch/test_data/plic.vtu";

  auto data_reader = vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
  data_reader->SetFileName(datafile.c_str());
  data_reader->Update();

  vtkRectilinearGrid* grid = data_reader->GetOutput();
  vtkCellData* cellData = grid->GetCellData();
  vtkDataArray* xCoords = grid->GetXCoordinates();
  vtkDataArray* yCoords = grid->GetYCoordinates();
  vtkDataArray* zCoords = grid->GetZCoordinates();
  vtkDataArray* vfrac_array =
      vtkDataArray::SafeDownCast(cellData->GetArray("VOF"));
  vtkDataArray* liq_bary_array =
      vtkDataArray::SafeDownCast(cellData->GetArray("liquid_bary"));
  int dims[3];
  grid->GetDimensions(dims);  // dims are for points, so cells are dims-1
  int cellDims[3] = {dims[0] - 1, dims[1] - 1, dims[2] - 1};
  auto interface_type = vtkSmartPointer<vtkIntArray>::New();
  interface_type->SetName("class");
  interface_type->SetNumberOfComponents(1);
  interface_type->SetNumberOfTuples(grid->GetNumberOfCells());
  for (int i = 0; i < grid->GetNumberOfCells(); i++) {
    interface_type->SetValue(i, -1);
  }

  std::cout << "Files read; processing data...\n";
  for (int i = 2; i < cellDims[0] - 2; ++i) {
    for (int j = 2; j < cellDims[1] - 2; ++j) {
      for (int k = 2; k < cellDims[2] - 2; ++k) {
        const int idx = i + j * cellDims[0] + k * cellDims[0] * cellDims[1];
        const double cell_vfrac = vfrac_array->GetComponent(idx, 0);
        if (cell_vfrac < IRL::global_constants::VF_LOW ||
            cell_vfrac > IRL::global_constants::VF_HIGH) {
          continue;
        }
        std::cout << "Processing cell "
                  << k + j * cellDims[2] + i * cellDims[2] * cellDims[1]
                  << " / " << cellDims[0] * cellDims[1] * cellDims[2] << "\r"
                  << std::flush;
        const double x0 = xCoords->GetComponent(i, 0);
        const double y0 = yCoords->GetComponent(j, 0);
        const double z0 = zCoords->GetComponent(k, 0);
        const double x1 = xCoords->GetComponent(i + 1, 0);
        const double y1 = yCoords->GetComponent(j + 1, 0);
        const double z1 = zCoords->GetComponent(k + 1, 0);
        const double dx = x1 - x0;
        const double dy = y1 - y0;
        const double dz = z1 - z0;
        const double cell_volume = std::abs(dx * dy * dz);
        IRL::Pt cell_center(0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1));
        for (int ii = 0; ii < 5; ++ii) {
          for (int jj = 0; jj < 5; ++jj) {
            for (int kk = 0; kk < 5; ++kk) {
              const int local_idx = (i + ii - 2) + (j + jj - 2) * cellDims[0] +
                                    (k + kk - 2) * cellDims[0] * cellDims[1];
              vfrac[ii + jj * N + kk * N * N] =
                  vfrac_array->GetComponent(local_idx, 0);
              IRL::Pt cell_liq_bary(liq_bary_array->GetComponent(local_idx, 0),
                                    liq_bary_array->GetComponent(local_idx, 1),
                                    liq_bary_array->GetComponent(local_idx, 2));
              cell_liq_bary -= cell_center;
              cell_liq_bary[0] /= dx;
              cell_liq_bary[1] /= dy;
              cell_liq_bary[2] /= dz;
              cell_liq_bary *= vfrac[ii + jj * N + kk * N * N];
              liq_bary[ii + jj * N + kk * N * N + 0 * N * N * N] =
                  cell_liq_bary[0];
              liq_bary[ii + jj * N + kk * N * N + 1 * N * N * N] =
                  cell_liq_bary[1];
              liq_bary[ii + jj * N + kk * N * N + 2 * N * N * N] =
                  cell_liq_bary[2];
            }
          }
        }
        interface_type->SetValue(
            idx, ml_classifier_fortran(vfrac.data(), liq_bary.data()));
      }
    }
  }

  int count_classes[6] = {0};
  for (int i = 0; i < grid->GetNumberOfCells(); i++) {
    int class_id = interface_type->GetValue(i);
    if (class_id >= 0 && class_id < 6) {
      count_classes[class_id]++;
    }
  }

  std::cout << "\nClass counts:\n";
  for (int i = 0; i < 6; ++i) {
    std::cout << "  " << classid_to_string(i) << "\t: " << count_classes[i]
              << "\n";
  }

  auto plic_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  plic_reader->SetFileName(plicfile.c_str());
  plic_reader->Update();
  vtkUnstructuredGrid* plic_grid = plic_reader->GetOutput();

  auto plic_interface_type = vtkSmartPointer<vtkIntArray>::New();
  plic_interface_type->SetName("class");
  plic_interface_type->SetNumberOfComponents(1);
  plic_interface_type->SetNumberOfTuples(plic_grid->GetNumberOfCells());
  for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
    plic_interface_type->SetValue(i, -1);
  }

  // Cell locator
  auto locator = vtkSmartPointer<vtkCellLocator>::New();
  locator->SetDataSet(grid);
  locator->BuildLocator();
  auto plic_cell_center_filter = vtkSmartPointer<vtkCellCenters>::New();
  plic_cell_center_filter->SetInputData(plic_grid);
  plic_cell_center_filter->Update();
  auto plic_cell_centers = plic_cell_center_filter->GetOutput()->GetPoints();

  // Convert interface_type to PLIC array
  for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
    auto type = interface_type->GetValue(
        locator->FindCell(plic_cell_centers->GetPoint(i)));
    plic_interface_type->SetValue(i, type);
  }

  // Write PLIC file
  plic_grid->GetCellData()->AddArray(plic_interface_type);
  auto plic_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  plic_writer->SetFileName("plic_with_classes.vtu");
  plic_writer->SetInputData(plic_grid);
  plic_writer->Write();

  return 0;
}