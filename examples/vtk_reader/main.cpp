#include <torch/torch.h>
#include <string>

#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/data_set.h"
#include "irl/ml_classification/net.h"

#include <vtkSmartPointer.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDataSet.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkInformation.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkRectilinearGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <iostream>
#include <vector>

int main() {
    auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
    reader->SetCaseFileName("/home/quirin/mlcfd/Repositories/jet/nga.case");
    reader->UpdateInformation();

    vtkInformation* info = reader->GetOutputInformation(0);
    if (!info->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
        std::cerr << "No time steps found." << std::endl;
        return 1;
    }

    int nSteps = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    std::vector<double> timeSteps(nSteps);
    info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), timeSteps.data());

    // Set to last time step
    reader->SetTimeValue(timeSteps[nSteps - 1]);
    reader->Update();

    // Get the output as multiblock
    auto mb = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
    if (!mb) {
        std::cerr << "Output is not a vtkMultiBlockDataSet." << std::endl;
        return 1;
    }

    // Assume data is in the first block (adjust if needed)
    auto block = mb->GetBlock(0);
    if (!block) {
        std::cerr << "No block found in multiblock dataset." << std::endl;
        return 1;
    }

    auto dataSet = vtkDataSet::SafeDownCast(block);
    if (!dataSet) {
        std::cerr << "Block is not a vtkDataSet." << std::endl;
        return 1;
    }

    // Access cell data (volume fraction is cell-based)
    vtkCellData* cellData = dataSet->GetCellData();
    if (!cellData) {
        std::cerr << "No cell data found." << std::endl;
        return 1;
    }

    vtkDataArray* vofArray = cellData->GetArray("VOF");
    if (!vofArray) {
        std::cerr << "VOF array not found in cell data." << std::endl;
        return 1;
    }

    // Cast to vtkRectilinearGrid to access dimensions
    auto grid = vtkRectilinearGrid::SafeDownCast(dataSet);
    if (!grid) {
        std::cerr << "Data is not a vtkRectilinearGrid." << std::endl;
        return 1;
    }

    int dims[3];  // point dimensions
    grid->GetDimensions(dims);

    // Convert point dims to cell dims (cells = points - 1)
    int nx = dims[0] - 1;
    int ny = dims[1] - 1;
    int nz = dims[2] - 1;

    std::cout << "Grid cell dimensions:" << std::endl;
    std::cout << "  nx = " << nx << std::endl;
    std::cout << "  ny = " << ny << std::endl;
    std::cout << "  nz = " << nz << std::endl;

    const int stencil_size = 3;

    int no_filled_cells = 0;

    // Loop through interior cells (excluding 1-cell boundary)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int k = 1; k < nz - 1; ++k) {

                vtkIdType centerCellId = i + j * nx + k * nx * ny;

                // Skip this cell if VOF = 0
                double vof_center = vofArray->GetComponent(centerCellId, 0);
                if (vof_center == 0) {
                    continue;
                }

                no_filled_cells ++;

                // 3D array to store VOF values in the 3x3x3 neighborhood
                double vfrac[stencil_size][stencil_size][stencil_size];

                // Fill the 3x3x3 neighborhood
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        for (int dk = -1; dk <= 1; ++dk) {
                            int ni = i + di;
                            int nj = j + dj;
                            int nk = k + dk;

                            vtkIdType cellId = ni + nj * nx + nk * nx * ny;
                            if (cellId >= vofArray->GetNumberOfTuples()) {
                                std::cerr << "Invalid cellId: " << cellId << std::endl;
                                continue;
                            }

                            vfrac[di + 1][dj + 1][dk + 1] = vofArray->GetComponent(cellId, 0);
                        }
                    }
                }

                // Flatten the 3D vector into a 1D vector
                std::vector<double> flattened_vfrac;
                for (int ii = 0; ii < stencil_size; ++ii) {
                    for (int jj = 0; jj < stencil_size; ++jj) {
                        for (int kk = 0; kk < stencil_size; ++kk) {
                            flattened_vfrac.push_back(vfrac[ii][jj][kk]);
                        }
                    }
                }


                // You could store flattened_vfrac here into a larger structure or dataset if needed.
            }
        }
    }
    std::cout << no_filled_cells << std::endl;
    return 0;
}


/*

#include <vtkSmartPointer.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkFieldData.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkRectilinearGrid.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>


int main (int argc, char* argv[]) {

  // To test if Torch works
  torch::Tensor tensor = torch::eye(3);
  std::cout << tensor << std::endl;


  // vtk reader


  auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
  reader->SetCaseFileName("/home/mlcfd/Repositories/jet/nga.case");  // Adjust the path as needed
  reader->UpdateInformation();  // Populate meta info like time steps
  vtkInformation* info = reader->GetOutputInformation(0);
  if (info->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
      int nSteps = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
      std::vector<double> timeSteps(nSteps);
      info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0]);


      std::cout << "Number of timesteps: " << nSteps << std::endl;
      
      
      nSteps = 50;

      for (int t = 49; t < nSteps; ++t) {
          reader->SetTimeValue(timeSteps[t]);
          reader->Update();

          auto output = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
          if (!output) {
              std::cerr << "Reader did not return a vtkMultiBlockDataSet." << std::endl;
              return 1;
          }
          
          std::cout << "Number of blocks: " << output->GetNumberOfBlocks() << std::endl;

          for (unsigned int i = 0; i < output->GetNumberOfBlocks(); ++i) {
              auto block = output->GetBlock(i);
              if (block) {
                  std::cout << "Block " << i << " is type: " << block->GetClassName() << std::endl;
              } else {
                  std::cout << "Block " << i << " is nullptr." << std::endl;
              }
          }

          auto block0 = vtkRectilinearGrid::SafeDownCast(reader->GetOutput()->GetBlock(0));
          if (!block0) {
              std::cerr << "Could not access block 0 of the dataset." << std::endl;
              return 1;
          }

          const auto dimensions = block0->GetDimensions();
          const int ncellsx = dimensions[0] - 1;
          const int ncellsy = dimensions[1] - 1;
          const int ncellsz = dimensions[2] - 1;

          std::cout << "Dimensions of grid are: " << ncellsx << ", " << ncellsy << ", " << ncellsz << std::endl;

          auto px = block0->GetXCoordinates();
          auto py = block0->GetYCoordinates();
          auto pz = block0->GetZCoordinates();
          std::vector<double> coordsx(ncellsx + 1, 0.0);
          std::vector<double> coordsy(ncellsy + 1, 0.0);
          std::vector<double> coordsz(ncellsz + 1, 0.0);
          for (int i = 0; i < ncellsx + 1; i++) coordsx[i] = px->GetComponent(i, 0);
          for (int j = 0; j < ncellsy + 1; j++) coordsy[j] = py->GetComponent(j, 0);
          for (int k = 0; k < ncellsz + 1; k++) coordsz[k] = pz->GetComponent(k, 0);

          auto cellData = block0->GetCellData();
          int numArrays = cellData->GetNumberOfArrays();
          std::cout << "Available cell data arrays: " << numArrays << std::endl;
          for (int i = 0; i < numArrays; ++i) {
              std::cout << "Array " << i << ": " << cellData->GetArrayName(i) << std::endl;
          }

          vtkDataArray* vofArray = block0->GetCellData()->GetArray("VOF");
          if (!vofArray) {
              std::cerr << "VOF array not found in the dataset." << std::endl;
              return 1;
          }
          std::cout << "Number of tuples = " << vofArray->GetNumberOfTuples() << std::endl;
          std::cout << "VOF values (first 10):" << std::endl;
          for (int i = 0; i < 10; ++i) {
              std::cout << vofArray->GetComponent(i, 0) << " ";
          }
          std::cout << std::endl;
      
          auto vfrac = std::vector<std::vector<std::vector<double>>>(ncellsx, std::vector<std::vector<double>>(ncellsy, std::vector<double>(ncellsz, 0.0)));

          for (int i = 0; i < ncellsx; ++i) {
              for (int j = 0; j < ncellsy; ++j) {
                  for (int k = 0; k < ncellsz; ++k) {
                      const int index = k * (ncellsy * ncellsx) + j * ncellsx + i;
                      vfrac[i][j][k] = static_cast<double>(vofArray->GetComponent(index, 0));
                      // Each z-slice has ncellsy * ncellsx elements, move forward k slices.
                      // Each row has ncellsx elements, so this skips over j rows
                      // Move to the correct column in a row
                  }
              }
          }

          WriteField(coordsx, coordsy, coordsz, vfrac, std::string("vfrac_")+std::to_string(t));
      
      }
  
  } else {
      std::cerr << "No time steps found in dataset." << std::endl;
  }

  return 0;
}
*/