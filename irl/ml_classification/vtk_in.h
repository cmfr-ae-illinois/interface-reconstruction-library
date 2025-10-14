#include <iostream>
#include <vector>
#include <string>

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
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkAppendFilter.h>
#include <vtkCellLocator.h>
#include <vtkPolygon.h>
#include <vtkCellCenters.h>

namespace IRL {

void classify_simulation(IRL::Classifier& classifier, const std::string& filename) {
    auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
    reader->SetCaseFileName(filename.c_str());
    //reader->SetCaseFileName("/home/quirin/mlcfd/Repositories/jet/nga.case");
    //reader->SetCaseFileName("/home/quirin/jet/nga.case");
    reader->UpdateInformation();

    vtkInformation* info = reader->GetOutputInformation(0);
    int nSteps = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    std::vector<double> timeSteps(nSteps);
    info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), timeSteps.data());
    reader->SetTimeValue(timeSteps[nSteps - 1]);
    reader->Update();
    auto mb = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
    auto block = mb->GetBlock(0);
    auto dataSet = vtkDataSet::SafeDownCast(block);
    vtkCellData* cellData = dataSet->GetCellData();
    vtkDataArray* vofArray = cellData->GetArray("VOF");
    auto grid = vtkRectilinearGrid::SafeDownCast(dataSet);

    int dims[3];
    grid->GetDimensions(dims);  // dims are for points, so cells are dims-1
    int cellDims[3] = {dims[0] - 1, dims[1] - 1, dims[2] - 1};

    int downsample[3] = {2, 2, 2};
    int newCellDims[3] = {
        cellDims[0] / downsample[0],
        cellDims[1] / downsample[1],
        cellDims[2] / downsample[2]
    };
    int newDims[3] = {
        newCellDims[0] + 1,
        newCellDims[1] + 1,
        newCellDims[2] + 1
    };

    // Create new coordinate arrays
    auto newXCoords = vtkSmartPointer<vtkDoubleArray>::New();
    auto newYCoords = vtkSmartPointer<vtkDoubleArray>::New();
    auto newZCoords = vtkSmartPointer<vtkDoubleArray>::New();

    auto oldX = grid->GetXCoordinates();
    auto oldY = grid->GetYCoordinates();
    auto oldZ = grid->GetZCoordinates();

    for (int i = 0; i < newDims[0]; ++i)
        newXCoords->InsertNextValue(oldX->GetComponent(i * downsample[0], 0));
    for (int j = 0; j < newDims[1]; ++j)
        newYCoords->InsertNextValue(oldY->GetComponent(j * downsample[1], 0));
    for (int k = 0; k < newDims[2]; ++k)
        newZCoords->InsertNextValue(oldZ->GetComponent(k * downsample[2], 0));

    // Create new grid
    auto downsampledGrid = vtkSmartPointer<vtkRectilinearGrid>::New();
    downsampledGrid->SetDimensions(newDims);
    downsampledGrid->SetXCoordinates(newXCoords);
    downsampledGrid->SetYCoordinates(newYCoords);
    downsampledGrid->SetZCoordinates(newZCoords);

    // Create new VOF array
    auto newVOF = vtkSmartPointer<vtkDoubleArray>::New();
    newVOF->SetName("VOF");
    newVOF->SetNumberOfComponents(1);
    newVOF->SetNumberOfTuples(newCellDims[0] * newCellDims[1] * newCellDims[2]);

    int idx = 0;
    for (int k = 0; k < newCellDims[2]; ++k) {
        for (int j = 0; j < newCellDims[1]; ++j) {
            for (int i = 0; i < newCellDims[0]; ++i) {
                double sum = 0.0;
                int count = 0;
                for (int dk = 0; dk < downsample[2]; ++dk) {
                    for (int dj = 0; dj < downsample[1]; ++dj) {
                        for (int di = 0; di < downsample[0]; ++di) {
                            int oldI = i * downsample[0] + di;
                            int oldJ = j * downsample[1] + dj;
                            int oldK = k * downsample[2] + dk;

                            int oldIdx = oldK * (cellDims[0] * cellDims[1]) + oldJ * cellDims[0] + oldI;
                            sum += vofArray->GetComponent(oldIdx, 0);
                            ++count;
                        }
                    }
                }
                newVOF->SetValue(idx++, sum / count);
            }
        }
    }

    downsampledGrid->GetCellData()->AddArray(newVOF);

    // Convert point dims to cell dims (cells = points - 1)
    int nx = newCellDims[0];
    int ny = newCellDims[1];
    int nz = newCellDims[2];

    std::cout << "Downsampled grid cell dimensions:" << std::endl;
    std::cout << "  nx = " << nx << std::endl;
    std::cout << "  ny = " << ny << std::endl;
    std::cout << "  nz = " << nz << std::endl;

    // Create new cell data array for interface type
    auto interface_type = vtkSmartPointer<vtkIntArray>::New();
    interface_type->SetName("InterfaceType");
    interface_type->SetNumberOfComponents(1);
    interface_type->SetNumberOfTuples(grid->GetNumberOfCells());
    for (int i = 0; i < grid->GetNumberOfCells(); i++) {
        interface_type->SetValue(i, -1);
    }
    // Create certainty array
    auto certainty = vtkSmartPointer<vtkFloatArray>::New();
    certainty->SetName("Certainty");
    certainty->SetNumberOfComponents(1);
    certainty->SetNumberOfTuples(grid->GetNumberOfCells());
    for (int i = 0; i < grid->GetNumberOfCells(); i++) {
        certainty->SetValue(i, 0.0);
    }

    // Read PLIC file for viz
    auto plic_reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
    // plic_reader->SetCaseFileName("/home/quirin/jet/plic.case");
    plic_reader->SetCaseFileName("/home/quirin/mlcfd/Repositories/jet/plic.case");
    plic_reader->UpdateInformation();

    vtkInformation* plic_info = plic_reader->GetOutputInformation(0);
    int plic_nSteps = plic_info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    std::vector<double> plic_timeSteps(plic_nSteps);
    plic_info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), plic_timeSteps.data());
    plic_reader->SetTimeValue(plic_timeSteps[plic_nSteps - 1]);
    plic_reader->Update();
    auto plic_mb = vtkMultiBlockDataSet::SafeDownCast(plic_reader->GetOutput());
    auto plic_block_merger = vtkSmartPointer<vtkAppendFilter>::New();
    for (int i = 0; i < plic_mb->GetNumberOfBlocks(); i++) {
        auto plic_block = plic_mb->GetBlock(i);
        if (plic_block) { plic_block_merger->AddInputData(plic_block);}
    }
    plic_block_merger->Update();
    auto plic_grid = vtkUnstructuredGrid::SafeDownCast(plic_block_merger->GetOutput());

    // Create new cell data array for interface type on PLIC surface
    auto plic_interface_type = vtkSmartPointer<vtkIntArray>::New();
    plic_interface_type->SetName("InterfaceType");
    plic_interface_type->SetNumberOfComponents(1);
    plic_interface_type->SetNumberOfTuples(plic_grid->GetNumberOfCells());
    for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
        plic_interface_type->SetValue(i, -1);
    }
    // Create certainty array on PLIC surface
    auto plic_certainty = vtkSmartPointer<vtkFloatArray>::New();
    plic_certainty->SetName("Certainty");
    plic_certainty->SetNumberOfComponents(1);
    plic_certainty->SetNumberOfTuples(plic_grid->GetNumberOfCells());
    for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
        plic_certainty->SetValue(i, 0.0);
    }

    int stencil_size_reader = classifier.getStencilSize(); // now configurable, must match training
    int half = stencil_size_reader / 2;

    int no_filled_cells = 0;
    int no_paraboloids = 0;
    int no_cylinders = 0;
    int no_spheres = 0;
    int no_sheets = 0;

    // Loop through interior cells (excluding half-cell boundary)
    for (int i = half; i < nx - half; ++i) {
        for (int j = half; j < ny - half; ++j) {
            for (int k = half; k < nz - half; ++k) {

                vtkIdType centerCellId = i + j * nx + k * nx * ny;

                // Skip if no interface
                double epsilon = 1e-10;
                double vof_center = newVOF->GetComponent(centerCellId, 0);
                if (vof_center <= epsilon || vof_center >= 1.0 - epsilon) {
                    continue;
                }

                no_filled_cells++;

                // Allocate vector for stencil
                std::vector<double> vfrac(stencil_size_reader * stencil_size_reader * stencil_size_reader, 0.0);

                // Fill stencil neighborhood
                for (int di = -half; di <= half; ++di) {
                    for (int dj = -half; dj <= half; ++dj) {
                        for (int dk = -half; dk <= half; ++dk) {
                            int ni = i + di;
                            int nj = j + dj;
                            int nk = k + dk;

                            vtkIdType cellId = ni + nj * nx + nk * nx * ny;
                            if (cellId >= newVOF->GetNumberOfTuples()) {
                                std::cerr << "Invalid cellId: " << cellId << std::endl;
                                continue;
                            }

                            int idx = (di + half) * stencil_size_reader * stencil_size_reader
                                    + (dj + half) * stencil_size_reader
                                    + (dk + half);
                            vfrac[idx] = newVOF->GetComponent(cellId, 0);
                        }
                    }
                }

                // Classify directly on vfrac
                std::vector<float> out_probs;
                int predicted_class = classifier.classify(vfrac, &out_probs);
                float max_prob = *std::max_element(out_probs.begin(), out_probs.end());
                certainty->SetValue(centerCellId, max_prob);

                switch (predicted_class) {
                    case 0:
                        no_paraboloids++;
                        interface_type->SetValue(centerCellId, 0);
                        break;
                    case 1:
                        no_cylinders++;
                        interface_type->SetValue(centerCellId, 1);
                        break;
                    case 2:
                        no_spheres++;
                        interface_type->SetValue(centerCellId, 2);
                        break;
                    case 3:
                        no_sheets++;
                        interface_type->SetValue(centerCellId, 3);
                        break;
                    default:
                        std::cerr << "Warning: unknown predicted_class = " << predicted_class << std::endl;
                        break;
                }
            }
        }
    }

  

    std::cout << "\n=== Classification Summary ===" << std::endl;
    std::cout << "Total cells with interface:   " << no_filled_cells << std::endl;
    std::cout << "Paraboloids:          " << no_paraboloids << std::endl;
    std::cout << "Cylinders:            " << no_cylinders << std::endl;
    std::cout << "Spheres:              " << no_spheres << std::endl;
    std::cout << "Sheets:               " << no_sheets << std::endl;

    // Write Grid file
    grid->GetCellData()->AddArray(interface_type);
    grid->GetCellData()->AddArray(certainty);
    auto grid_writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
    grid_writer->SetFileName("grid.vtr");
    grid_writer->SetInputData(downsampledGrid);
    grid_writer->Write();

    // Convert interface_type to PLIC array
    auto locator = vtkSmartPointer<vtkCellLocator>::New();
    locator->SetDataSet(downsampledGrid);
    locator->BuildLocator();
    auto cell_center_filter = vtkSmartPointer<vtkCellCenters>::New();
    cell_center_filter->SetInputData(plic_grid);
    cell_center_filter->Update();
    auto cell_centers = cell_center_filter->GetOutput()->GetPoints();

    for (int i = 0; i < plic_grid->GetNumberOfCells(); i++) {
        auto type = interface_type->GetValue(locator->FindCell(cell_centers->GetPoint(i)));
        plic_interface_type->SetValue(i, type);
        auto cert = certainty->GetValue(locator->FindCell(cell_centers->GetPoint(i)));
        plic_certainty->SetValue(i, cert);
    }

    // Write PLIC file
    plic_grid->GetCellData()->AddArray(plic_interface_type);
    plic_grid->GetCellData()->AddArray(plic_certainty);
    auto plic_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    plic_writer->SetFileName("plic.vtu");
    plic_writer->SetInputData(plic_grid);
    plic_writer->Write();
}
} // namespace IRL