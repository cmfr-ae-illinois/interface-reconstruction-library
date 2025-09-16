// #include <torch/torch.h>
#include <iostream>
#include <vector>
#include <string>

#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/data_set.h"
#include "irl/ml_classification/net.h"
#include "irl/ml_classification/trainer.h"

#include <random> //REMOVE, was just for testing

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
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkAppendFilter.h>
#include <vtkCellLocator.h>
#include <vtkPolygon.h>
#include <vtkCellCenters.h>

int main (int argc, char* argv[]) {

  // To test if Torch works
  torch::Tensor tensor = torch::eye(3);
  std::cout << tensor << std::endl;

  int stencil_size = 3;
  int hidden_size1 = 128;
  int hidden_size2 = 64;
  int hidden_size3 = 32;
  int output_size = 4; // 0 = plane (currently not used), 1 = paraboloid, 2 = cylinder, 3 = 
  int batch_size = 64;
  double learning_rate = 0.001; //was 0.01 for SGD optimizer
  int no_batches = 2048; // Should be divisible by batch size
  int epochs = 20;
  std::vector<double> lossVector;

  // Declare vectors of states and labels
  std::vector<std::vector<double>> statesV;
  std::vector<int> labelsV;

  IRL::Data_gen data_gen;
  data_gen.generate_Data(&statesV, &labelsV, no_batches*batch_size, stencil_size, output_size);

  // Split the data into training, test, and validation
  int total_samples = statesV.size();
  int train_end = total_samples * 0.7;
  int test_end = total_samples * 0.85; // Next 15%

  std::vector<std::vector<double>> train_states(statesV.begin(), statesV.begin() + train_end);
  std::vector<int> train_labels(labelsV.begin(), labelsV.begin() + train_end);

  std::vector<std::vector<double>> test_states(statesV.begin() + train_end, statesV.begin() + test_end);
  std::vector<int> test_labels(labelsV.begin() + train_end, labelsV.begin() + test_end);

  std::vector<std::vector<double>> val_states(statesV.begin() + test_end, statesV.end());
  std::vector<int> val_labels(labelsV.begin() + test_end, labelsV.end());

  IRL::MyDataset train_dataset(&train_states, &train_labels);
  IRL::MyDataset test_dataset(&test_states, &test_labels);
  IRL::MyDataset val_dataset(&val_states, &val_labels);

  auto train_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
      std::move(train_dataset).map(torch::data::transforms::Stack<>()), batch_size);

  auto test_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
      std::move(test_dataset).map(torch::data::transforms::Stack<>()), batch_size);

  auto val_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
      std::move(val_dataset).map(torch::data::transforms::Stack<>()), batch_size);

  IRL::Net net(stencil_size*stencil_size*stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);

  // Create Trainer and run
  torch::optim::Adam optimizer(net.parameters(), torch::optim::AdamOptions(learning_rate));
  IRL::Trainer<IRL::Net> trainer(net, optimizer, train_loader.get(), test_loader.get(), val_loader.get(), epochs);
  trainer.train();

  // vtk reader

  auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
  reader->SetCaseFileName("/home/quirin/mlcfd/Repositories/jet/nga.case");
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

  const int stencil_size_reader = 3;
  int no_filled_cells = 0;
  int no_paraboloids = 0;
  int no_cylinders = 0;
  int no_spheres = 0;
  int no_sheets = 0;

  // Loop through interior cells (excluding 1-cell boundary)
  for (int i = 1; i < nx - 1; ++i) {
      for (int j = 1; j < ny - 1; ++j) {
          for (int k = 1; k < nz - 1; ++k) {

              vtkIdType centerCellId = i + j * nx + k * nx * ny;

              // Skip this cell if there is no surface
              double epsilon = 1e-10;
              double vof_center = newVOF->GetComponent(centerCellId, 0);
              if (vof_center <= epsilon || vof_center >= 1.0 - epsilon) {
                  continue;
              }

              no_filled_cells ++;

              // 3D array to store VOF values in the 3x3x3 neighborhood
              double vfrac[stencil_size_reader][stencil_size_reader][stencil_size_reader];

              // Fill the 3x3x3 neighborhood
              for (int di = -1; di <= 1; ++di) {
                  for (int dj = -1; dj <= 1; ++dj) {
                      for (int dk = -1; dk <= 1; ++dk) {
                          int ni = i + di;
                          int nj = j + dj;
                          int nk = k + dk;

                          vtkIdType cellId = ni + nj * nx + nk * nx * ny;
                          if (cellId >= newVOF->GetNumberOfTuples()) {
                              std::cerr << "Invalid cellId: " << cellId << std::endl;
                              continue;
                          }

                          vfrac[di + 1][dj + 1][dk + 1] = newVOF->GetComponent(cellId, 0);
                      }
                  }
              }

              // Flatten the 3D vector into a 1D vector
              std::vector<double> flattened_vfrac;
              for (int ii = 0; ii < stencil_size_reader; ++ii) {
                  for (int jj = 0; jj < stencil_size_reader; ++jj) {
                      for (int kk = 0; kk < stencil_size_reader; ++kk) {
                          flattened_vfrac.push_back(vfrac[ii][jj][kk]);
                      }
                  }
              }

              // Classify flattened vfrac vector
              int predicted_class = net.forward(torch::tensor(flattened_vfrac)).argmax().item<int>();
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
  }

  // Write PLIC file
  plic_grid->GetCellData()->AddArray(plic_interface_type);
  auto plic_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  plic_writer->SetFileName("plic.vtu");
  plic_writer->SetInputData(plic_grid);
  plic_writer->Write();


  return 0;
}