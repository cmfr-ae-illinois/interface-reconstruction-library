// #include <torch/torch.h>
#include <iostream>
#include <vector>
#include <string>

#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/data_set.h"
#include "irl/ml_classification/net.h"

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

template <class CellData, class CoordData>
void WriteField(const CoordData& coordsx, const CoordData& coordsy, const CoordData& coordsz,
                const CellData& field, const std::string& field_name) {
  const auto file_name = "./" + field_name + ".vtr"; 

  double start = 0.0;
  double end = 1.0;

  const int ncellsx = coordsx.size() - 1;
  const int ncellsy = coordsy.size() - 1;
  const int ncellsz = coordsz.size() - 1;

  FILE* file;
  file = fopen(file_name.c_str(), "w");

  fprintf(file, "<VTKFile type=\"RectilinearGrid\">\n");
  fprintf(file, "<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0,
          ncellsx, 0, ncellsy, 0, ncellsz);
  fprintf(file, "<Piece Extent=\"%d %d %d %d %d %d\">\n", 0, ncellsx, 0, ncellsy,
          0, ncellsz);

  fprintf(file, "<Coordinates>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < ncellsx + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(coordsx[i]));
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < ncellsy + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(coordsy[i]));
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < ncellsz + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(coordsz[i]));
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file, "</Coordinates>\n");

  fprintf(file, "<PointData>\n</PointData>\n");

  fprintf(file, "<CellData Scalars=\"");
  fprintf(file, "%s ", field_name.c_str());
  fprintf(file, "\" >\n");
  fprintf(file,
          "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n",
          field_name.c_str());
  for (int k = 0; k < ncellsz; ++k) {
    for (int j = 0; j < ncellsy; ++j) {
      for (int i = 0; i < ncellsx; ++i) {
        fprintf(file, "%15.8E ", static_cast<double>(field[i][j][k]));
      }
    }
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n</RectilinearGrid>\n</VTKFile>\n");
  fclose(file);
}

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
  int no_batches = 1024; // Should be divisible by batch size
  int epochs = 20;
  std::vector<double> lossVector;

  // Declare vectors of states and labels
  std::vector<std::vector<double>> statesV;
  std::vector<int> labelsV;

  IRL::Data_gen data_gen;
  data_gen.generate_Data(&statesV, &labelsV, no_batches*batch_size, stencil_size, output_size); //if include planes = true, dont forget to adjust output size.

  //output:
  //std::cout << "States:" << std::endl;
  // Iterate over the outer vector
  //for (const auto& innerVec : statesV) {
      // Iterate over the inner vector
      //for (double val : innerVec) {
          //std::cout << val << " ";  // Print each element of the inner vector
      //}
      //std::cout << std::endl;  // Print a newline after each inner vector
  //}
  //std::cout << "Labels:" << std::endl;
  //for (const int& num : labelsV) {
      //std::cout << num << " ";
  //}

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

  
  
  
  //IRL::MyDataset dataset(&statesV, &labelsV);
  // Accessing the first sample
  //auto example = dataset.get(0);  // Get the first sample
  // Print the state and label for the first sample
  //std::cout << "State: " << example.data << std::endl;
  //std::cout << "Label: " << example.target << std::endl;

  IRL::Net net(stencil_size*stencil_size*stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);

  //auto data_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(

  //std::move(dataset).map(torch::data::transforms::Stack<>()), batch_size);

  //torch::optim::SGD optimizer(net.parameters(), learning_rate);
  torch::optim::Adam optimizer(net.parameters(), torch::optim::AdamOptions(learning_rate));

  for (int epoch=1; epoch<=epochs; ++epoch) {
    int batch_index = 0;
    double total_loss = 0.0;
    int correct_predictions = 0;  // To count the number of correct predictions
    int total_samples = 0;        // To count the total number of samples processed

    // Iterate data loader to yield batches from the dataset

    for (auto& batch: *train_loader) {
      // Reset gradients
			optimizer.zero_grad();

      // Time measurement before
      //auto start_time = std::chrono::high_resolution_clock::now();

      // Execute model
      torch::Tensor prediction = net.forward(batch.data);

      // Time measurement after
      //auto current_time = std::chrono::high_resolution_clock::now();
      //std::cout << "Program has been running for " << std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count() << " milliseconds" << std::endl;
      
      // Compute loss value
      torch::Tensor loss = torch::cross_entropy_loss(prediction, batch.target);

      // To use negative log likelyhood loss, implement softmax layer. Alternatively above use cross entropy loss, has softmax built in
      //torch::Tensor loss = torch::nll_loss(prediction, batch.target);

      // Backward pass: Compute gradients and update the parameters
			loss.backward();
			optimizer.step();
      total_loss += loss.item<double>();

      // Calculate accuracy for the current batch
      // Get predicted class by taking the index of the max logit (log probability)
      auto predicted_classes = prediction.argmax(1); // For multi-class classification
      auto correct = predicted_classes.eq(batch.target); // Check if predictions match the true labels
      correct_predictions += correct.sum().item().toInt();  // Sum of correct predictions in the batch
      total_samples += batch.target.size(0);  // Total number of samples in the batch
    }
    
    lossVector.push_back(total_loss);

    //double accuracy = static_cast<double>(correct_predictions) / total_samples;
    //std::cout << "Epoch " << epoch << " - Loss: " << total_loss 
    //          << " - Accuracy: " << accuracy * 100 << "%" << std::endl;
    //std::cout << "Loss:" << std::endl;

    // Training metrics
    std::cout << "Epoch [" << epoch << "/" << epochs << "] , Accuracy: "
              << static_cast<double>(correct_predictions) / total_samples << std::endl;

    // --- Test Accuracy Evaluation ---
    int test_correct = 0;
    int test_total = 0;

    for (auto& batch : *test_loader) {
      auto prediction = net.forward(batch.data);
      auto predicted = prediction.argmax(1);
      test_correct += predicted.eq(batch.target).sum().item<int>();
      test_total += batch.target.size(0);
    }

    double test_accuracy = static_cast<double>(test_correct) / test_total;
    std::cout << "Test Accuracy: " << test_accuracy << std::endl;
  }

  std::cout << "Loss vector: " << std::endl;
  for (const int& num : lossVector) {
    std::cout << num << std::endl;
  }

  // --- Final Validation Accuracy ---
  int val_correct = 0;
  int val_total = 0;

  for (auto& batch : *val_loader) {
      auto prediction = net.forward(batch.data);
      auto predicted = prediction.argmax(1);
      val_correct += predicted.eq(batch.target).sum().item<int>();
      val_total += batch.target.size(0);
  }

  double val_accuracy = static_cast<double>(val_correct) / val_total;
  std::cout << "Final Validation Accuracy: " << val_accuracy << std::endl;



  // vtk reader

  auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
//   reader->SetCaseFileName("/home/quirin/mlcfd/Repositories/jet/nga.case");
  reader->SetCaseFileName("/home/quirin/jet/nga.case");
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
  plic_reader->SetCaseFileName("/home/quirin/jet/plic.case");
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

              // Skip this cell if VOF = 0
              double vof_center = vofArray->GetComponent(centerCellId, 0);
              if (vof_center == 0) {
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
  std::cout << "Total filled cells:   " << no_filled_cells << std::endl;
  std::cout << "Paraboloids:          " << no_paraboloids << std::endl;
  std::cout << "Cylinders:            " << no_cylinders << std::endl;
  std::cout << "Spheres:              " << no_spheres << std::endl;
  std::cout << "Sheets:               " << no_sheets << std::endl;

  // Write Grid file
  grid->GetCellData()->AddArray(interface_type);
  auto grid_writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
  grid_writer->SetFileName("grid.vtr");
  grid_writer->SetInputData(grid);
  grid_writer->Write();

  // Convert interface_type to PLIC array
  auto locator = vtkSmartPointer<vtkCellLocator>::New();
  locator->SetDataSet(grid);
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