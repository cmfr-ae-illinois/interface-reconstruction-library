#include <torch/torch.h>
#include <iostream>
#include <vector>
#include <string>

#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/data_set.h"
#include "irl/ml_classification/net.h"

#include <random> //REMOVE, was just for testing

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

  // make a prediction based on a 3x3 vector
  std::vector<double> fractions1 = data_gen.generate_Cylinder(3);
  // Whithin generate_Cylinder, it is already compressed into a 1 dimensional vector. Go into that function to find out how

  int predicted_class1 = net.forward(torch::tensor(fractions1)).argmax().item<int>();
  // Current correct classes: 0=paraboloid, 1=cylinder, 2=sphere, 3=sheet

  std::cout << "Predicted class: " << predicted_class1 << std::endl;

    // make a prediction based on a 3x3 vector
  std::vector<double> fractions2 = data_gen.generate_Cylinder(3);
  // Whithin generate_Cylinder, it is already compressed into a 1 dimensional vector. Go into that function to find out how

  int predicted_class2 = net.forward(torch::tensor(fractions1)).argmax().item<int>();
  // Current correct classes: 0=paraboloid, 1=cylinder, 2=sphere, 3=sheet

  std::cout << "Predicted class: " << predicted_class2 << std::endl;

  return 0;
}