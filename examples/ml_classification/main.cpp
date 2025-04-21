#include <torch/torch.h>
#include <iostream>
#include <vector>
#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/data_set.h"
#include "irl/ml_classification/net.h"

#include <random> //REMOVE, was just for testing

int main () {
    // To test if Torch works
    torch::Tensor tensor = torch::eye(3);
    std::cout << tensor << std::endl;

    int stencil_size = 3;
    int hidden_size1 = 100;
    int hidden_size2 = 100;
    int hidden_size3 = 100;
    int output_size = 4; // Two numbers 0 = plane, 1 = paraboloid, 2 = cylinder
    int batch_size = 64;
    double learning_rate = 0.01;
    int no_batches = 1024; // Should be divisible by batch size?
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


    IRL::MyDataset dataset(&statesV, &labelsV);
    // Accessing the first sample
    auto example = dataset.get(0);  // Get the first sample

    // Print the state and label for the first sample
    std::cout << "State: " << example.data << std::endl;
    std::cout << "Label: " << example.target << std::endl;

    IRL::Net net(stencil_size*stencil_size*stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);

    auto data_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
        std::move(dataset).map(torch::data::transforms::Stack<>()), 
        batch_size);


    torch::optim::SGD optimizer(net.parameters(), learning_rate);
    
    for (int epoch=1; epoch<=epochs; ++epoch) {
        int batch_index = 0;
        double total_loss = 0.0;
        int correct_predictions = 0;  // To count the number of correct predictions
        int total_samples = 0;        // To count the total number of samples processed

        // Iterate data loader to yield batches from the dataset

        for (auto& batch: *data_loader) {
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
            correct_predictions += correct.sum().item<int>();  // Sum of correct predictions in the batch
            total_samples += batch.target.size(0);  // Total number of samples in the batch
        }
        lossVector.push_back(total_loss);

        double accuracy = static_cast<double>(correct_predictions) / total_samples;
        std::cout << "Epoch " << epoch << " - Loss: " << total_loss 
                  << " - Accuracy: " << accuracy * 100 << "%" << std::endl;
    }

    
    std::cout << "Loss:" << std::endl;
    for (const int& num : lossVector) {
        std::cout << num << std::endl;
    }

    return 0;
}