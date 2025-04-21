#include <torch/torch.h>

namespace IRL {
    struct Net : torch::nn::Module {
        // Define the four linear layers: three hidden layers and one output layer
        torch::nn::Linear fc1, fc2, fc3, fc4;

        // Constructor to initialize the layers
        Net(int input_size, int hidden_size1, int hidden_size2, int hidden_size3, int output_size)
            : fc1(input_size, hidden_size1),  // First hidden layer
            fc2(hidden_size1, hidden_size2), // Second hidden layer
            fc3(hidden_size2, hidden_size3), // Third hidden layer
            fc4(hidden_size3, output_size)   // Output layer
        {
            // Register the layers as sub-modules
            register_module("fc1", fc1);
            register_module("fc2", fc2);
            register_module("fc3", fc3);
            register_module("fc4", fc4);
        }

        // Forward pass definition with ReLU activation after each hidden layer
        torch::Tensor forward(torch::Tensor x) {
            x = torch::relu(fc1(x));  // Apply first hidden layer + ReLU activation
            x = torch::relu(fc2(x));  // Apply second hidden layer + ReLU activation
            x = torch::relu(fc3(x));  // Apply third hidden layer + ReLU activation
            x = fc4(x);               // Apply output layer
            return x;
        }
    };
}