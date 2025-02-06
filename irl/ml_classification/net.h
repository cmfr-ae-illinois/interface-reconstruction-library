#include <torch/torch.h>

namespace IRL {
    struct Net : torch::nn::Module {
        // Define the two linear layers
        torch::nn::Linear fc1, fc2;

        // Constructor to initialize the layers
        Net(int input_size, int hidden_size, int output_size)
            : fc1(input_size, hidden_size),  // First layer
              fc2(hidden_size, output_size)   // Second layer
        {
            // Register the layers as sub-modules
            register_module("fc1", fc1);
            register_module("fc2", fc2);
        }

        // Forward pass definition with ReLU activation
        torch::Tensor forward(torch::Tensor x) {
            x = torch::relu(fc1(x));  // Apply first linear layer + ReLU activation
            x = fc2(x);               // Apply second linear layer (output)
            return x;
        }
    };
}