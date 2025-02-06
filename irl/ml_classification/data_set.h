#include <torch/torch.h>
#include <vector>

namespace IRL {
    // OLD, REMOVE SOON Read the data and directly store it as tensor
    torch::Tensor read_data(int states_or_labels) {
        if (states_or_labels == 0) {
            // Define multiple vectors
            std::vector<double> v1 = {1.01, 2.33};
            std::vector<double> v2 = {2.01, 5.536};
            std::vector<double> v3 = {1.0001, 8.0};

            // Create a 2D vector (vector of vectors)
            std::vector<std::vector<double>> v = {v1, v2, v3};

            std::cout << v.size() << std::endl;
            std::cout << v[0].size() << std::endl;

            // Flatten the 2D vector into a 1D vector
            std::vector<double> flattened;
            for (const auto& vec : v) {
                flattened.insert(flattened.end(), vec.begin(), vec.end());
            }

            // Convert the flattened vector to a 1D tensor
            torch::Tensor tStates = torch::tensor(flattened, torch::kDouble);

            // Now reshape the tensor to match the shape of v (3 rows, 2 columns)
            tStates = tStates.view({v.size(), v[0].size()});
            return tStates;

        } else if (states_or_labels == 1) {
            std::vector<int> vStates = {0, 1, 0};
            // Convert the std::vector<int> to a torch::Tensor
            torch::Tensor tLabels = torch::tensor(vStates, torch::kInt32);
            return tLabels;        

        } else {
            std::cout << "Didnt properly specify if state (0) or label (1)" << std::endl; 
            torch::Tensor error_tensor = torch::eye(1);
            return error_tensor;
        }
    };

    torch::Tensor read_states(std::vector<std::vector<double>>* statesV) {
        // Flatten the 2D vector into a 1D vector
        std::vector<double> flattened;
        for (const auto& vec : *statesV) {
            flattened.insert(flattened.end(), vec.begin(), vec.end());
        }

        // Convert the flattened vector to a 1D tensor, use float instead of double, otherwise: mat1 and mat2 must have the same dtype, but got Double and Float
        torch::Tensor tStates = torch::tensor(flattened, torch::kFloat);

        // Now reshape the tensor to match the shape of v (3 rows, 2 columns)
        tStates = tStates.view({statesV->size(), statesV->at(0).size()});
        return tStates;
    };

    torch::Tensor read_labels(std::vector<int>* labelsV) {
        torch::Tensor tLabels = torch::tensor(*labelsV, torch::kLong); // Use LongTensor as labels, it is equivalent to int. Otherwise no work 
        return tLabels;    
    };

    class MyDataset : public torch::data::Dataset<MyDataset>
    {
        private:
            torch::Tensor states, labels;

        public:
            explicit MyDataset(std::vector<std::vector<double>>* statesV, std::vector<int>* labelsV) 
                : states(read_states(statesV)),
                labels(read_labels(labelsV)) {   };

            //torch::data::Example<> get(size_t index) override;

            // Implement the `get()` function to return an example
            torch::data::Example<> get(size_t index) override {
                return {states[index], labels[index]}; // Return state and label at the given index
            }

            // Implement the `size()` function to return the number of samples
            torch::optional<size_t> size() const override {
                return states.size(0);  // Number of rows in the `states` tensor (number of samples)
            }
    };

    //torch::data::Example<> MyDataset::get(size_t index) {
        //return {states[index], labels[index]};
    //} 
}