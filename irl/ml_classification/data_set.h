#include <torch/torch.h>
#include <vector>

namespace IRL {
    torch::Tensor read_states(std::vector<std::vector<float>>* statesV) {
        // Flatten the 2D vector into a 1D vector
        std::vector<float> flattened;
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
            explicit MyDataset(std::vector<std::vector<float>>* statesV, std::vector<int>* labelsV) 
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