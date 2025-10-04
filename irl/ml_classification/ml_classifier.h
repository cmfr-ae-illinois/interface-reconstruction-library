#pragma once
#include "classifier.h"
#include "data_gen.h"
#include "data_set.h"
#include "trainer.h"
#include "net.h"

#include <torch/torch.h>
#include <iostream>

namespace IRL {

class MLClassifier : public Classifier {
private:
    int hidden1, hidden2, hidden3, output_size;
    int batch_size, epochs;
    double learning_rate;

    Net net;
    std::unique_ptr<torch::optim::Optimizer> optimizer;

    // Store training data
    std::vector<std::vector<double>> statesV;
    std::vector<int> labelsV;

public:
    MLClassifier(int stencil, int input_size, int h1, int h2, int h3, int out_size)
        : Classifier(stencil),
          hidden1(h1), hidden2(h2), hidden3(h3), output_size(out_size),
          net(input_size, h1, h2, h3, out_size) {
    }

    void generateDataset(int no_batches, int batch_size=64) {
        IRL::Data_gen data_gen;
        data_gen.generate_Data(&statesV, &labelsV, no_batches * batch_size,
                               stencil_size, output_size);
    }

    void trainModel(int epochs=20, double learning_rate=0.001, int batch_size=64) {
        optimizer = std::make_unique<torch::optim::Adam>(net.parameters(), torch::optim::AdamOptions(learning_rate));

        int total_samples = statesV.size();
        int train_end = total_samples * 0.7;
        int test_end = total_samples * 0.85;

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

        using Loader = std::remove_reference_t<decltype(*train_loader)>;
        IRL::Trainer<IRL::Net, Loader> trainer(net, *optimizer, train_loader.get(), test_loader.get(), val_loader.get(), epochs);
        trainer.train();
    }

    int classify(const std::vector<double>& flattened_state) override {
        
        torch::Tensor input = torch::from_blob(
            const_cast<double*>(flattened_state.data()),
            {1, (long)(flattened_state.size())},
            torch::kFloat32);

        torch::Tensor output = net.forward(input);
        int predicted_class = output.argmax(1).item<int>();
        
        //int predicted_class = net.forward(torch::tensor(flattened_state)).argmax().item<int>();
        return predicted_class;
    }
};

} // namespace IRL