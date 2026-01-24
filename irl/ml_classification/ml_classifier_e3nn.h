#pragma once
#include "ml_classifier.h"
#include "e3nn.h"   // defines e3nnImpl and TORCH_MODULE(e3nn)

#include <torch/torch.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>

namespace IRL {

class MLClassifier_E3NN : public MLClassifier {
public:
    // Equivariant network holder (ModuleHolder around e3nnImpl)
    e3nn equivariant_net;

    // Track how long training takes (in seconds)
    double training_time = 0.0;
    double final_validation_accuracy = 0.0;

    MLClassifier_E3NN(int stencil = 3, int h1 = 128, int h2 = 64, int h3 = 32, int out = 4)
        : MLClassifier(
              /*stencil=*/stencil,
              /*input=*/stencil * stencil * stencil * 4, // vof + 3 barycenter coords per cell
              /*h1=*/h1, /*h2=*/h2, /*h3=*/h3, /*out=*/out),
          equivariant_net(stencil, h1, h2, h3, out)
    {
        std::cout << "Initialized E(3)-equivariant classifier with stencil size "
                  << stencil << std::endl;
    }

    int classify(const std::vector<float>& flattened_state,
                 std::vector<float>* out_probs = nullptr) override
    {
        torch::NoGradGuard no_grad;

        auto input = torch::tensor(flattened_state,
                                   torch::TensorOptions().dtype(torch::kFloat32))
                         .unsqueeze(0); // [1, input_size]

        // Send input to same device as network
        auto params = equivariant_net->parameters();
        if (!params.empty()) {
            auto device = params.front().device();
            auto dtype  = params.front().dtype();
            input = input.to(device, dtype);
        }

        torch::Tensor logits = equivariant_net->forward(input);
        torch::Tensor probs  = torch::softmax(logits, 1).to(torch::kCPU).squeeze(0);

        if (out_probs) {
            out_probs->resize(probs.size(0));
            for (int i = 0; i < probs.size(0); ++i)
                (*out_probs)[i] = probs[i].item<float>();
        }

        return probs.argmax().item<int>();
    }

    void trainModel() {
        optimizer = std::make_unique<torch::optim::Adam>(
            equivariant_net->parameters(), torch::optim::AdamOptions(learning_rate));

        const int total_samples = static_cast<int>(statesV.size());
        const int train_end = static_cast<int>(total_samples * 0.7);
        const int test_end  = static_cast<int>(total_samples * 0.85);

        std::vector<std::vector<float>> train_states(statesV.begin(), statesV.begin() + train_end);
        std::vector<int>                 train_labels(labelsV.begin(), labelsV.begin() + train_end);

        std::vector<std::vector<float>> test_states(statesV.begin() + train_end, statesV.begin() + test_end);
        std::vector<int>                 test_labels(labelsV.begin() + train_end, labelsV.begin() + test_end);

        std::vector<std::vector<float>> val_states(statesV.begin() + test_end, statesV.end());
        std::vector<int>                 val_labels(labelsV.begin() + test_end, labelsV.end());

        IRL::MyDataset train_dataset(&train_states, &train_labels);
        IRL::MyDataset test_dataset (&test_states,  &test_labels);
        IRL::MyDataset val_dataset  (&val_states,   &val_labels);

        auto train_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
            std::move(train_dataset).map(torch::data::transforms::Stack<>()), batch_size);
        auto test_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
            std::move(test_dataset).map(torch::data::transforms::Stack<>()), batch_size);
        auto val_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
            std::move(val_dataset).map(torch::data::transforms::Stack<>()), batch_size);

        using Loader = std::remove_reference_t<decltype(*train_loader)>;

        IRL::Trainer<IRL::e3nnImpl, Loader> trainer(
            *equivariant_net, *optimizer,
            train_loader.get(), test_loader.get(), val_loader.get(), epochs);

        // Measure training time
        using namespace std::chrono;
        auto start_time = high_resolution_clock::now();

        final_validation_accuracy = trainer.train();

        auto end_time = high_resolution_clock::now();
        training_time = duration_cast<duration<double>>(end_time - start_time).count();

        std::cout << "⏳ Training completed in " << training_time << " seconds." << std::endl;
    }

    void saveModel(const std::string& path, bool save_optimizer_state = false) {
        namespace fs = std::filesystem;
        if (!path.empty()) {
            fs::create_directories(fs::path(path).parent_path());
        }

        torch::serialize::OutputArchive archive;
        equivariant_net->save(archive);
        archive.save_to(path + "ml_model_e3nn.pt");
        std::cout << "💾 Saved E3NN model to " << path << std::endl;

        if (save_optimizer_state && optimizer) {
            std::string opt_path = path + ".opt";
            torch::save(*optimizer, opt_path);
            std::cout << "💾 Saved optimizer state to " << opt_path << std::endl;
        }

        std::string param_file = path + "_Parameters.txt";
        std::ofstream param_out(param_file);
        if (!param_out)
            throw std::runtime_error("Failed to open " + param_file);

        param_out << "[Network]\n";
        param_out << "stencil_size " << stencil_size << "\n";
        param_out << "hidden_size1 " << hidden_size1 << "\n";
        param_out << "hidden_size2 " << hidden_size2 << "\n";
        param_out << "hidden_size3 " << hidden_size3 << "\n";
        param_out << "output_size " << output_size << "\n\n";

        param_out << "[Training]\n";
        param_out << "learning_rate " << learning_rate << "\n";
        param_out << "batch_size " << batch_size << "\n";
        param_out << "epochs " << epochs << "\n";
        param_out << "training_time_seconds " << training_time << "\n";
        param_out << "final_validation_accuracy " << final_validation_accuracy << "\n";

        param_out.close();
        std::cout << "📝 Saved model parameters (including training time) to " << param_file << std::endl;
    }

    void loadModel(const std::string& path, bool load_optimizer_state = false) {
        torch::serialize::InputArchive archive;
        archive.load_from(path);
        equivariant_net->load(archive);
        std::cout << "📂 Loaded trained E3NN model from " << path << std::endl;

        if (load_optimizer_state) {
            if (!optimizer)
                optimizer = std::make_unique<torch::optim::Adam>(
                    equivariant_net->parameters(), torch::optim::AdamOptions(learning_rate));
            std::string opt_path = path + ".opt";
            torch::load(*optimizer, opt_path);
            std::cout << "📂 Loaded optimizer state from " << opt_path << std::endl;
        }
    }
};

} // namespace IRL
