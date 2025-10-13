#pragma once
#include "classifier.h"
#include "data_gen.h"
#include "data_set.h"
#include "trainer.h"
#include "net.h"

#include <torch/torch.h>
#include <iostream>
#include <fstream>
#include <filesystem>

namespace IRL {

class MLClassifier : public Classifier {
private:
    int hidden1, hidden2, hidden3, output_size;
    int batch_size, epochs;
    double learning_rate;

    Net net;
    std::unique_ptr<torch::optim::Optimizer> optimizer;
    std::string dataset_path = "";

    // Store training data
    std::vector<std::vector<double>> statesV;
    std::vector<int> labelsV;

public:
    MLClassifier(int stencil, int input_size, int h1, int h2, int h3, int out_size)
        : Classifier(stencil),
          hidden1(h1), hidden2(h2), hidden3(h3), output_size(out_size),
          net(input_size, h1, h2, h3, out_size) {}

    void generateDataset(int no_batches, int batch_size = 64,
                         int include_Moments = 0,
                         double paraboloid_coeff_stddev = 0.1,
                         double sheet_coeff_stddev = 0.1,
                         double max_sheet_thickness = 0.5,
                         double sheet_thickness_stddev = 0.0,
                         double max_cylinder_radius = 0.5,
                         double cylinder_radius_stddev = 0.0,
                         double max_sphere_radius = 0.5,
                         double sphere_radius_stddev = 0.0)
    {
        IRL::Data_gen data_gen;
        data_gen.generate_Data(&statesV, &labelsV,
                               no_batches * batch_size,
                               stencil_size, output_size,
                               include_Moments,
                               paraboloid_coeff_stddev,
                               sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev,
                               max_cylinder_radius, cylinder_radius_stddev,
                               max_sphere_radius, sphere_radius_stddev);

        // Save dataset immediately
        saveDataset("data",
                    no_batches, batch_size, include_Moments,
                    paraboloid_coeff_stddev, sheet_coeff_stddev,
                    max_sheet_thickness, sheet_thickness_stddev,
                    max_cylinder_radius, cylinder_radius_stddev,
                    max_sphere_radius, sphere_radius_stddev);
    }

    void saveDataset(const std::string& dir_path,
                     int no_batches, int batch_size, int include_Moments,
                     double paraboloid_coeff_stddev,
                     double sheet_coeff_stddev,
                     double max_sheet_thickness, double sheet_thickness_stddev,
                     double max_cylinder_radius, double cylinder_radius_stddev,
                     double max_sphere_radius, double sphere_radius_stddev) const
    {
        namespace fs = std::filesystem;
        fs::create_directories(dir_path);

        // --- Save data binary ---
        std::string data_file = dir_path + "/data.bin";
        std::ofstream data_out(data_file, std::ios::binary);
        if (!data_out) {
            throw std::runtime_error("Failed to open " + data_file);
        }

        size_t num_samples = statesV.size();
        size_t input_size = (num_samples > 0) ? statesV[0].size() : 0;

        data_out.write(reinterpret_cast<const char*>(&num_samples), sizeof(size_t));
        data_out.write(reinterpret_cast<const char*>(&input_size), sizeof(size_t));

        for (const auto& s : statesV) {
            data_out.write(reinterpret_cast<const char*>(s.data()), s.size() * sizeof(double));
        }
        data_out.write(reinterpret_cast<const char*>(labelsV.data()), labelsV.size() * sizeof(int));
        data_out.close();

        // --- Save metadata txt ---
        std::ofstream meta_out(dir_path + "/Data_Parameters.txt");
        meta_out << "stencil_size " << stencil_size << "\n";
        meta_out << "output_size " << output_size << "\n";
        meta_out << "no_batches " << no_batches << "\n";
        meta_out << "batch_size " << batch_size << "\n";
        meta_out << "include_Moments " << include_Moments << "\n";
        meta_out << "paraboloid_coeff_stddev " << paraboloid_coeff_stddev << "\n";
        meta_out << "sheet_coeff_stddev " << sheet_coeff_stddev << "\n";
        meta_out << "max_sheet_thickness " << max_sheet_thickness << "\n";
        meta_out << "sheet_thickness_stddev " << sheet_thickness_stddev << "\n";
        meta_out << "max_cylinder_radius " << max_cylinder_radius << "\n";
        meta_out << "cylinder_radius_stddev " << cylinder_radius_stddev << "\n";
        meta_out << "max_sphere_radius " << max_sphere_radius << "\n";
        meta_out << "sphere_radius_stddev " << sphere_radius_stddev << "\n";
        meta_out.close();

        std::cout << "💾 Dataset and parameters saved to " << dir_path << std::endl;
    }

    void loadDataset(const std::string& dir_path)
    {
        std::string data_file = dir_path /*+ "/data.bin"*/;
        std::ifstream data_in(data_file, std::ios::binary);
        if (!data_in) {
            throw std::runtime_error("Failed to open " + data_file);
        }

        size_t num_samples, input_size;
        data_in.read(reinterpret_cast<char*>(&num_samples), sizeof(size_t));
        data_in.read(reinterpret_cast<char*>(&input_size), sizeof(size_t));

        statesV.resize(num_samples, std::vector<double>(input_size));
        labelsV.resize(num_samples);

        for (size_t i = 0; i < num_samples; ++i) {
            data_in.read(reinterpret_cast<char*>(statesV[i].data()), input_size * sizeof(double));
        }
        data_in.read(reinterpret_cast<char*>(labelsV.data()), num_samples * sizeof(int));
        data_in.close();

        dataset_path = dir_path;

        std::cout << "📂 Loaded " << num_samples << " samples from " << dir_path << std::endl;
    }

    void trainModel(int epochs = 20, double learning_rate = 0.001, int batch_size = 64)
    {
        optimizer = std::make_unique<torch::optim::Adam>(
            net.parameters(), torch::optim::AdamOptions(learning_rate));

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
        IRL::Trainer<IRL::Net, Loader> trainer(
            net, *optimizer, train_loader.get(), test_loader.get(), val_loader.get(), epochs);
        trainer.train();
    }

    void saveModel(const std::string& path, bool save_optimizer = false) {
        namespace fs = std::filesystem;
        fs::create_directories(fs::path(path).parent_path());

        // Use OutputArchive to save custom nn::Module
        torch::serialize::OutputArchive archive;
        net.save(archive);
        archive.save_to(path);
        std::cout << "💾 Saved trained model to " << path << std::endl;

        // Optionally save optimizer
        if (save_optimizer && optimizer) {
            std::string opt_path = path + ".opt";
            torch::save(*optimizer, opt_path);
            std::cout << "💾 Saved optimizer state to " << opt_path << std::endl;
        }

        // Save dataset path info if available
        if (!dataset_path.empty()) {
            std::string ds_info_file = path + "_dataset_path.txt";
            std::ofstream ds_out(ds_info_file);
            ds_out << dataset_path << std::endl;
            ds_out.close();
            std::cout << "📝 Saved dataset path info to " << ds_info_file << std::endl;
        }
    }

    void loadModel(const std::string& path, bool load_optimizer = false) {
        torch::serialize::InputArchive archive;
        archive.load_from(path);
        net.load(archive);
        std::cout << "📂 Loaded trained model from " << path << std::endl;

        if (load_optimizer) {
            if (!optimizer) {
                optimizer = std::make_unique<torch::optim::Adam>(
                    net.parameters(), torch::optim::AdamOptions(learning_rate));
            }
            std::string opt_path = path + ".opt";
            torch::load(*optimizer, opt_path);
            std::cout << "📂 Loaded optimizer state from " << opt_path << std::endl;
        }
    }

    int classify(const std::vector<double>& flattened_state,
                 std::vector<float>* out_probs = nullptr) override
    {
        torch::NoGradGuard no_grad;
        auto input = torch::tensor(flattened_state,
                                   torch::TensorOptions().dtype(torch::kFloat32))
                         .unsqueeze(0);
        auto device = net.fc1->weight.device();
        input = input.to(device);

        torch::Tensor logits = net.forward(input);
        torch::Tensor probs =
            torch::softmax(logits, 1).to(torch::kCPU).squeeze(0);

        if (out_probs) {
            out_probs->resize(probs.size(0));
            for (int i = 0; i < probs.size(0); ++i)
                (*out_probs)[i] = probs[i].item<float>();
        }

        return probs.argmax().item<int>();
    }
};

} // namespace IRL
