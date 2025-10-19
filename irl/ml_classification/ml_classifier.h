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
    //Network Parameters
    int stencil_size = 3;
    int input_size = stencil_size * stencil_size * stencil_size; // 27 if stencil_size=3 and only vof
    int hidden_size1 = 128;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 4;

    //Training Parameters
    double learning_rate = 0.001; // was 0.01 for SGD
    int batch_size = 64;
    int epochs = 20;

    //Data Parameters
    int no_batches = 256;
    int include_Moments = 0;
    double paraboloid_coeff_stddev = 0.1;
    double sheet_coeff_stddev = 0.1;
    double max_sheet_thickness = 0.5;
    double sheet_thickness_stddev = 0.0;
    double max_cylinder_radius = 0.5;
    double cylinder_radius_stddev = 0.0;
    double max_sphere_radius = 0.5;
    double sphere_radius_stddev = 0.0;

    //Internals
    Net net;
    std::unique_ptr<torch::optim::Optimizer> optimizer;
    std::string dataset_path = "";
    std::vector<std::vector<double>> statesV;
    std::vector<int> labelsV;

public:
    MLClassifier(int stencil = 3, int input = 27,
                 int h1 = 128, int h2 = 64, int h3 = 32, int out = 4)
        : Classifier(stencil),
          stencil_size(stencil),
          input_size(input),
          hidden_size1(h1), hidden_size2(h2), hidden_size3(h3),
          output_size(out),
          net(input, h1, h2, h3, out) {}

    void updateTrainingParameters(double lr, int bs, int ep) {
        learning_rate = lr;
        batch_size = bs;
        epochs = ep;
    }

    void updateDataParameters(int nb, int incMoments,
                              double parab_std, double sheet_std,
                              double max_sheet_th, double sheet_th_std,
                              double max_cyl_r, double cyl_r_std,
                              double max_sph_r, double sph_r_std) {
        no_batches = nb;
        include_Moments = incMoments;
        paraboloid_coeff_stddev = parab_std;
        sheet_coeff_stddev = sheet_std;
        max_sheet_thickness = max_sheet_th;
        sheet_thickness_stddev = sheet_th_std;
        max_cylinder_radius = max_cyl_r;
        cylinder_radius_stddev = cyl_r_std;
        max_sphere_radius = max_sph_r;
        sphere_radius_stddev = sph_r_std;
    }

    void generateDataset() {
        IRL::Data_gen data_gen;
        data_gen.generate_Data(&statesV, &labelsV,
                               no_batches * batch_size,
                               stencil_size, output_size,
                               include_Moments,
                               paraboloid_coeff_stddev,
                               sheet_coeff_stddev,
                               max_sheet_thickness, sheet_thickness_stddev,
                               max_cylinder_radius, cylinder_radius_stddev,
                               max_sphere_radius, sphere_radius_stddev);
        saveDataset("data");
    }

    void appendDataset(const std::string& existing_path, bool save_combined = false) {
        loadDataset(existing_path);
        //size_t old_count = statesV.size();

        std::vector<std::vector<double>> new_states;
        std::vector<int> new_labels;

        IRL::Data_gen data_gen;
        data_gen.generate_Data(&new_states, &new_labels,
                               no_batches * batch_size,
                               stencil_size, output_size,
                               include_Moments,
                               paraboloid_coeff_stddev,
                               sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev,
                               max_cylinder_radius, cylinder_radius_stddev,
                               max_sphere_radius, sphere_radius_stddev);

        statesV.insert(statesV.end(), new_states.begin(), new_states.end());
        labelsV.insert(labelsV.end(), new_labels.begin(), new_labels.end());

        std::cout << "✅ Appended " << new_states.size()
                  << " new samples. Total: " << statesV.size() << std::endl;

        if (save_combined)
            saveDataset("data");
    }

    void saveDataset(const std::string& dir_path) const {
        namespace fs = std::filesystem;
        fs::create_directories(dir_path);

        std::string data_file = dir_path + "/data.bin";
        std::ofstream data_out(data_file, std::ios::binary);
        if (!data_out)
            throw std::runtime_error("Failed to open " + data_file);

        size_t num_samples = statesV.size();
        size_t in_size = (num_samples > 0) ? statesV[0].size() : 0;
        data_out.write(reinterpret_cast<const char*>(&num_samples), sizeof(size_t));
        data_out.write(reinterpret_cast<const char*>(&in_size), sizeof(size_t));
        for (const auto& s : statesV)
            data_out.write(reinterpret_cast<const char*>(s.data()), s.size() * sizeof(double));
        data_out.write(reinterpret_cast<const char*>(labelsV.data()), labelsV.size() * sizeof(int));
        data_out.close();

        std::ofstream meta(dir_path + "/Data_Parameters.txt");
        meta << "stencil_size " << stencil_size << "\n";
        meta << "output_size " << output_size << "\n";
        meta << "no_batches " << no_batches << "\n";
        meta << "batch_size " << batch_size << "\n";
        meta << "include_Moments " << include_Moments << "\n";
        meta << "paraboloid_coeff_stddev " << paraboloid_coeff_stddev << "\n";
        meta << "sheet_coeff_stddev " << sheet_coeff_stddev << "\n";
        meta << "max_sheet_thickness " << max_sheet_thickness << "\n";
        meta << "sheet_thickness_stddev " << sheet_thickness_stddev << "\n";
        meta << "max_cylinder_radius " << max_cylinder_radius << "\n";
        meta << "cylinder_radius_stddev " << cylinder_radius_stddev << "\n";
        meta << "max_sphere_radius " << max_sphere_radius << "\n";
        meta << "sphere_radius_stddev " << sphere_radius_stddev << "\n";
        meta.close();

        std::cout << "💾 Dataset and parameters saved to " << dir_path << std::endl;
    }

    void loadDataset(const std::string& dir_path) {
        std::string data_file = dir_path;
        std::ifstream data_in(data_file, std::ios::binary);
        if (!data_in)
            throw std::runtime_error("Failed to open " + data_file);

        size_t num_samples, in_size;
        data_in.read(reinterpret_cast<char*>(&num_samples), sizeof(size_t));
        data_in.read(reinterpret_cast<char*>(&in_size), sizeof(size_t));
        statesV.resize(num_samples, std::vector<double>(in_size));
        labelsV.resize(num_samples);
        for (size_t i = 0; i < num_samples; ++i)
            data_in.read(reinterpret_cast<char*>(statesV[i].data()), in_size * sizeof(double));
        data_in.read(reinterpret_cast<char*>(labelsV.data()), num_samples * sizeof(int));
        data_in.close();
        dataset_path = dir_path;

        std::cout << "📂 Loaded " << num_samples << " samples from " << dir_path << std::endl;
    }

    void trainModel() {
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

    void saveModel(const std::string& path, bool save_optimizer_state = false) {
        namespace fs = std::filesystem;
        fs::create_directories(fs::path(path).parent_path());

        torch::serialize::OutputArchive archive;
        net.save(archive);
        archive.save_to(path + "ml_model.pt");
        std::cout << "💾 Saved trained model to " << path << std::endl;

        if (save_optimizer_state && optimizer) {
            std::string opt_path = path + ".opt";
            torch::save(*optimizer, opt_path);
            std::cout << "💾 Saved optimizer state to " << opt_path << std::endl;
        }

        std::string param_file = path + "_Parameters.txt";
        std::ofstream param_out(param_file);
        if (!param_out) {
            throw std::runtime_error("Failed to open " + param_file);
        }
        param_out << "[Dataset]\n";

        if (!dataset_path.empty()) {
            param_out << dataset_path << std::endl;
        } else {
            param_out << "no_batches " << no_batches << "\n";
            param_out << "include_Moments " << include_Moments << "\n";
            param_out << "paraboloid_coeff_stddev " << paraboloid_coeff_stddev << "\n";
            param_out << "sheet_coeff_stddev " << sheet_coeff_stddev << "\n";
            param_out << "max_sheet_thickness " << max_sheet_thickness << "\n";
            param_out << "sheet_thickness_stddev " << sheet_thickness_stddev << "\n";
            param_out << "max_cylinder_radius " << max_cylinder_radius << "\n";
            param_out << "cylinder_radius_stddev " << cylinder_radius_stddev << "\n";
            param_out << "max_sphere_radius " << max_sphere_radius << "\n";
            param_out << "sphere_radius_stddev " << sphere_radius_stddev << "\n";
        }

        param_out << "[Network]\n";
        param_out << "stencil_size " << stencil_size << "\n";
        param_out << "input_size " << input_size << "\n";
        param_out << "hidden_size1 " << hidden_size1 << "\n";
        param_out << "hidden_size2 " << hidden_size2 << "\n";
        param_out << "hidden_size3 " << hidden_size3 << "\n";
        param_out << "output_size " << output_size << "\n\n";

        param_out << "[Training]\n";
        param_out << "learning_rate " << learning_rate << "\n";
        param_out << "batch_size " << batch_size << "\n";
        param_out << "epochs " << epochs << "\n\n";

        param_out.close();
        std::cout << "📝 Saved model parameters to " << param_file << std::endl;
    }

    void loadModel(const std::string& path, bool load_optimizer_state = false) {
        torch::serialize::InputArchive archive;
        archive.load_from(path);
        net.load(archive);
        std::cout << "📂 Loaded trained model from " << path << std::endl;

        if (load_optimizer_state) {
            if (!optimizer)
                optimizer = std::make_unique<torch::optim::Adam>(
                    net.parameters(), torch::optim::AdamOptions(learning_rate));
            std::string opt_path = path + ".opt";
            torch::load(*optimizer, opt_path);
            std::cout << "📂 Loaded optimizer state from " << opt_path << std::endl;
        }
    }

    int classify(const std::vector<double>& flattened_state,
                 std::vector<float>* out_probs = nullptr) override {
        torch::NoGradGuard no_grad;
        auto input = torch::tensor(flattened_state,
                                   torch::TensorOptions().dtype(torch::kFloat32))
                         .unsqueeze(0);
        auto device = net.fc1->weight.device();
        input = input.to(device);

        torch::Tensor logits = net.forward(input);
        torch::Tensor probs = torch::softmax(logits, 1).to(torch::kCPU).squeeze(0);

        if (out_probs) {
            out_probs->resize(probs.size(0));
            for (int i = 0; i < probs.size(0); ++i)
                (*out_probs)[i] = probs[i].item<float>();
        }

        return probs.argmax().item<int>();
    }
};

} // namespace IRL
