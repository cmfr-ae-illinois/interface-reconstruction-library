#pragma once
#include "classifier.h"
#include "data_gen.h"
#include "data_set.h"
#include "trainer.h"
#include "net.h"
#include "stencil_rotator.h"

#include <torch/torch.h>
#include <iostream>
#include <fstream>
#include <filesystem>

namespace IRL {

class MLClassifier : public Classifier {
protected:
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
    int reduce_lr_patience = 4;
    int early_stop_patience = 8;

    //Data Parameters
    int no_batches = 256;
    int include_Moments = 0;
    double paraboloid_coeff_stddev = 0.1;
    double sheet_coeff_stddev = 0.1;
    double max_sheet_thickness = 0.5;
    double sheet_thickness_stddev = 0.0;
    double max_cylinder_radius = 0.5;
    double cylinder_radius_stddev = 0.0;
    bool include_truncated_cylinder = false;
    double max_sphere_radius = 0.5;
    double sphere_radius_stddev = 0.0;
    bool exact_2nd_moment = false;

    //Internals
    Net net;
    std::unique_ptr<torch::optim::Optimizer> optimizer;
    std::string dataset_path = "";
    std::vector<std::vector<float>> statesV;
    std::vector<int> labelsV;
    size_t no_samples = batch_size * no_batches;
    double generation_time = 0.0;

    //Training logs
    std::vector<double> train_loss;
    std::vector<double> val_loss;
    std::vector<double> test_accuracy_vec;
    std::vector<double> val_accuracy_vec;

    double final_test_loss = 0.0;
    double final_test_accuracy = 0.0;

    std::vector<std::vector<int64_t>> confusion_matrix;

public:
    MLClassifier(int stencil = 3, int input = 27,
                 int h1 = 128, int h2 = 64, int h3 = 32, int out = 4)
        : Classifier(stencil),
          stencil_size(stencil),
          input_size(input),
          hidden_size1(h1), hidden_size2(h2), hidden_size3(h3),
          output_size(out),
          net(input, h1, h2, h3, out) {}

    void updateTrainingParameters(double lr, int bs, int ep, int reduce_lr_pat = 4, int early_stop_pat = 8) {
        learning_rate = lr;
        batch_size = bs;
        epochs = ep;
        reduce_lr_patience = reduce_lr_pat;
        early_stop_patience = early_stop_pat;
    }

    void updateDataParameters(int nb, int incMoments,
                              double parab_std, double sheet_std,
                              double max_sheet_th, double sheet_th_std,
                              double max_cyl_r, double cyl_r_std, bool incl_trunc_cyl,
                              double max_sph_r, double sph_r_std,
                              bool exact_2nd_mom) {
        no_batches = nb;
        include_Moments = incMoments;
        paraboloid_coeff_stddev = parab_std;
        sheet_coeff_stddev = sheet_std;
        max_sheet_thickness = max_sheet_th;
        sheet_thickness_stddev = sheet_th_std;
        max_cylinder_radius = max_cyl_r;
        cylinder_radius_stddev = cyl_r_std;
        include_truncated_cylinder = incl_trunc_cyl;
        max_sphere_radius = max_sph_r;
        sphere_radius_stddev = sph_r_std;
        exact_2nd_moment = exact_2nd_mom;
    }

    void generateDataset() {
        IRL::Data_gen data_gen;
        using namespace std::chrono;
        // Record start time
        auto start_time = high_resolution_clock::now();
        
        data_gen.generateData(&statesV, &labelsV,
                               no_batches * batch_size,
                               stencil_size, output_size,
                               include_Moments,
                               paraboloid_coeff_stddev,
                               sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev,
                               max_cylinder_radius, cylinder_radius_stddev, include_truncated_cylinder,
                               max_sphere_radius, sphere_radius_stddev);
        /*
        void generateData (std::vector<std::vector<double>>* statesV, std::vector<int>* labelsV, int no_datapoints, int stencil_size = 3, int no_datapoint_types_in = 4, int include_Moments = 0,
                                        double paraboloid_coeff_stddev = 0.1,
                                        double sheet_coeff_stddev = 0.1, double max_sheet_thickness = 0.5, double sheet_thickness_stddev = 0.0,
                                        double max_cylinder_radius = 0.5, double cylinder_radius_stddev = 0.0, bool include_truncated_cylinder = false,
                                        double max_sphere_radius = 0.5, double sphere_radius_stddev = 0.0)

        */
        /*
        data_gen.generateData(&statesV, &labelsV,
                               no_batches * batch_size,
                               stencil_size);
        */
        // Record end time
        auto end_time = high_resolution_clock::now();
        generation_time = duration_cast<seconds>(end_time - start_time).count();

        saveDataset("data");
    }

    void appendDataset(const std::string& existing_path, bool save_combined = false) {
        loadDataset(existing_path);

        std::vector<std::vector<float>> new_states;
        std::vector<int> new_labels;

        IRL::Data_gen data_gen;

        using namespace std::chrono;
        // Record start time
        auto start_time = high_resolution_clock::now();

        data_gen.generateData(&statesV, &labelsV,
                               no_batches * batch_size,
                               stencil_size, output_size,
                               include_Moments,
                               paraboloid_coeff_stddev,
                               sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev,
                               max_cylinder_radius, cylinder_radius_stddev, include_truncated_cylinder,
                               max_sphere_radius, sphere_radius_stddev);
        // Record end time
        auto end_time = high_resolution_clock::now();
        generation_time = duration_cast<seconds>(end_time - start_time).count();

        statesV.insert(statesV.end(), new_states.begin(), new_states.end());
        labelsV.insert(labelsV.end(), new_labels.begin(), new_labels.end());

        std::cout << "✅ Appended " << new_states.size()
                  << " new samples. Total: " << statesV.size() << std::endl;

        if (save_combined)
            saveDataset("data");
    }

    void saveDataset(const std::string& dir_path, double time=0.0) const {
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
            data_out.write(reinterpret_cast<const char*>(s.data()), s.size() * sizeof(float));
        data_out.write(reinterpret_cast<const char*>(labelsV.data()), labelsV.size() * sizeof(int));
        data_out.close();

        std::ofstream meta(dir_path + "/Data_Parameters.txt");
        meta << "stencil_size " << stencil_size << "\n";
        meta << "output_size " << output_size << "\n";
        meta << "no_batches generated " << no_batches << "\n";
        meta << "batch_size " << batch_size << "\n";
        meta << "no_samples generated " << no_batches*batch_size << "\n";
        meta << "generation_time " << generation_time << "\n";
        meta << "no_batches loaded " << (no_batches - num_samples/batch_size) << "\n";
        meta << "total samples " << num_samples << "\n";
        meta << "include_Moments " << include_Moments << "\n";
        meta << "paraboloid_coeff_stddev " << paraboloid_coeff_stddev << "\n";
        meta << "sheet_coeff_stddev " << sheet_coeff_stddev << "\n";
        meta << "max_sheet_thickness " << max_sheet_thickness << "\n";
        meta << "sheet_thickness_stddev " << sheet_thickness_stddev << "\n";
        meta << "max_cylinder_radius " << max_cylinder_radius << "\n";
        meta << "cylinder_radius_stddev " << cylinder_radius_stddev << "\n";
        meta << "include_truncated_cylinder " << include_truncated_cylinder << "\n";
        meta << "max_sphere_radius " << max_sphere_radius << "\n";
        meta << "sphere_radius_stddev " << sphere_radius_stddev << "\n";
        meta << "exact_2nd_moment " << exact_2nd_moment << "\n";
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
        statesV.resize(num_samples, std::vector<float>(in_size));
        labelsV.resize(num_samples);
        for (size_t i = 0; i < num_samples; ++i)
            data_in.read(reinterpret_cast<char*>(statesV[i].data()), in_size * sizeof(float));
        data_in.read(reinterpret_cast<char*>(labelsV.data()), num_samples * sizeof(int));
        data_in.close();
        dataset_path = dir_path;

        std::cout << "📂 Loaded " << num_samples << " samples from " << dir_path << std::endl;
    }

    void canonicalize_data(int no_symmetries) { //OLD kept for compatability
        if (statesV.empty()) {
            std::cerr << "⚠ No data loaded. Cannot canonicalize." << std::endl;
            return;
        }

        std::cout << "🔄 Canonicalizing " << statesV.size()
                << " samples using " << no_symmetries << " symmetries..." << std::endl;

        for (size_t sample = 0; sample < statesV.size(); sample++) {
            auto& flat = statesV[sample];

            IRL::preprocess_stencil(flat, stencil_size, no_symmetries, include_Moments, false);
        }
        std::cout << "Length of flattened state: " << statesV[0].size() << std::endl;

        std::cout << "✅ Canonicalization complete!" << std::endl;
    }

    void preprocess_data(int no_canonical_symmetries, bool add_noise = false) {
        if (statesV.empty()) {
            std::cerr << "⚠ No data loaded. Cannot preprocess." << std::endl;
            return;
        }

        std::cout << "🔧 Preprocessing " << statesV.size() << " samples..." << std::endl;

        for (size_t sample = 0; sample < statesV.size(); sample++) {
            auto& flat = statesV[sample];

            IRL::preprocess_stencil(flat, stencil_size, no_canonical_symmetries, include_Moments, add_noise);
        }
        std::cout << "Length of flattened state: " << statesV[0].size() << std::endl;

        std::cout << "✅ Preprocessing complete!" << std::endl;
    }


    void trainModel() {
        optimizer = std::make_unique<torch::optim::AdamW>(
            net.parameters(),
            torch::optim::AdamWOptions(learning_rate)
                .weight_decay(1e-4)     // start here
                .betas(std::make_tuple(0.9, 0.999))
                .eps(1e-8)
        );

        int total_samples = statesV.size();
        int train_end = total_samples * 0.7;
        int test_end = total_samples * 0.85;

        IRL::MyDataset train_dataset(&statesV, &labelsV, 0, train_end);
        IRL::MyDataset test_dataset(&statesV, &labelsV, train_end, test_end);
        IRL::MyDataset val_dataset(&statesV, &labelsV, test_end, total_samples);

        auto train_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
            std::move(train_dataset).map(torch::data::transforms::Stack<>()), batch_size);
        auto test_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
            std::move(test_dataset).map(torch::data::transforms::Stack<>()), batch_size);
        auto val_loader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(
            std::move(val_dataset).map(torch::data::transforms::Stack<>()), batch_size);

        
        using Loader = std::remove_reference_t<decltype(*train_loader)>;
        IRL::Trainer<IRL::Net, Loader> trainer(
            net, *optimizer,
            train_loader.get(), test_loader.get(), val_loader.get(),
            epochs,
            reduce_lr_patience,
            early_stop_patience,
            &train_loss,
            &val_loss,
            &test_accuracy_vec,
            &val_accuracy_vec,
            &final_test_loss,
            &final_test_accuracy,
            &confusion_matrix
        );

        trainer.train();
    }

    void outputTrainingResults() const {
        std::cout << "=== Training Results ===" << std::endl;

        std::cout << "train_loss:" << std::endl;
        for (double v : train_loss) std::cout << v << std::endl;

        std::cout << "val_loss:" << std::endl;
        for (double v : val_loss) std::cout << v << std::endl;

        std::cout << "test_accuracy_vec:" << std::endl;
        for (double v : test_accuracy_vec) std::cout << v << std::endl;

        std::cout << "val_accuracy_vec:" << std::endl;
        for (double v : val_accuracy_vec) std::cout << v << std::endl;

        std::cout << "final_test_loss:" << std::endl;
        std::cout << final_test_loss << std::endl;

        std::cout << "final_test_accuracy:" << std::endl;
        std::cout << final_test_accuracy << std::endl;

        std::cout << "confusion_matrix:" << std::endl;

        for (size_t i = 0; i < confusion_matrix.size(); ++i) {
            for (size_t j = 0; j < confusion_matrix[i].size(); ++j) {
                std::cout << confusion_matrix[i][j];
                if (j + 1 < confusion_matrix[i].size())
                    std::cout << " ";
            }
            std::cout << std::endl;
        }
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

        param_out << "[TrainingParameters]\n";
        param_out << "learning_rate " << learning_rate << "\n";
        param_out << "batch_size " << batch_size << "\n";
        param_out << "epochs " << epochs << "\n\n";

        param_out << "[TrainingResults]\n";
        param_out << "train_loss:\n";
        for (double v : train_loss)
            param_out << v << "\n";
        param_out << "val_loss:\n";
        for (double v : val_loss)
            param_out << v << "\n";
        param_out << "test_accuracy_vec:\n";
        for (double v : test_accuracy_vec)
            param_out << v << "\n";
        param_out << "val_accuracy_vec:\n";
        for (double v : val_accuracy_vec)
            param_out << v << "\n";
        param_out << "final_test_loss:\n";
        param_out << final_test_loss << "\n";
        param_out << "final_test_accuracy:\n";
        param_out << final_test_accuracy << "\n";
        param_out << "confusion_matrix:\n";
        for (size_t i = 0; i < confusion_matrix.size(); ++i) {
            for (size_t j = 0; j < confusion_matrix[i].size(); ++j) {
                param_out << confusion_matrix[i][j];
                if (j + 1 < confusion_matrix[i].size())
                    param_out << " ";
            }
            param_out << "\n";
        }

        param_out << "\n";

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

    int classify(const std::vector<float>& flattened_state,
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
