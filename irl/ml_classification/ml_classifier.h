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
#include <cmath>

namespace IRL {

class MLClassifier : public Classifier {
protected:
    // Helper for determining network input size from the selected data layout
    static int calculateInputSize(
        int stencil,
        int include_moments,
        bool include_surface_area,
        bool include_eigenvalues)
    {
        // Per-cell layout:
        // include_moments = 0: [vfrac]
        // include_moments >= 1: [vfrac, mx, my, mz]
        // include_surface_area: + [area]
        const int per_cell_stride =
            (include_moments >= 1 ? 4 : 1)
            + (include_surface_area ? 1 : 0);

        int size =
            stencil * stencil * stencil * per_cell_stride;

        // Stencil-wide second moment tensor
        if (include_moments >= 2) {
            size += 6;
        }

        // Stencil-wide eigenvalues
        if (include_eigenvalues) {
            size += 3;
        }

        return size;
    }

    // Network Parameters
    int input_size = calculateInputSize(stencil_size, include_Moments, include_Surface_Area, include_Eigenvalues);
    int hidden_size1 = 256;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 6;


    // Training Parameters
    double learning_rate = 0.001;
    int batch_size = 64;
    int epochs = 50;
    int reduce_lr_patience = 4;
    int early_stop_patience = 8;

    // Data Parameters
    int no_batches = 4096 * 4;
    int no_datapoints = no_batches * batch_size;

    double paraboloid_coeff_stddev = 0.1;
    double hyperbolic_cylinder_opening_angle_stddev = 20.0;
    double sheet_coeff_stddev = 0.1;
    double sheet_thickness_stddev = 0.0;
    double cylinder_radius_stddev = 0.0;
    double radius_circle_min = 2.5;
    double radius_circle_max = 10.0;
    double sphere_radius_stddev = 0.0;
    double ellipsoid_subgrid_stddev = 0.7;
    double min_long_ellipsoid_axis = 3.0;
    double max_long_ellipsoid_axis = 5.0;

    bool exact_2nd_moment = false;
    bool visualize = false;

    double machineZero = 1e-12;
    double lower_limit_subgrid = machineZero;
    double upper_limit_subgrid = std::sqrt(3.0);
    double class0_max_characteristic = 2.5;

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
    MLClassifier(
        int stencil = 5,
        int include_moments = 1,
        bool include_surface_area = false,
        bool include_eigenvalues = false,
        int h1 = 256,
        int h2 = 64,
        int h3 = 32,
        int out = 6)
        : Classifier(stencil),
        input_size(calculateInputSize(
            stencil,
            include_moments,
            include_surface_area,
            include_eigenvalues)),
        hidden_size1(h1),
        hidden_size2(h2),
        hidden_size3(h3),
        output_size(out),
        net(input_size, h1, h2, h3, out)
    {
        include_Moments = include_moments;
        include_Surface_Area = include_surface_area;
        include_Eigenvalues = include_eigenvalues;
    }

    void updateTrainingParameters(
        double lr = 0.001,
        int bs = 64,
        int ep = 50,
        int reduce_lr_pat = 4,
        int early_stop_pat = 8)
    {
        learning_rate = lr;
        batch_size = bs;
        epochs = ep;
        reduce_lr_patience = reduce_lr_pat;
        early_stop_patience = early_stop_pat;

        // Number of generated datapoints depends on batch size
        no_datapoints = no_batches * batch_size;
    }

    void updateDataParameters(
        int no_batches_in,
        double paraboloid_coeff_stddev_in,
        double hyperbolic_cylinder_opening_angle_stddev_in,
        double sheet_coeff_stddev_in,
        double sheet_thickness_stddev_in,
        double cylinder_radius_stddev_in,
        double radius_circle_min_in,
        double radius_circle_max_in,
        double sphere_radius_stddev_in,
        double ellipsoid_subgrid_stddev_in,
        double min_long_ellipsoid_axis_in,
        double max_long_ellipsoid_axis_in,
        bool exact_2nd_mom = false,
        bool visualize_in = false,
        double machineZero_in = 1e-12,
        double lower_limit_subgrid_in = 1e-12,
        double upper_limit_subgrid_in = std::sqrt(3.0),
        double class0_max_characteristic_in = 2.5)
    {
        no_batches = no_batches_in;
        no_datapoints = no_batches * batch_size;

        paraboloid_coeff_stddev = paraboloid_coeff_stddev_in;
        hyperbolic_cylinder_opening_angle_stddev =
            hyperbolic_cylinder_opening_angle_stddev_in;
        sheet_coeff_stddev = sheet_coeff_stddev_in;
        sheet_thickness_stddev = sheet_thickness_stddev_in;
        cylinder_radius_stddev = cylinder_radius_stddev_in;
        radius_circle_min = radius_circle_min_in;
        radius_circle_max = radius_circle_max_in;
        sphere_radius_stddev = sphere_radius_stddev_in;
        ellipsoid_subgrid_stddev = ellipsoid_subgrid_stddev_in;
        min_long_ellipsoid_axis = min_long_ellipsoid_axis_in;
        max_long_ellipsoid_axis = max_long_ellipsoid_axis_in;

        exact_2nd_moment = exact_2nd_mom;
        visualize = visualize_in;
        machineZero = machineZero_in;
        lower_limit_subgrid = lower_limit_subgrid_in;
        upper_limit_subgrid = upper_limit_subgrid_in;
        class0_max_characteristic = class0_max_characteristic_in;
    }

    void generateDataset() {
        IRL::Data_gen data_gen;
        data_gen.set_stencil_size(stencil_size);
        data_gen.updateDataParameters(
            no_datapoints,
            include_Moments,
            include_Surface_Area,
            include_Eigenvalues,
            paraboloid_coeff_stddev,
            hyperbolic_cylinder_opening_angle_stddev,
            sheet_coeff_stddev,
            sheet_thickness_stddev,
            cylinder_radius_stddev,
            radius_circle_min,
            radius_circle_max,
            sphere_radius_stddev,
            ellipsoid_subgrid_stddev,
            min_long_ellipsoid_axis,
            max_long_ellipsoid_axis,
            exact_2nd_moment,
            visualize,
            machineZero,
            lower_limit_subgrid,
            upper_limit_subgrid,
            class0_max_characteristic
        );
        using namespace std::chrono;
        // Record start time
        auto start_time = high_resolution_clock::now();
        
        data_gen.generateData(&statesV, &labelsV);
        
        // Record end time
        auto end_time = high_resolution_clock::now();
        generation_time = duration_cast<seconds>(end_time - start_time).count();
    }

    void appendDataset(const std::string& existing_path, bool save_combined = false) {
        loadDataset(existing_path);

        std::vector<std::vector<float>> new_states;
        std::vector<int> new_labels;

        IRL::Data_gen data_gen;
        data_gen.updateDataParameters(
            no_datapoints,
            include_Moments,
            include_Surface_Area,
            include_Eigenvalues,
            paraboloid_coeff_stddev,
            hyperbolic_cylinder_opening_angle_stddev,
            sheet_coeff_stddev,
            sheet_thickness_stddev,
            cylinder_radius_stddev,
            radius_circle_min,
            radius_circle_max,
            sphere_radius_stddev,
            ellipsoid_subgrid_stddev,
            min_long_ellipsoid_axis,
            max_long_ellipsoid_axis,
            exact_2nd_moment,
            visualize,
            machineZero,
            lower_limit_subgrid,
            upper_limit_subgrid,
            class0_max_characteristic
        );

        using namespace std::chrono;
        auto start_time = high_resolution_clock::now();

        data_gen.generateData(&new_states, &new_labels);

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
        meta << "no_samples generated " << no_datapoints << "\n";
        meta << "include_Moments " << include_Moments << "\n";
        meta << "include_Surface_Area " << include_Surface_Area << "\n";
        meta << "include_Eigenvalues " << include_Eigenvalues << "\n";
        meta << "paraboloid_coeff_stddev " << paraboloid_coeff_stddev << "\n";
        meta << "hyperbolic_cylinder_opening_angle_stddev " << hyperbolic_cylinder_opening_angle_stddev << "\n";
        meta << "sheet_coeff_stddev " << sheet_coeff_stddev << "\n";
        meta << "sheet_thickness_stddev " << sheet_thickness_stddev << "\n";
        meta << "cylinder_radius_stddev " << cylinder_radius_stddev << "\n";
        meta << "radius_circle_min " << radius_circle_min << "\n";
        meta << "radius_circle_max " << radius_circle_max << "\n";
        meta << "sphere_radius_stddev " << sphere_radius_stddev << "\n";
        meta << "ellipsoid_subgrid_stddev " << ellipsoid_subgrid_stddev << "\n";
        meta << "min_long_ellipsoid_axis " << min_long_ellipsoid_axis << "\n";
        meta << "max_long_ellipsoid_axis " << max_long_ellipsoid_axis << "\n";
        meta << "exact_2nd_moment " << exact_2nd_moment << "\n";
        meta << "visualize " << visualize << "\n";
        meta << "machineZero " << machineZero << "\n";
        meta << "lower_limit_subgrid " << lower_limit_subgrid << "\n";
        meta << "upper_limit_subgrid " << upper_limit_subgrid << "\n";
        meta << "class0_max_characteristic " << class0_max_characteristic << "\n";
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

    void checkStatesForNaNOrInf(size_t max_print = 20, bool remove_bad_samples = false) {
        size_t total_bad_cases = 0;
        size_t printed_cases = 0;

        if (statesV.empty()) {
            std::cout << "No states data available to check." << std::endl;
            return;
        }

        if (statesV.size() != labelsV.size()) {
            std::cerr << "Warning: statesV size (" << statesV.size()
                    << ") does not match labelsV size (" << labelsV.size()
                    << "). Labels may be unavailable for some samples."
                    << std::endl;
        }

        std::vector<std::vector<float>> cleaned_states;
        std::vector<int> cleaned_labels;

        if (remove_bad_samples) {
            cleaned_states.reserve(statesV.size());
            cleaned_labels.reserve(labelsV.size());
        }

        for (size_t sample = 0; sample < statesV.size(); ++sample) {
            const auto& state = statesV[sample];

            bool bad_sample = false;
            size_t bad_feature = 0;
            float bad_value = 0.0f;

            for (size_t feature = 0; feature < state.size(); ++feature) {
                if (!std::isfinite(state[feature])) {
                    bad_sample = true;
                    bad_feature = feature;
                    bad_value = state[feature];
                    break;
                }
            }

            if (bad_sample) {
                ++total_bad_cases;

                if (printed_cases < max_print) {
                    std::cout << "\n=== Bad stencil/state found ===" << std::endl;
                    std::cout << "sample_index=" << sample << std::endl;

                    if (sample < labelsV.size()) {
                        std::cout << "label=" << labelsV[sample] << std::endl;
                    } else {
                        std::cout << "label=N/A" << std::endl;
                    }

                    std::cout << "bad_feature_index=" << bad_feature << std::endl;
                    std::cout << "bad_value=" << bad_value << std::endl;
                    std::cout << "state_size=" << state.size() << std::endl;

                    std::cout << "full_state=[" << std::endl;
                    for (size_t i = 0; i < state.size(); ++i) {
                        std::cout << "  [" << i << "] " << state[i];

                        if (!std::isfinite(state[i])) {
                            std::cout << "  <-- NaN/Inf";
                        }

                        std::cout << std::endl;
                    }
                    std::cout << "]" << std::endl;

                    ++printed_cases;
                }

                if (remove_bad_samples) {
                    continue;
                }
            }

            if (remove_bad_samples) {
                cleaned_states.push_back(state);

                if (sample < labelsV.size()) {
                    cleaned_labels.push_back(labelsV[sample]);
                }
            }
        }

        if (remove_bad_samples) {
            statesV = std::move(cleaned_states);
            labelsV = std::move(cleaned_labels);
            no_samples = statesV.size();
        }

        std::cout << "\nTotal samples containing NaN or Inf: "
                << total_bad_cases << std::endl;

        if (remove_bad_samples) {
            std::cout << "Removed samples: "
                    << total_bad_cases << std::endl;

            std::cout << "Remaining samples: "
                    << statesV.size() << std::endl;
        }
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

            preprocess_stencil(flat, stencil_size, no_symmetries, include_Moments, include_Surface_Area, include_Eigenvalues);
        }
        std::cout << "Length of flattened state: " << statesV[0].size() << std::endl;

        std::cout << "✅ Canonicalization complete!" << std::endl;
    }

    void preprocess_data(int no_canonical_symmetries, float noise_stddev = 0.0f) {
        if (statesV.empty()) {
            std::cerr << "⚠ No data loaded. Cannot preprocess." << std::endl;
            return;
        }

        std::cout << "🔧 Preprocessing " << statesV.size() << " samples..." << std::endl;
        std::cout << "Length of flattened state before: " << statesV[0].size() << std::endl;
        std::cout << "including Eigenvalues: " << include_Eigenvalues << ", including Moments: " << include_Moments << std::endl;
        std::cout << "including Surface Area: " << include_Surface_Area << std::endl;
        for (size_t sample = 0; sample < statesV.size(); sample++) {
            auto& flat = statesV[sample];

            preprocess_stencil(flat, stencil_size, no_canonical_symmetries, include_Moments, include_Surface_Area, include_Eigenvalues, noise_stddev);
        }
        std::cout << "Length of flattened state: " << statesV[0].size() << std::endl;

        std::cout << "✅ Preprocessing complete!" << std::endl;
    }


    double trainModel() {
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
        //return final_test_accuracy;
        return *std::min_element(val_loss.begin(), val_loss.end());
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
            param_out << "stencil_size " << stencil_size << "\n";
            param_out << "output_size " << output_size << "\n";
            param_out << "no_batches generated " << no_batches << "\n";
            param_out << "batch_size " << batch_size << "\n";
            param_out << "no_samples generated " << no_datapoints << "\n";
            param_out << "include_Moments " << include_Moments << "\n";
            param_out << "include_Surface_Area " << include_Surface_Area << "\n";
            param_out << "include_Eigenvalues " << include_Eigenvalues << "\n";
            param_out << "paraboloid_coeff_stddev " << paraboloid_coeff_stddev << "\n";
            param_out << "hyperbolic_cylinder_opening_angle_stddev " << hyperbolic_cylinder_opening_angle_stddev << "\n";
            param_out << "sheet_coeff_stddev " << sheet_coeff_stddev << "\n";
            param_out << "sheet_thickness_stddev " << sheet_thickness_stddev << "\n";
            param_out << "cylinder_radius_stddev " << cylinder_radius_stddev << "\n";
            param_out << "radius_circle_min " << radius_circle_min << "\n";
            param_out << "radius_circle_max " << radius_circle_max << "\n";
            param_out << "sphere_radius_stddev " << sphere_radius_stddev << "\n";
            param_out << "ellipsoid_subgrid_stddev " << ellipsoid_subgrid_stddev << "\n";
            param_out << "min_long_ellipsoid_axis " << min_long_ellipsoid_axis << "\n";
            param_out << "max_long_ellipsoid_axis " << max_long_ellipsoid_axis << "\n";
            param_out << "exact_2nd_moment " << exact_2nd_moment << "\n";
            param_out << "visualize " << visualize << "\n";
            param_out << "machineZero " << machineZero << "\n";
            param_out << "lower_limit_subgrid " << lower_limit_subgrid << "\n";
            param_out << "upper_limit_subgrid " << upper_limit_subgrid << "\n";
            param_out << "class0_max_characteristic " << class0_max_characteristic << "\n";
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

        if (flattened_state.size() != static_cast<size_t>(input_size)) {
            throw std::runtime_error(
                "MLClassifier::classify: expected input size "
                + std::to_string(input_size)
                + ", but received "
                + std::to_string(flattened_state.size())
            );
        }

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

    void exportRuntimeWeightsAndBiasesHeader(const std::string& filename = "ml_classifier_weights_and_biases.h") const {
        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Failed to open " + filename);
        }

        auto tensor_to_cpu_double = [](const torch::Tensor& tensor) {
            return tensor.detach()
                        .to(torch::kCPU)
                        .to(torch::kFloat64)
                        .contiguous();
        };

        auto write_vector = [&out, &tensor_to_cpu_double](
                                const std::string& name,
                                const torch::Tensor& tensor,
                                int size) {
            torch::Tensor cpu = tensor_to_cpu_double(tensor);

            if (cpu.dim() != 1 || cpu.size(0) != size) {
                throw std::runtime_error("Unexpected shape for " + name);
            }

            const double* data = cpu.data_ptr<double>();

            out << "static constexpr std::array<double, "
                << size << "> " << name << " {{\n";

            out << std::setprecision(std::numeric_limits<double>::max_digits10)
                << std::scientific;

            for (int i = 0; i < size; ++i) {
                if (!std::isfinite(data[i])) {
                    throw std::runtime_error("Non-finite value found in " + name);
                }

                out << "    " << data[i];

                if (i + 1 < size) {
                    out << ",";
                }

                out << "\n";
            }

            out << "}};\n\n";
        };

        auto write_matrix = [&out, &tensor_to_cpu_double](
                                const std::string& name,
                                const torch::Tensor& tensor,
                                int rows,
                                int cols) {
            torch::Tensor cpu = tensor_to_cpu_double(tensor);

            if (cpu.dim() != 2 || cpu.size(0) != rows || cpu.size(1) != cols) {
                throw std::runtime_error("Unexpected shape for " + name);
            }

            const double* data = cpu.data_ptr<double>();

            out << "static constexpr std::array<std::array<double, "
                << cols << ">, " << rows << "> " << name << " {{\n";

            out << std::setprecision(std::numeric_limits<double>::max_digits10)
                << std::scientific;

            for (int r = 0; r < rows; ++r) {
                out << "    {{";

                for (int c = 0; c < cols; ++c) {
                    const double value = data[r * cols + c];

                    if (!std::isfinite(value)) {
                        throw std::runtime_error("Non-finite value found in " + name);
                    }

                    out << value;

                    if (c + 1 < cols) {
                        out << ", ";
                    }
                }

                out << "}}";

                if (r + 1 < rows) {
                    out << ",";
                }

                out << "\n";
            }

            out << "}};\n\n";
        };

        out << "// This file was generated from a trained LibTorch model.\n";
        out << "// Do not edit by hand unless you know what you are doing.\n\n";

        out << "#ifndef IRL_ML_CLASSIFIER_WEIGHTS_AND_BIASES_H_\n";
        out << "#define IRL_ML_CLASSIFIER_WEIGHTS_AND_BIASES_H_\n\n";

        out << "#include <array>\n\n";

        out << "namespace IRL {\n";
        out << "namespace mlclassifier {\n\n";

        // Input/data layout
        out << "static constexpr int stencil_size = " << stencil_size << ";\n";
        out << "static constexpr int include_Moments = " << include_Moments << ";\n";
        out << "static constexpr bool include_Surface_Area = " << (include_Surface_Area ? "true" : "false") << ";\n";
        out << "static constexpr bool include_Eigenvalues = " << (include_Eigenvalues ? "true" : "false") << ";\n\n";

        // Network architecture
        out << "static constexpr int input_size = " << input_size << ";\n";
        out << "static constexpr int hidden_size1 = " << hidden_size1 << ";\n";
        out << "static constexpr int hidden_size2 = " << hidden_size2 << ";\n";
        out << "static constexpr int hidden_size3 = " << hidden_size3 << ";\n";
        out << "static constexpr int output_size = " << output_size << ";\n\n";

        write_matrix("fc1_weight", net.fc1->weight, hidden_size1, input_size);
        write_vector("fc1_bias",   net.fc1->bias,   hidden_size1);

        write_matrix("fc2_weight", net.fc2->weight, hidden_size2, hidden_size1);
        write_vector("fc2_bias",   net.fc2->bias,   hidden_size2);

        write_matrix("fc3_weight", net.fc3->weight, hidden_size3, hidden_size2);
        write_vector("fc3_bias",   net.fc3->bias,   hidden_size3);

        write_matrix("fc4_weight", net.fc4->weight, output_size, hidden_size3);
        write_vector("fc4_bias",   net.fc4->bias,   output_size);

        out << "}  // namespace mlclassifier\n";
        out << "}  // namespace IRL\n\n";

        out << "#endif  // IRL_ML_CLASSIFIER_WEIGHTS_AND_BIASES_H_\n";

        std::cout << "Exported compile-time runtime weights to "
                << filename << std::endl;
    }

    void retain_only_0th_moments() {
        if (statesV.empty()) {
            std::cerr
                << "⚠ No data loaded. Cannot retain only 0th moments."
                << std::endl;
            return;
        }

        std::cout
            << "🔧 Reducing " << statesV.size()
            << " states to volume fractions only..."
            << std::endl;

        std::cout
            << "Length of flattened state before: "
            << statesV.front().size()
            << std::endl;

        // These values describe the layout before the states are reduced.
        const int old_include_moments = include_Moments;
        const bool old_include_surface_area = include_Surface_Area;
        const bool old_include_eigenvalues = include_Eigenvalues;

        for (auto& flattened_state : statesV) {
            IRL::retain_only_0th_moments(
                flattened_state,
                stencil_size,
                old_include_moments,
                old_include_surface_area,
                old_include_eigenvalues
            );
        }

        // The states now contain only one volume fraction per stencil cell.
        include_Moments = 0;
        include_Surface_Area = false;
        include_Eigenvalues = false;

        input_size =
            stencil_size * stencil_size * stencil_size;

        std::cout
            << "Length of flattened state after: "
            << statesV.front().size()
            << std::endl;

        std::cout
            << "include_Moments: " << include_Moments
            << ", include_Surface_Area: " << include_Surface_Area
            << ", include_Eigenvalues: " << include_Eigenvalues
            << std::endl;

        std::cout
            << "✅ States now contain only 0th moments."
            << std::endl;
    }
};

} // namespace IRL
