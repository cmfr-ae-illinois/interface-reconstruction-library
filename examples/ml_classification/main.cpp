#include <vtkCellCenters.h>
#include "testcases.h"
#include "irl/ml_classification/ml_classifier.h"
#include "irl/ml_classification/vtk_in.h"
#include "irl/ml_classification/inertia_classifier.h"
#include "irl/ml_classification/hybid_classifier.h"
#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/ml_classifier_e3nn.h"
#include "irl/ml_classification/ml_classifier_notorch.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;


void find_dataset_size() {
    int stencil_size = 5;

    //Data parameters
    int no_batches;
    int include_Moments = 1;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = true;
    double paraboloid_coeff_stddev = 0.1;
    double hyperbolic_cylinder_opening_angle_stddev = 20; //degrees
    double sheet_coeff_stddev = 0.1;
    double sheet_thickness_stddev = 0.0;
    double cylinder_radius_stddev = 0.0;
    double radius_circle_min = 2.5;
    double radius_circle_max = 10.0;
    double sphere_radius_stddev = 0.0;
    double ellipsoid_subgrid_stddev = 0.7;
    double min_long_ellipsoid_axis = 3.0;
    double max_long_ellipsoid_axis = 5.0;
    bool exact_2nd_moment = false;  // enable calculation of exact 2nd moments for data generation
    bool visualize = false; // if true, print centroids and / or write surfaces
    double machineZero = 1e-12;
    double lower_limit_subgrid = machineZero;
    double upper_limit_subgrid = std::sqrt(3.0);
    double class0_max_characteristic = 2.5;
    float epsilon_connectivity = 1e-12f;

    // Net Parameters
    int input_size = stencil_size * stencil_size * stencil_size 
    * (include_Moments >= 1 ? (include_Surface_Area ? 5 : 4) : 1)  // 4 if include_Moments >= 1 because we have vfrac + (mx,my,mz) per cell, otherwise just vfrac
    + (include_Moments >= 2 ? 6 : 0)  // +6 if include_Moments >= 2 because we have (xx, yy, zz, xy, xz, yz) components of the 2nd moment tensor; otherwise none
    + (include_Eigenvalues ? 3 : 0); // +3 if include_Eigenvalues because we add the 3 eigenvalues of the inertia matrix; otherwise none
    int hidden_size1 = 256;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 6; //CHANGED 4 to 6

    // Training parameters
    double learning_rate = 0.001;
    int batch_size = 64;
    int max_epochs = 50;
    int reduce_lr_patience = 4;
    int early_stop_patience = 6;

    // Preprocessing
    int canonicalize_symmetries = 48;
    float noise_stddev = 0.0f;

    // Study parameters
    const int small_batch_increment = 16;       
    const int large_batch_increment = 128;
    const int switch_to_large_after_no_batches = 128;
    const int max_total_large_increments = 20; //= total 4*1024k *64 = 262k
                                        
    const int small_steps =
        switch_to_large_after_no_batches / small_batch_increment;

    const int max_steps =
        small_steps + max_total_large_increments -1;

    const double min_improvement = 0.005; // absolute accuracy improvement

    const std::string base_dir =
        "/home/quirin/mlcfd/Datasets/SixClasses/Thesis/DatasetSizeStudy/";
    fs::create_directories(base_dir);

    std::vector<double> accuracies;
    std::vector<int> batch_counts;
    std::vector<long long> sample_counts;

    std::string previous_dataset_file;

    for (int step = 1; step <= max_steps; ++step) {
        int total_batches;

        if (step <= small_steps) {
            total_batches = step * small_batch_increment;
            no_batches = small_batch_increment;
        } else {
            int large_step = step - small_steps;
            total_batches =
                switch_to_large_after_no_batches
                + large_step * large_batch_increment;
            no_batches = large_batch_increment;
        }

        long long total_samples =
            static_cast<long long>(total_batches) * batch_size;

        std::string dataset_name =
            "s5size" + std::to_string(total_batches) + "x" + std::to_string(batch_size);
        const std::string dataset_dir = base_dir + dataset_name;
        const std::string dataset_file = dataset_dir + "/data.bin";

        std::cout << "\n========================================\n";
        std::cout << "Step " << step << " / " << max_steps << "\n";
        std::cout << "Dataset       : " << dataset_name << "\n";
        std::cout << "Total batches : " << total_batches << "\n";
        std::cout << "Total samples : " << total_samples << "\n";
        std::cout << "========================================\n";

        double final_acc = 0.0;

        {
            IRL::MLClassifier ml(
                stencil_size, input_size,
                hidden_size1, hidden_size2, hidden_size3, output_size
            );

            // For generation/appending, no_batches is only the increment
            ml.updateDataParameters(
                no_batches,
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

            ml.updateTrainingParameters(
                learning_rate, batch_size, max_epochs,
                reduce_lr_patience, early_stop_patience
            );

            if (step == 1) {
                ml.generateDataset();        // generates 2048 * 64 samples
            } else {
                ml.appendDataset(previous_dataset_file, false); // load old + append 2048*64
            }

            ml.saveDataset(dataset_dir);     // save current full dataset under new name
            previous_dataset_file = dataset_file;

            ml.preprocess_data(canonicalize_symmetries, noise_stddev);
            final_acc = ml.trainModel();

            std::cout << "Final validation loss: " << final_acc << "\n";
        } // ml destroyed here, freeing dataset/model memory

        accuracies.push_back(final_acc);
        batch_counts.push_back(total_batches);
        sample_counts.push_back(total_samples);

        {
            std::ofstream out(base_dir + "dataset_size_results.txt");
            out << std::fixed << std::setprecision(6);
            out << "step total_batches total_samples final_test_accuracy\n";
            for (size_t i = 0; i < accuracies.size(); ++i) {
                out << (i + 1) << " "
                    << batch_counts[i] << " "
                    << sample_counts[i] << " "
                    << accuracies[i] << "\n";
            }
        }

        if (accuracies.size() >= 2) {
            const double improvement =
                accuracies.back() - accuracies[accuracies.size() - 2];

            std::cout << "Improvement over previous size: "
                      << improvement << "\n";
            /*
            if (improvement < min_improvement) {
                std::cout << "Stopping criterion reached.\n";
                break;
            }*/
        }
    }

    std::cout << "\n=== Dataset size study finished ===\n";
    for (size_t i = 0; i < accuracies.size(); ++i) {
        std::cout << "Batches: " << batch_counts[i]
                  << ", Samples: " << sample_counts[i]
                  << ", Final validation loss: " << accuracies[i] << "\n"; // changed from accuracy to loss in ml.train() return, variable names in this example are not adjusted
    }
}

// Fraction of equal entries between two prediction vectors
static double agreement_fraction(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("Prediction vectors have different sizes; align by cellId instead.");
    }
    if (a.empty()) return 0.0;

    size_t same = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        same += (a[i] == b[i]);
    }
    return static_cast<double>(same) / static_cast<double>(a.size());
}

// Pick the run whose predictions have the highest mean agreement with all other runs
static int pick_most_agreeing_model(const std::vector<std::vector<int>>& predictions, double* out_mean_agreement = nullptr) {
    const int R = static_cast<int>(predictions.size());
    if (R == 0) {
        if (out_mean_agreement) *out_mean_agreement = 0.0;
        return -1;
    }
    if (R == 1) {
        if (out_mean_agreement) *out_mean_agreement = 1.0;
        return 0;
    }

    const size_t N = predictions[0].size();
    for (int r = 1; r < R; ++r) {
        if (predictions[r].size() != N) {
            throw std::runtime_error("Prediction vectors differ in length!");
        }
    }

    int most_agreeing_run = 0;
    double best_score = -1.0;

    for (int r = 0; r < R; ++r) {
        double sum_agree = 0.0;
        for (int q = 0; q < R; ++q) {
            if (q == r) continue;
            sum_agree += agreement_fraction(predictions[r], predictions[q]);
        }
        const double mean_agree = sum_agree / static_cast<double>(R - 1);

        if (mean_agree > best_score) {
            best_score = mean_agree;
            most_agreeing_run = r;
        }
    }

    if (out_mean_agreement) *out_mean_agreement = best_score;
    return most_agreeing_run;
}

// Compute and print mean per-cell instability (can be removed if not needed)
static double compute_mean_instability(const std::vector<std::vector<int>>& predictions) {
    const int R = static_cast<int>(predictions.size());
    if (R == 0) return 0.0;

    const size_t N = predictions[0].size();
    for (int r = 1; r < R; ++r) {
        if (predictions[r].size() != N) {
            throw std::runtime_error("Prediction vectors differ in length!");
        }
    }

    double sum_u = 0.0;
    for (size_t t = 0; t < N; ++t) {
        std::array<int,4> counts{0,0,0,0};
        for (int r = 0; r < R; ++r) {
            int c = predictions[r][t];
            if (c < 0 || c >= 4) continue; // should not happen
            counts[c]++;
        }
        int max_count = std::max(std::max(counts[0], counts[1]), std::max(counts[2], counts[3]));
        double u = 1.0 - static_cast<double>(max_count) / static_cast<double>(R);
        sum_u += u;
    }

    return (N > 0) ? (sum_u / static_cast<double>(N)) : 0.0;
}
/*
void stable_classification() {
    int stencil_size = 5;

    //Data parameters
    int no_batches;
    int include_Moments = 1;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = false;
    double paraboloid_coeff_stddev = 0.1;
    double hyperbolic_cylinder_opening_angle_stddev = 20; //degrees
    double sheet_coeff_stddev = 0.1;
    double sheet_thickness_stddev = 0.0;
    double cylinder_radius_stddev = 0.0;
    double radius_circle_min = 2.5;
    double radius_circle_max = 10.0;
    double sphere_radius_stddev = 0.0;
    double ellipsoid_subgrid_stddev = 0.7;
    double min_long_ellipsoid_axis = 3.0;
    double max_long_ellipsoid_axis = 5.0;
    bool exact_2nd_moment = false;  // enable calculation of exact 2nd moments for data generation
    bool visualize = false; // if true, print centroids and / or write surfaces
    double machineZero = 1e-12;
    double lower_limit_subgrid = machineZero;
    double upper_limit_subgrid = std::sqrt(3.0);
    double class0_max_characteristic = 2.5;
    float epsilon_connectivity = 1e-12f;

    // Net Parameters
    int input_size = stencil_size * stencil_size * stencil_size 
    * (include_Moments >= 1 ? (include_Surface_Area ? 5 : 4) : 1)  // 4 if include_Moments >= 1 because we have vfrac + (mx,my,mz) per cell, otherwise just vfrac
    + (include_Moments >= 2 ? 6 : 0)  // +6 if include_Moments >= 2 because we have (xx, yy, zz, xy, xz, yz) components of the 2nd moment tensor; otherwise none
    + (include_Eigenvalues ? 3 : 0); // +3 if include_Eigenvalues because we add the 3 eigenvalues of the inertia matrix; otherwise none
    int hidden_size1 = 256;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 6; //CHANGED 4 to 6

    //Training parameters
    double learning_rate = 0.001; //was 0.01 for SGD optimizer
    int batch_size = 64;
    int max_epochs = 50;
    int reduce_lr_patience = 4;
    int early_stop_patience = 8;

    // Classification parameters
    int canonicalize_symmetries = 48;
    const int no_runs = 10;
    float noise_stddev = 0.0f;
    int downsample_factor = 2; // set to 2 to only classify every 2nd cell in each dimension (8x fewer cells, faster but less accurate)

    // Simulation file
    std::string filenameNGA = "/home/quirin/mlcfd/Repositories/jet/nga.case";
    std::string filenamePlic = "/home/quirin/mlcfd/Repositories/jet/plic.case";

    // Dataset
    //std::string dataset_path = "/home/quirin/mlcfd/Datasets/SixClasses/FirstMoment/s5_1M/data/data.bin";
    std::string dataset_path = "/home/quirin/mlcfd/Datasets/SixClasses/FirstMomentEigenv/s5_262k/data/data.bin";

    // Output directory for this whole experiment (distinct per call)
    // Example: stable_run_models/2026-02-25_153012/
    const auto now = std::chrono::system_clock::now();
    const std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d_%02d%02d%02d",
                  tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
                  tm.tm_hour, tm.tm_min, tm.tm_sec);

    fs::path experiment_dir = fs::path("stable_run_models") / buf;
    fs::create_directories(experiment_dir);

    // Store only predictions in RAM, not ml models or datasets, to save memory. Each entry is a vector of predicted classes for all interface cells.
    std::vector<std::vector<int>> predictions;
    predictions.reserve(no_runs);

    // Store model directories so we can load them later
    std::vector<fs::path> run_dirs;
    run_dirs.reserve(no_runs);

    // Train + classify multiple times
    for (int r = 0; r < no_runs; ++r) {
        std::cout << "\n=== Run " << r << " / " << no_runs << " ===\n";

        // Each run in its own distinct folder
        fs::path run_dir = experiment_dir / ("run_" + std::to_string(r));
        fs::create_directories(run_dir);
        run_dirs.push_back(run_dir);

        // Create a fresh classifier each run (only one dataset in RAM at a time)
        IRL::MLClassifier ml(stencil_size, input_size, hidden_size1, hidden_size2, hidden_size3, output_size);

        ml.updateDataParameters(
            no_batches,
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

        ml.loadDataset(dataset_path);
        //ml.canonicalize_data(canonicalize_symmetries);
        ml.preprocess_data(canonicalize_symmetries, noise_stddev);
        ml.updateTrainingParameters(learning_rate, batch_size, max_epochs, reduce_lr_patience, early_stop_patience);

        // In the future, could use seed setter
        // ml.setSeed(1234 + r);

        ml.trainModel();

        std::vector<int> savedClasses;
        //IRL::classify_simulation(ml, filenameNGA, filenamePlic, canonicalize_symmetries, include_Moments, include_Surface_Area, include_Eigenvalues, noise_stddev, 1e-10f, &savedClasses);
        IRL::classify_simulation(ml, filenameNGA, filenamePlic, canonicalize_symmetries, include_Moments, include_Surface_Area, include_Eigenvalues, noise_stddev, epsilon_connectivity, &savedClasses, downsample_factor);
    
        if (r > 0 && savedClasses.size() != predictions[0].size()) {
            throw std::runtime_error("Saved class vector size differs between runs!");
        }

        predictions.push_back(std::move(savedClasses));
        std::cout << "Classified interface cells: " << predictions.back().size() << "\n";

        // Save trained model into its folder.
        ml.saveModel(run_dir.string() + "/", false);

        // ml goes out of scope here => dataset RAM freed
    }

    // Print mean instability across ensemble
    const double mean_instability = compute_mean_instability(predictions);
    std::cout << "\nMean per-cell instability = " << mean_instability << "\n";

    // Pick the most agreeing model
    double most_agreeing_mean_agreement = 0.0;
    const int most_agreeing_run = pick_most_agreeing_model(predictions, &most_agreeing_mean_agreement);
    if (most_agreeing_run < 0) {
        throw std::runtime_error("No runs were executed; cannot select most agreeing model.");
    }

    fs::path most_agreeing_model_dir  = run_dirs[most_agreeing_run];
    fs::path most_agreeing_model_path = most_agreeing_model_dir / "ml_model.pt";

    // Save and print model selection file
    fs::path selection_file = experiment_dir / "model_selection.txt";
    {
        std::ofstream out(selection_file);
        if (!out) throw std::runtime_error("Failed to open " + selection_file.string());

        out << "simulation_file " << filenameNGA << "\n";
        out << "no_runs " << no_runs << "\n";
        out << "most_agreeing_model_run " << most_agreeing_run << "\n";
        out << "most_agreeing_model_path " << most_agreeing_model_path.string() << "\n";
        out << "most_agreeing_mean_agreement " << most_agreeing_mean_agreement << "\n";
        out << "mean_per_cell_instability " << mean_instability << "\n";
    }

    std::cout << "\n=== Model selection ===\n";
    std::cout << "Selection file: " << selection_file.string() << "\n";
    std::cout << "simulation_file: " << filenameNGA << "\n";
    std::cout << "no_runs: " << no_runs << "\n";
    std::cout << "most_agreeing_model_run: " << most_agreeing_run << "\n";
    std::cout << "most_agreeing_model_path: " << most_agreeing_model_path.string() << "\n";
    std::cout << "most_agreeing_mean_agreement: " << most_agreeing_mean_agreement << "\n";
    std::cout << "mean_per_cell_instability: " << mean_instability << "\n";

    // Reload selected model and classify simulation again
    {
        IRL::MLClassifier most_agreeing(stencil_size, input_size, hidden_size1, hidden_size2, hidden_size3, output_size);

        most_agreeing.loadModel(most_agreeing_model_path.string(), false);

        // Classify again using the most agreeing model
        //IRL::classify_simulation(most_agreeing, filenameNGA, filenamePlic, canonicalize_symmetries, include_Moments, include_Surface_Area, include_Eigenvalues, noise_stddev, 1e-10f, nullptr);
        IRL::classify_simulation(most_agreeing, filenameNGA, filenamePlic, canonicalize_symmetries, include_Moments, include_Surface_Area, include_Eigenvalues, noise_stddev, epsilon_connectivity, nullptr, downsample_factor);
    
    }

}
*/

void stable_classification() {
    int stencil_size = 5;

    // Data parameters
    int no_batches;
    int include_Moments = 0;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = false;  // true because dataset path says FirstMomentEigenv
    double paraboloid_coeff_stddev = 0.1;
    double hyperbolic_cylinder_opening_angle_stddev = 20; // degrees
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
    float epsilon_connectivity = 1e-12f;

    // Net Parameters
    int input_size =
        stencil_size * stencil_size * stencil_size
        * (include_Moments >= 1 ? (include_Surface_Area ? 5 : 4) : 1)
        + (include_Moments >= 2 ? 6 : 0)
        + (include_Eigenvalues ? 3 : 0);

    int hidden_size1 = 256;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 6;

    // Training parameters
    double learning_rate = 0.001;
    int batch_size = 64;
    int max_epochs = 50;
    int reduce_lr_patience = 4;
    int early_stop_patience = 8;

    // Classification parameters
    int canonicalize_symmetries = 48;
    const int no_runs = 10;
    float noise_stddev = 0.0f;
    //int downsample_factor = 2;

    // Simulations used for stable model selection
    struct SimulationFiles {
        std::string name;
        std::string filenameNGA;
        std::string filenamePlic;
        int downsample_factor;
    };

    std::vector<SimulationFiles> simulations = {
        {
            "Jet",
            "/home/quirin/mlcfd/Repositories/jetvtr/data.000031.vtr",
            "/home/quirin/mlcfd/Repositories/jetvtr/plic.000031.vtu",
            2   // downsample factor for sim_1
        },
        {
            "Bag",
            "/home/quirin/mlcfd/Repositories/bagvtr/data.000001.vtr",
            "/home/quirin/mlcfd/Repositories/jetvtr/plic.000031.vtu",
            1   // downsample factor for sim_2
        }
    };

    // Dataset
    std::string dataset_path =
        "/home/quirin/mlcfd/Datasets/SixClasses/Thesis2/ZerothMoments/data/data.bin";

    // Output directory for this whole experiment
    const auto now = std::chrono::system_clock::now();
    const std::time_t t = std::chrono::system_clock::to_time_t(now);

    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif

    char buf[64];
    std::snprintf(
        buf,
        sizeof(buf),
        "%04d-%02d-%02d_%02d%02d%02d",
        tm.tm_year + 1900,
        tm.tm_mon + 1,
        tm.tm_mday,
        tm.tm_hour,
        tm.tm_min,
        tm.tm_sec
    );

    fs::path experiment_dir = fs::path("stable_run_models") / buf;
    fs::create_directories(experiment_dir);

    // Store only predictions in RAM, not models or datasets.
    // Each entry corresponds to one trained model.
    // The vector contains the concatenated predictions from both simulations.
    std::vector<std::vector<int>> predictions;
    predictions.reserve(no_runs);

    // Store model directories so we can load the selected model later
    std::vector<fs::path> run_dirs;
    run_dirs.reserve(no_runs);

    // Reference number of classified cells for each simulation.
    // This ensures that every model classified the same cells in the same simulation.
    std::vector<std::size_t> reference_cell_counts;

    // Train + classify multiple times
    for (int r = 0; r < no_runs; ++r) {
        std::cout << "\n=== Run " << r << " / " << no_runs << " ===\n";

        // Each run in its own folder
        fs::path run_dir = experiment_dir / ("run_" + std::to_string(r));
        fs::create_directories(run_dir);
        run_dirs.push_back(run_dir);

        // Create a fresh classifier each run
        IRL::MLClassifier ml(
            stencil_size,
            input_size,
            hidden_size1,
            hidden_size2,
            hidden_size3,
            output_size
        );

        ml.updateDataParameters(
            no_batches,
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

        ml.loadDataset(dataset_path);

        // ml.canonicalize_data(canonicalize_symmetries);
        ml.preprocess_data(canonicalize_symmetries, noise_stddev);

        ml.updateTrainingParameters(
            learning_rate,
            batch_size,
            max_epochs,
            reduce_lr_patience,
            early_stop_patience
        );

        // In the future, could use seed setter
        // ml.setSeed(1234 + r);

        ml.trainModel();

        // Classify both simulations and concatenate the class vectors.
        // This combined vector is what is used for model agreement.
        std::vector<int> savedClassesCombined;
        std::vector<std::size_t> current_cell_counts;

        for (const auto& sim : simulations) {
            std::vector<int> savedClassesThisSimulation;

            std::cout << "\nClassifying " << sim.name << "...\n";

            IRL::classify_simulation(
                ml,
                sim.filenameNGA,
                sim.filenamePlic,
                canonicalize_symmetries,
                include_Moments,
                include_Surface_Area,
                include_Eigenvalues,
                noise_stddev,
                epsilon_connectivity,
                &savedClassesThisSimulation,
                sim.downsample_factor
            );

            current_cell_counts.push_back(savedClassesThisSimulation.size());

            std::cout << sim.name
                      << " classified interface cells: "
                      << savedClassesThisSimulation.size()
                      << "\n";

            savedClassesCombined.insert(
                savedClassesCombined.end(),
                savedClassesThisSimulation.begin(),
                savedClassesThisSimulation.end()
            );
        }

        // Store reference classified-cell counts from the first run
        if (r == 0) {
            reference_cell_counts = current_cell_counts;
        } else {
            if (current_cell_counts.size() != reference_cell_counts.size()) {
                throw std::runtime_error(
                    "Number of classified simulations differs between runs!"
                );
            }

            for (std::size_t s = 0; s < current_cell_counts.size(); ++s) {
                if (current_cell_counts[s] != reference_cell_counts[s]) {
                    throw std::runtime_error(
                        "Classified cell count differs between runs for simulation " +
                        simulations[s].name
                    );
                }
            }
        }

        if (r > 0 && savedClassesCombined.size() != predictions[0].size()) {
            throw std::runtime_error(
                "Combined saved class vector size differs between runs!"
            );
        }

        predictions.push_back(std::move(savedClassesCombined));

        std::cout << "Combined classified interface cells: "
                  << predictions.back().size()
                  << "\n";

        // Save trained model into its folder
        ml.saveModel(run_dir.string() + "/", false);

        // ml goes out of scope here, so dataset RAM is freed
    }

    // Print mean instability across ensemble
    const double mean_instability = compute_mean_instability(predictions);

    std::cout << "\nMean per-cell instability = "
              << mean_instability
              << "\n";

    // Pick the most agreeing model based on both simulations combined
    double most_agreeing_mean_agreement = 0.0;

    const int most_agreeing_run =
        pick_most_agreeing_model(predictions, &most_agreeing_mean_agreement);

    if (most_agreeing_run < 0) {
        throw std::runtime_error(
            "No runs were executed; cannot select most agreeing model."
        );
    }

    fs::path most_agreeing_model_dir = run_dirs[most_agreeing_run];
    fs::path most_agreeing_model_path = most_agreeing_model_dir / "ml_model.pt";

    // Save model selection file
    fs::path selection_file = experiment_dir / "model_selection.txt";

    {
        std::ofstream out(selection_file);

        if (!out) {
            throw std::runtime_error(
                "Failed to open " + selection_file.string()
            );
        }

        out << "dataset_path " << dataset_path << "\n";
        out << "no_runs " << no_runs << "\n";
        
        out << "canonicalize_symmetries " << canonicalize_symmetries << "\n";
        out << "include_Moments " << include_Moments << "\n";
        out << "include_Surface_Area " << include_Surface_Area << "\n";
        out << "include_Eigenvalues " << include_Eigenvalues << "\n";

        out << "simulations_used_for_selection "
            << simulations.size()
            << "\n";

        for (std::size_t s = 0; s < simulations.size(); ++s) {
            out << "simulation_" << s << "_name "
                << simulations[s].name
                << "\n";

            out << "simulation_" << s << "_nga "
                << simulations[s].filenameNGA
                << "\n";

            out << "simulation_" << s << "_downsample_factor "
                << simulations[s].downsample_factor
                << "\n";

            out << "simulation_" << s << "_plic "
                << simulations[s].filenamePlic
                << "\n";

            out << "simulation_" << s << "_classified_cells "
                << reference_cell_counts[s]
                << "\n";
        }

        out << "most_agreeing_model_run "
            << most_agreeing_run
            << "\n";

        out << "most_agreeing_model_path "
            << most_agreeing_model_path.string()
            << "\n";

        out << "most_agreeing_mean_agreement "
            << most_agreeing_mean_agreement
            << "\n";

        out << "mean_per_cell_instability "
            << mean_instability
            << "\n";
    }

    std::cout << "\n=== Model selection ===\n";
    std::cout << "Selection file: "
              << selection_file.string()
              << "\n";

    std::cout << "dataset_path: "
              << dataset_path
              << "\n";

    std::cout << "no_runs: "
              << no_runs
              << "\n";

    std::cout << "simulations_used_for_selection: "
              << simulations.size()
              << "\n";

    for (std::size_t s = 0; s < simulations.size(); ++s) {
        std::cout << "simulation_" << s << "_name: "
                  << simulations[s].name
                  << "\n";

        std::cout << "simulation_" << s << "_nga: "
                  << simulations[s].filenameNGA
                  << "\n";

        std::cout << "simulation_" << s << "_plic: "
                  << simulations[s].filenamePlic
                  << "\n";

        std::cout << "simulation_" << s << "_classified_cells: "
                  << reference_cell_counts[s]
                  << "\n";
    }

    std::cout << "most_agreeing_model_run: "
              << most_agreeing_run
              << "\n";

    std::cout << "most_agreeing_model_path: "
              << most_agreeing_model_path.string()
              << "\n";

    std::cout << "most_agreeing_mean_agreement: "
              << most_agreeing_mean_agreement
              << "\n";

    std::cout << "mean_per_cell_instability: "
              << mean_instability
              << "\n";

    // Reload selected model and classify both simulations again
    {
        IRL::MLClassifier most_agreeing(
            stencil_size,
            input_size,
            hidden_size1,
            hidden_size2,
            hidden_size3,
            output_size
        );

        most_agreeing.loadModel(most_agreeing_model_path.string(), false);

        for (const auto& sim : simulations) {
            std::cout << "\nReclassifying "
                      << sim.name
                      << " with most agreeing model...\n";

            IRL::classify_simulation(
                most_agreeing,
                sim.filenameNGA,
                sim.filenamePlic,
                canonicalize_symmetries,
                include_Moments,
                include_Surface_Area,
                include_Eigenvalues,
                noise_stddev,
                epsilon_connectivity,
                nullptr,
                sim.downsample_factor
            );
        }
    }
}

int main (int argc, char* argv[]) {
    
    int stencil_size = 5;

    //Data parameters
    int no_batches = 4096 * 4 * 2;
    int include_Moments = 1;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = false;
    double paraboloid_coeff_stddev = 0.1;
    double hyperbolic_cylinder_opening_angle_stddev = 20; //degrees
    double sheet_coeff_stddev = 0.1;
    double sheet_thickness_stddev = 0.0;
    double cylinder_radius_stddev = 0.0;
    double radius_circle_min = 2.5;
    double radius_circle_max = 10.0;
    double sphere_radius_stddev = 0.0;
    double ellipsoid_subgrid_stddev = 0.7;
    double min_long_ellipsoid_axis = 3.0;
    double max_long_ellipsoid_axis = 5.0;
    bool exact_2nd_moment = false;  // enable calculation of exact 2nd moments for data generation
    bool visualize = false; // if true, print centroids and / or write surfaces
    double machineZero = 1e-12;
    double lower_limit_subgrid = machineZero;
    double upper_limit_subgrid = std::sqrt(3.0);
    double class0_max_characteristic = 2.5;
    float epsilon_connectivity = 1e-12f;

    // Net Parameters
    int input_size = stencil_size * stencil_size * stencil_size 
    * (include_Moments >= 1 ? /*(include_Surface_Area ? 5 : 4) ---> remove the 4 again if uncommenting ->*/ 4 : 1)  // 4 if include_Moments >= 1 because we have vfrac + (mx,my,mz) per cell, otherwise just vfrac
    + (include_Moments >= 2 ? 6 : 0)  // +6 if include_Mome´nts >= 2 because we have (xx, yy, zz, xy, xz, yz) components of the 2nd moment tensor; otherwise none
    + (include_Eigenvalues ? 3 : 0);
    int hidden_size1 = 256;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 6; //CHANGED 4 to 6

    //Training parameters
    double learning_rate = 0.001; //was 0.01 for SGD optimizer
    int batch_size = 64;
    int max_epochs = 50;
    int reduce_lr_patience = 4;
    int early_stop_patience = 8;

    IRL::InertiaClassifier inertia_classifier(stencil_size, 1, 0.85, 1.5);
    //IRL::MLClassifier_E3NN ml(stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    IRL::MLClassifier ml(stencil_size, input_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    ml.updateDataParameters(
            no_batches,
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
    ml.updateTrainingParameters(learning_rate, batch_size, max_epochs, reduce_lr_patience, early_stop_patience);
    //ml.generateDataset();
    //ml.loadDataset("/home/quirin/mlcfd/Datasets/SixClasses/Thesis2/FirstMoments/data/data.bin");
    //ml.loadDataset("/home/quirin/mlcfd/Datasets/SixClasses/Transition_new/s5_262k_new/ZerothMoments/data/data.bin");

    //ml.appendDataset("/home/quirin/mlcfd/Datasets/SixClasses/Transition/s5_262k/data/data.bin", false);
    //ml.saveDataset("data");
    /*//Below: Used to remove 1st moments from dataset/
    include_Moments = 0;
    ml.updateDataParameters(
            no_batches,
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
        */

    int canonicalize_symmetries = 48;
    float noise_stddev = 0.0f;
    //ml.preprocess_data(canonicalize_symmetries, noise_stddev);
    //ml.saveDataset("data");

    //ml.checkStatesForNaNOrInf();
    
    //ml.trainModel();
    //ml.outputTrainingResults();
    //ml.saveModel("model/");
    //ml.loadModel("/home/quirin/mlcfd/Datasets/SixClasses/Thesis2/ZerothMoments/stable_run_models/2026-07-01_133021/run_8/ml_model.pt"); //Thesis zeroth most agreeing
    ml.loadModel("/home/quirin/mlcfd/Datasets/SixClasses/Thesis2/FirstMoments/stable_run_models/2026-07-01_014127/run_4/ml_model.pt"); // Thesis first moments most agreeing
    //ml.loadModel("/home/quirin/mlcfd/Datasets/SixClasses/Thesis2/FirstMomentsEigv/stable_run_models/2026-06-30_203023/run_0/ml_model.pt"); // Thesis first moments eigenvalues most agreeing
    ml.exportRuntimeWeightsAndBiasesHeader();

    //IRL::MLClassifierNoTorch ml_no_torch(stencil_size);
    //ml_no_torch.loadWeights("/home/quirin/mlcfd/Datasets/SixClasses/ZerothMoment/model/runtime_weights.bin");

    // vtk reader
    std::string filenameNGA = "/home/quirin/mlcfd/Repositories/jetvtr/data.000031.vtr";
    std::string filenamePlic = "/home/quirin/mlcfd/Repositories/jetvtr/plic.000031.vtu";
    //std::string filenameNGA = "/home/quirin/mlcfd/Repositories/bagvtr/data.000001.vtr";
    //std::string filenamePlic = "/home/quirin/mlcfd/Repositories/bagvtr/plic.000001.vtu";
    int downsample_factor = 2;
    //IRL::classify_simulation(ml, filenameNGA, filenamePlic, canonicalize_symmetries, include_Moments, include_Surface_Area, include_Eigenvalues, noise_stddev, epsilon_connectivity, nullptr, downsample_factor);
    
    //stable_classification();

    //IRL::Data_gen gen;

    //gen.generateState(2,5,1,false,0.1,0.1,0.5,0.0,0.5,0.0,0.5,0.0,true);

    //find_dataset_size();
    //shell_testcase(ml);
    return 0;
}