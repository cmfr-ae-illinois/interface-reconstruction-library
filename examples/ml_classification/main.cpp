#include <vtkCellCenters.h>
//#include "testcases.h"
#include "irl/ml_classification/ml_classifier.h"
#include "irl/ml_classification/vtk_in.h"
#include "irl/ml_classification/inertia_classifier.h"
#include "irl/ml_classification/hybid_classifier.h"
#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/ml_classifier_e3nn.h"

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

    // Data parameters
    int include_Moments = 1;
    bool include_Surface_Area = false;
    double paraboloid_coeff_stddev = 0.1;
    double sheet_coeff_stddev = 0.1;
    double max_sheet_thickness = 1.0;
    double sheet_thickness_stddev = 0.0;
    double max_cylinder_radius = 0.5;
    double cylinder_radius_stddev = 0.0;
    bool include_truncated_cylinder = true;
    double max_sphere_radius = 0.5;
    double sphere_radius_stddev = 0.0;
    bool exact_2nd_moment = true;

    // Net parameters
    int input_size = stencil_size * stencil_size * stencil_size
        * (include_Moments >= 1 ? 4 : 1)
        + (include_Moments >= 2 ? 6 : 0)
        + (include_Moments >= 3 ? 3 : 0);

    int hidden_size1 = 256;
    int hidden_size2 = 128;
    int hidden_size3 = 64;
    int output_size = 6;

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
    const int batch_increment = 1024;
    const int max_steps = 20;
    const double min_improvement = 0.005; // absolute accuracy improvement

    const std::string base_dir =
        "/home/quirin/mlcfd/Datasets/SixClasses/FirstMoment/DatasetSizeStudy/";
    fs::create_directories(base_dir);

    std::vector<double> accuracies;
    std::vector<int> batch_counts;
    std::vector<long long> sample_counts;

    std::string previous_dataset_file;

    for (int step = 1; step <= max_steps; ++step) {
        const int total_batches = step * batch_increment;
        const long long total_samples =
            static_cast<long long>(total_batches) * batch_size;

        const std::string dataset_name =
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
                batch_increment, include_Moments, include_Surface_Area,
                paraboloid_coeff_stddev,
                sheet_coeff_stddev,
                max_sheet_thickness, sheet_thickness_stddev,
                max_cylinder_radius, cylinder_radius_stddev, include_truncated_cylinder,
                max_sphere_radius, sphere_radius_stddev,
                exact_2nd_moment
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

            std::cout << "Final test accuracy: " << final_acc << "\n";
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

            if (improvement < min_improvement) {
                std::cout << "Stopping criterion reached.\n";
                break;
            }
        }
    }

    std::cout << "\n=== Dataset size study finished ===\n";
    for (size_t i = 0; i < accuracies.size(); ++i) {
        std::cout << "Batches: " << batch_counts[i]
                  << ", Samples: " << sample_counts[i]
                  << ", Final test accuracy: " << accuracies[i] << "\n";
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

void stable_classification() {
    int stencil_size = 5;

    //Data parameters
    int no_batches = 4096;
    int include_Moments = 1;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = true;
    double paraboloid_coeff_stddev = 0.1;
    double sheet_coeff_stddev = 0.1;
    double max_sheet_thickness = 1.0;
    double sheet_thickness_stddev = 0.0;
    double max_cylinder_radius = 0.5;
    double cylinder_radius_stddev = 0.0;
    double max_sphere_radius = 0.5;
    double sphere_radius_stddev = 0.0;
    bool exact_2nd_moment = false;  // enable calculation of exact 2nd moments for data generation
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

        ml.updateDataParameters(no_batches, include_Moments, include_Surface_Area, include_Eigenvalues,
                            paraboloid_coeff_stddev,
                            sheet_coeff_stddev,
                            max_sheet_thickness, sheet_thickness_stddev,
                            max_cylinder_radius, cylinder_radius_stddev,
                            max_sphere_radius, sphere_radius_stddev,
                            exact_2nd_moment);           

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

int main (int argc, char* argv[]) {
    
    int stencil_size = 5;

    //Data parameters
    int no_batches = 4096;
    int include_Moments = 1;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = true;
    double paraboloid_coeff_stddev = 0.1;
    double sheet_coeff_stddev = 0.1;
    double max_sheet_thickness = 1.0;
    double sheet_thickness_stddev = 0.0;
    double max_cylinder_radius = 0.5;
    double cylinder_radius_stddev = 0.0;
    double max_sphere_radius = 0.5;
    double sphere_radius_stddev = 0.0;
    bool exact_2nd_moment = false;  // enable calculation of exact 2nd moments for data generation
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

    IRL::InertiaClassifier inertia_classifier(stencil_size, 1, 0.85, 1.5);
    //IRL::MLClassifier_E3NN ml(stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    IRL::MLClassifier ml(stencil_size, input_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    
    ml.updateDataParameters(no_batches, include_Moments, include_Surface_Area, include_Eigenvalues,
                            paraboloid_coeff_stddev,
                            sheet_coeff_stddev,
                            max_sheet_thickness, sheet_thickness_stddev,
                            max_cylinder_radius, cylinder_radius_stddev,
                            max_sphere_radius, sphere_radius_stddev,
                            exact_2nd_moment);                    
    
    //ml.generateDataset();
    //ml.loadDataset("/home/quirin/mlcfd/Datasets/SixClasses/SecondApproxEigenvSurfaces/s5_524k/data/data.bin");
    //ml.loadDataset("/home/quirin/mlcfd/Datasets/SixClasses/SecondApproxEigenvSurfaces/s5_524k/data/data.bin");

    //ml.appendDataset("/home/quirin/mlcfd/Datasets/float/From1/s5_2M/data/data.bin", false);
    //ml.saveDataset("data");
    int canonicalize_symmetries = 48;
    float noise_stddev = 0.0f;
    //ml.preprocess_data(canonicalize_symmetries, noise_stddev);

    //ml.checkStatesForNaNOrInf();

    ml.updateTrainingParameters(learning_rate, batch_size, max_epochs, reduce_lr_patience, early_stop_patience);
    //ml.trainModel();
    //ml.outputTrainingResults();
    //ml.saveModel("model/");
    //ml.loadModel("/home/quirin/mlcfd/Datasets/SixClasses/FirstMomentEigenv/s5_262k/model/ml_model.pt");
    ml.loadModel("/home/quirin/mlcfd/Datasets/SixClasses/FirstMomentEigenv/s5_262k/stable_run_models/2026-05-09_110350/most agreeing/ml_model.pt");

    // vtk reader
    std::string filenameNGA = "/home/quirin/mlcfd/Repositories/bagvtr/data.000005new.vtr";
    std::string filenamePlic = "/home/quirin/mlcfd/Repositories/bagvtr/plic.000005.vtu";
    int downsample_factor = 1;
    IRL::classify_simulation(ml, filenameNGA, filenamePlic, canonicalize_symmetries, include_Moments, include_Surface_Area, include_Eigenvalues, noise_stddev, epsilon_connectivity, nullptr, downsample_factor);
    
    //stable_classification();

    //IRL::Data_gen gen;

    //gen.generateState(2,5,1,false,0.1,0.1,0.5,0.0,0.5,0.0,0.5,0.0,true);

    //find_dataset_size();
    //shell_testcase(ml);
    return 0;
}