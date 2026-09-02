#include <vtkCellCenters.h>
#include "irl/ml_classification/ml_classifier.h"
#include "irl/ml_classification/vtk_in.h"
#include "irl/ml_classification/inertia_classifier.h"
#include "irl/ml_classification/data_gen.h"

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

            ml.preprocess_data();
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

int main (int argc, char* argv[]) {
    
    int stencil_size = 5;

    //Data parameters
    int no_batches = 4096 * 4;
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
    int max_epochs = 10;
    int reduce_lr_patience = 4;
    int early_stop_patience = 8;

    int canonicalize_symmetries = 48;
    float noise_stddev = 0.0f;

    IRL::InertiaClassifier inertia_classifier(stencil_size, 1, 0.85, 1.5);
    //IRL::MLClassifier_E3NN ml(stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    IRL::MLClassifier ml(stencil_size, epsilon_connectivity, canonicalize_symmetries, noise_stddev, include_Moments, include_Surface_Area, include_Eigenvalues, hidden_size1, hidden_size2, hidden_size3, output_size);

    ml.updateDataParameters(
        no_batches,
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
    //ml.loadDataset("/home/quirin/mlcfd/Datasets/SixClasses/Thesis2/ZerothMoments/data/data.bin");
    //ml.loadDataset("/home/quirin/mlcfd/Datasets/SixClasses/NoEllipsoidLigTips/s5_2M/data/data.bin");
    //ml.appendDataset("/home/quirin/mlcfd/Datasets/SixClasses/NoEllipsoidLigTips/s5_1M/data/data.bin", false);
    //ml.retain_only_0th_moments();
    //ml.saveDataset("data");


    //ml.preprocess_data();
    //ml.saveDataset("data");

    //ml.checkStatesForNaNOrInf();
    
    //ml.trainModel();
    //ml.outputTrainingResults();
    //ml.saveModel("model/");
    ml.loadModel("/home/quirin/mlcfd/Datasets/SixClasses/NoEllipsoidLigTips/s5_2M/stable_run_models/2026-08-29_032729/run_9/ml_model.pt"); 
    //ml.loadModel("/home/quirin/mlcfd/Datasets/SixClasses/Thesis2/ZerothMoments/stable_run_models/2026-07-15_164648/run_2/ml_model.pt"); //Thesis zeroth most agreeing
    //ml.loadModel("/home/quirin/mlcfd/Datasets/SixClasses/Thesis2/FirstMoments/stable_run_models/2026-07-01_014127/run_4/ml_model.pt"); // Thesis first moments most agreeing
    //ml.loadModel("/home/quirin/mlcfd/Datasets/SixClasses/Thesis2/FirstMomentsEigv/stable_run_models/2026-06-30_203023/run_0/ml_model.pt"); // Thesis first moments eigenvalues most agreeing
    //ml.exportRuntimeWeightsAndBiasesHeader();

    //IRL::MLClassifierNoTorch ml_no_torch(stencil_size);
    //ml_no_torch.loadWeights("/home/quirin/mlcfd/Datasets/SixClasses/ZerothMoment/model/runtime_weights.bin");

    // vtk reader
    //std::string filenameNGA = "/home/quirin/mlcfd/Repositories/jetvtr/data.000031.vtr";
    //std::string filenamePlic = "/home/quirin/mlcfd/Repositories/jetvtr/plic.000031.vtu";
    // std::string filenameNGA = "/home/quirin/mlcfd/Repositories/bagvtr/data.000001.vtr";
    // std::string filenamePlic = "/home/quirin/mlcfd/Repositories/bagvtr/plic.000001.vtu";
    // int downsample_factor = 1;
    // double pdistribution_step = 0.0;
    // IRL::classify_simulation(ml, filenameNGA, filenamePlic, canonicalize_symmetries, include_Moments, include_Surface_Area, include_Eigenvalues, noise_stddev, epsilon_connectivity, nullptr, downsample_factor, pdistribution_step);
    
    //std::string dataDirectory = "/home/quirin/mlcfd/Repositories/bag-breakup/data";
    //std::string plicDirectory = "/home/quirin/mlcfd/Repositories/bag-breakup/plic";
    std::string dataDirectory = "/home/quirin/mlcfd/Repositories/round-jet/data";
    std::string plicDirectory = "/home/quirin/mlcfd/Repositories/round-jet/plic";
    std::string outputDirectory = ".";

    int downsample_factor = 2;
    double pdistribution_step = 0.0;
    
    IRL::classify_simulation(
        ml,
        dataDirectory,
        plicDirectory,
        outputDirectory,
        nullptr,
        downsample_factor,
        pdistribution_step);

    return 0;
}