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

int main (int argc, char* argv[]) {
    
    // Classifier parameters
    int stencil_size = 5;
    int canonicalize_symmetries = 48;
    float noise_stddev = 0.0f;
    int include_Moments = 1;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = false;
    float epsilon_connectivity = 1e-12f;

    // Preprocessing parameters
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

    const std::string base_dir = "/home/quirin/mlcfd/Datasets/SixClasses/Thesis/DatasetSizeStudy/";
    
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
            IRL::MLClassifier ml(stencil_size, epsilon_connectivity, canonicalize_symmetries, noise_stddev, include_Moments, include_Surface_Area, include_Eigenvalues, hidden_size1, hidden_size2, hidden_size3, output_size);

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
    
    return 0;
}