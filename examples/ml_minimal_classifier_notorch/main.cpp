#include "irl/ml_classification/stencil_rotator.h"
#include "irl/ml_classification/ml_classifier_notorch.h"
#include "irl/ml_classification/inertia_calc.h"

#include <vector>
#include <iostream>

int main (int argc, char* argv[]) {
    
    int stencil_size = 5;
    int include_Moments = 1;
    bool include_Surface_Area = false;
    bool include_Eigenvalues = true;
    int canonicalize_symmetries = 48;
    float noise_stddev = 0.0f;
    float epsilon_connectivity = 1e-12f;

    // Allocate stencil containers
    std::vector<std::vector<std::vector<double>>> vfrac(
        stencil_size,
        std::vector<std::vector<double>>(
            stencil_size,
            std::vector<double>(stencil_size, 0.0)));

    std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoment(
        stencil_size,
        std::vector<std::vector<Eigen::Vector3d>>(
            stencil_size,
            std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())));

    std::vector<std::vector<std::vector<double>>> surfaceArea(
                    stencil_size,
                    std::vector<std::vector<double>>(
                        stencil_size,
                        std::vector<double>(stencil_size, 0.0)));

    IRL::MLClassifierNoTorch classifier(stencil_size);


    // Flatten stencil into 1D vector
    std::vector<double> flattened_state;
    for (int si = 0; si < stencil_size; ++si) {
        for (int sj = 0; sj < stencil_size; ++sj) {
            for (int sk = 0; sk < stencil_size; ++sk) {
                if (include_Moments >= 0) {
                    flattened_state.push_back(vfrac[si][sj][sk]);
                }
                if (include_Moments >= 1) {
                    flattened_state.push_back(firstMoment[si][sj][sk].x());
                    flattened_state.push_back(firstMoment[si][sj][sk].y());
                    flattened_state.push_back(firstMoment[si][sj][sk].z());
                }
                if (include_Surface_Area) {
                    flattened_state.push_back(surfaceArea[si][sj][sk]);
                }
            }
        }
    }

    // Make flattened_state a float vector
    std::vector<float> flattened_state_float(flattened_state.begin(), flattened_state.end());

    if (include_Moments >= 2) {
        Eigen::Matrix3d secondMoment = IRL::compute2ndMoment(flattened_state_float, stencil_size, 1); // from_ith_moment=1
        flattened_state_float.push_back(static_cast<float>(secondMoment(0, 0))); // xx
        flattened_state_float.push_back(static_cast<float>(secondMoment(1, 1))); // yy
        flattened_state_float.push_back(static_cast<float>(secondMoment(2, 2))); // zz
        flattened_state_float.push_back(static_cast<float>(secondMoment(0, 1))); // xy
        flattened_state_float.push_back(static_cast<float>(secondMoment(0, 2))); // xz
        flattened_state_float.push_back(static_cast<float>(secondMoment(1, 2))); // yz
    }

    if (include_Eigenvalues) {
        IRL::appendInertiaEigenvalues(flattened_state_float, stencil_size, 1);
    }

    //Preprocess stencil
    IRL::preprocess_stencil(flattened_state_float,
            stencil_size,
            canonicalize_symmetries,
            include_Moments,
            include_Surface_Area,
            include_Eigenvalues,
            noise_stddev,
            epsilon_connectivity);

    // Classify
    std::vector<float> out_probs;
    int predicted_class = classifier.classify(flattened_state_float, &out_probs);
    float max_prob = *std::max_element(out_probs.begin(), out_probs.end());

    std::cout << "Predicted class: " << predicted_class << "\n";
    
    return 0;
}