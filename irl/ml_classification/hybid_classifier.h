#pragma once
#include "ml_classifier.h"
#include "inertia_calc.h"
#include <Eigen/Dense>
#include <algorithm>
#include <torch/torch.h>
#include <iostream>

namespace IRL {

class HybridClassifier : public MLClassifier {
private:
    int from_ith_moment = 0;
public:
    // Reuse MLClassifier’s constructor
    using MLClassifier::MLClassifier;

    // Override classify()
    int classify(const std::vector<float>& flattened_state,
                 std::vector<float>* out_probs = nullptr) override {
        
        // Convert flattened_state to double vector
        std::vector<double> flattened_state_double(flattened_state.begin(), flattened_state.end());
        
        // Compute inertia tensor
        Eigen::Matrix3d I = computeInertiaTensor(flattened_state_double,
                                                 stencil_size,
                                                 from_ith_moment);

        // Get eigenvalues
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(I);
        Eigen::Vector3d evals = solver.eigenvalues();

        // Sort eigenvalues descending: I1 >= I2 >= I3
        std::sort(evals.data(), evals.data() + 3, std::greater<double>());
        double I1 = evals[0], I2 = evals[1], I3 = evals[2];
        std::vector<double> evals_vec = {I1, I2, I3};
        
        torch::NoGradGuard no_grad;
        auto input = torch::tensor(evals_vec,
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