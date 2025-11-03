#pragma once
#include "classifier.h"
#include "inertia_calc.h" // for computeInertiaTensor
#include <Eigen/Dense>
#include <algorithm>

namespace IRL {

class InertiaClassifier : public Classifier {
private:
    int from_ith_moment;
    double sphere_tol;
    double gap;

public:
    explicit InertiaClassifier(int stencil,
                               int from_ith_moment_ = 1,
                               double sphere_tol_ = 0.85,
                               double gap_ = 1.5)
        : Classifier(stencil),
          from_ith_moment(from_ith_moment_),
          sphere_tol(sphere_tol_),
          gap(gap_) {}

    int classify(const std::vector<double>& flattened_state, std::vector<float>* out_probs = nullptr) override {
        // Compute inertia tensor

        Eigen::Matrix3d I = computeInertiaTensor(flattened_state,
                                                 stencil_size,
                                                 from_ith_moment);

        // Get eigenvalues
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(I);
        Eigen::Vector3d evals = solver.eigenvalues();

        // Sort eigenvalues descending: I1 >= I2 >= I3
        std::sort(evals.data(), evals.data() + 3, std::greater<double>());
        double I1 = evals[0], I2 = evals[1], I3 = evals[2];

        if (out_probs) {
            out_probs->resize(4, 1.0f);
        }

        // Avoid divide-by-zero
        if (I1 <= 1e-14) return 0;

        // --- Classification rules ---
        if (I1 > gap * I2 && I2 > I3) {
            return 3; // sheet
        }
        else if (I2 > gap * I3) {
            return 1; // ligament
        }
        else if ((I3 / I1) > sphere_tol) {
            return 2; // sphere
        }
        else {
            return 0; // paraboloid
        }
    }
};

} // namespace IRL