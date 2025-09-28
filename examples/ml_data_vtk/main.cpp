#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <random>
#include <Eigen/Dense>
#include <vector>

#include "irl/ml_classification/data_gen.h"


// Compute inertia tensor from stencil data with first moments
Eigen::Matrix3d computeInertiaTensor(const std::vector<double>& flattened_state,
                                     bool include_firstMoment,
                                     int stencil_size,
                                     double cell_volume = 1.0)
{
    struct Cell {
        double alpha;
        Eigen::Vector3d firstMoment; // stores V * centroid
    };

    // Rebuild stencil cells from flattened_state
    std::vector<Cell> cells;
    size_t idx = 0;
    for (int i = 0; i < stencil_size; ++i) {
        for (int j = 0; j < stencil_size; ++j) {
            for (int k = 0; k < stencil_size; ++k) {
                Cell c;
                c.alpha = flattened_state[idx++];
                if (include_firstMoment) {
                    double mx = flattened_state[idx++];
                    double my = flattened_state[idx++];
                    double mz = flattened_state[idx++];
                    c.firstMoment = Eigen::Vector3d(mx, my, mz);
                } else {
                    c.firstMoment = Eigen::Vector3d::Zero();
                }
                cells.push_back(c);
            }
        }
    }

    // Compute neighborhood liquid centroid from first moments
    double m_total = 0.0;
    Eigen::Vector3d c_liq(0,0,0);
    for (auto& c : cells) {
        double vol = c.alpha * cell_volume;
        m_total += vol;
        c_liq += c.firstMoment;
    }
    if (m_total > 1e-14) {
        c_liq /= m_total;
    }

    // Assemble inertia tensor about c_liq
    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    for (auto& c : cells) {
        double vol = c.alpha * cell_volume;
        if (vol < 1e-14) continue;

        // Get the subcell centroid back from first moment
        Eigen::Vector3d centroid = c.firstMoment / vol;

        Eigen::Vector3d r = centroid - c_liq;
        double x = r[0], y = r[1], z = r[2];

        I(0,0) += vol * (y*y + z*z);
        I(1,1) += vol * (x*x + z*z);
        I(2,2) += vol * (x*x + y*y);
        I(0,1) -= vol * x * y;
        I(1,0) -= vol * x * y;
        I(0,2) -= vol * x * z;
        I(2,0) -= vol * x * z;
        I(1,2) -= vol * y * z;
        I(2,1) -= vol * y * z;
    }

    return I;
}

int classifyCellViaInertia(const Eigen::Vector3d& evals_unsorted, double gap = 1.5, double sphere_tol = 0.85) {
    // Sort eigenvalues in descending order: I1 >= I2 >= I3
    Eigen::Vector3d evals = evals_unsorted;
    std::sort(evals.data(), evals.data()+3, std::greater<double>());
    double I1 = evals[0], I2 = evals[1], I3 = evals[2];

    // Avoid divide-by-zero
    if (I1 <= 1e-14) return 0;

    // Ratios
    double R21 = I2 / I1;
    double R31 = I3 / I1;

    // --- Classification rules ---
    // Sheet: gap between I1 and I2
    if (I1 > gap * I2 && I2 > I3) {
        return 3; // sheet
    }
    // Ligament: gap between I2 and I3
    else if (I2 > gap * I3) {
        return 1; // ligament
    }
    // Sphere: nearly isotropic
    else if (R21 > sphere_tol && R31 > sphere_tol) {
        return 2; // sphere
    }
    // Else: paraboloid / well-resolved interface
    else {
        return 0; // paraboloid
    }
}



int main(int argc, char* argv[]) {
    int stencil_size = 3;
    IRL::Data_gen data_gen;

    // Map integers to class names
    std::map<int, std::string> classNames = {
        {0, "paraboloid"},
        {1, "ligament"},
        {2, "sphere"},
        {3, "sheet"}
    };

    // ---- Part 1: Single example per class with visualize ----
    for (int trueClass = 0; trueClass < 4; ++trueClass) {
        std::vector<double> flattened_state =
            data_gen.generate_State(trueClass, stencil_size, true, true);

        Eigen::Matrix3d I = computeInertiaTensor(flattened_state, true, stencil_size);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(I);
        Eigen::Vector3d evals = es.eigenvalues();

        std::cout << "Inertia tensor:\n" << I << "\n";
        std::cout << "Eigenvalues:\n" << evals.transpose() << "\n";

        int detectedClass = classifyCellViaInertia(evals);

        std::cout << "True Class: " << classNames[trueClass]
                  << ", Detected Class: " << classNames[detectedClass] << "\n";
    }

    std::cout << "--------------------------------------\n";

    // ---- Part 2: 10 random examples per class without visualize ----
    int numSamples = 100;
    for (int trueClass = 0; trueClass < 4; ++trueClass) {
        // Counters
        std::map<int, int> counts = { {0,0}, {1,0}, {2,0}, {3,0} };

        for (int s = 0; s < numSamples; ++s) {
            std::vector<double> flattened_state =
                data_gen.generate_State(trueClass, stencil_size, true, false);

            Eigen::Matrix3d I = computeInertiaTensor(flattened_state, true, stencil_size);

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(I);
            Eigen::Vector3d evals = es.eigenvalues();

            int detectedClass = classifyCellViaInertia(evals);
            counts[detectedClass]++;
        }

        std::cout << "True Class: " << classNames[trueClass] << ", Detected classes: "
                  << counts[0] << " paraboloids, "
                  << counts[1] << " ligaments, "
                  << counts[2] << " spheres, "
                  << counts[3] << " sheets.\n";
    }
                  
    //data_gen.generate_State(2, stencil_size, true, true);

    return 0;
}