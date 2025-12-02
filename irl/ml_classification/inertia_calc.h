#pragma once
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>

namespace IRL {

// Compute inertia tensor from stencil data
// flattened_state interpretation depends on 'from_ith_moment':
//   from_ith_moment = 0 → use geometric cell centers
//   from_ith_moment = 1 → read first moments from flattened_state (current behavior)
//   from_ith_moment = 3 → RESERVED: transported moment immediately (to be implemented)
Eigen::Matrix3d computeInertiaTensor(const std::vector<double>& flattened_state,
                                     int stencil_size,
                                     int from_ith_moment = 1,
                                     double machineZero = 1e-14,
                                     double cell_volume = 1.0)
{
    struct Cell {
        double alpha;
        Eigen::Vector3d firstMoment; // stores V * centroid
    };

    std::vector<Cell> cells;
    cells.reserve(stencil_size * stencil_size * stencil_size);

    size_t idx = 0;
    for (int i = 0; i < stencil_size; ++i) {
        for (int j = 0; j < stencil_size; ++j) {
            for (int k = 0; k < stencil_size; ++k) {
                Cell c;
                c.alpha = flattened_state[idx++];

                if (from_ith_moment == 1) {
                    // --- transported first moment from input ---
                    double mx = flattened_state[idx++];
                    double my = flattened_state[idx++];
                    double mz = flattened_state[idx++];
                    c.firstMoment = Eigen::Vector3d(mx, my, mz);
                }
                else if (from_ith_moment == 0) {
                    // --- assume centroid = geometric cell center ---
                    double cx = (i + 0.5);
                    double cy = (j + 0.5);
                    double cz = (k + 0.5);
                    double vol = c.alpha * cell_volume;
                    c.firstMoment = vol * Eigen::Vector3d(cx, cy, cz);
                }
                else if (from_ith_moment == 3) {
                    // WIP
                    // placeholder
                    throw std::invalid_argument("Invalid mode for computeInertiaTensor");
                }
                else {
                    throw std::invalid_argument("Invalid mode for computeInertiaTensor");
                }

                cells.push_back(c);
            }
        }
    }

    // Compute neighborhood liquid centroid
    double m_total = 0.0;
    Eigen::Vector3d c_liq = Eigen::Vector3d::Zero();
    for (auto& c : cells) {
        double vol = c.alpha * cell_volume;
        m_total += vol;
        c_liq += c.firstMoment;
    }
    if (m_total > machineZero) {
        c_liq /= m_total;
    }

    // Assemble inertia tensor about c_liq
    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    for (auto& c : cells) {
        double vol = c.alpha * cell_volume;
        if (vol < machineZero) continue;

        Eigen::Vector3d centroid = c.firstMoment / vol;
        Eigen::Vector3d r = centroid - c_liq;
        double x = r[0], y = r[1], z = r[2];

        I(0,0) += vol * (y*y + z*z);
        I(1,1) += vol * (x*x + z*z);
        I(2,2) += vol * (x*x + y*y);

        I(0,1) -= vol * x * y; I(1,0) -= vol * x * y;
        I(0,2) -= vol * x * z; I(2,0) -= vol * x * z;
        I(1,2) -= vol * y * z; I(2,1) -= vol * y * z;
    }

    return I;
}

} // namespace IRL