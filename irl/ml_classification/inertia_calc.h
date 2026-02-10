#pragma once
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>

namespace IRL {

// Compute 2nd moment (central) from stencil data
// flattened_state interpretation depends on 'from_ith_moment':
//   from_ith_moment = 0 → use geometric cell centers
//   from_ith_moment = 1 → read first moments from flattened_state (current behavior)
//   from_ith_moment = 3 → RESERVED: transported moment immediately (to be implemented)
//
// Returns the central 2nd moment tensor about the neighborhood liquid centroid:
//   mu = ∑_cells vol * (r r^T), where r = (cell_centroid - c_liq)
inline Eigen::Matrix3d compute2ndMoment(const std::vector<double>& flattened_state,
                                       int stencil_size,
                                       int from_ith_moment = 1,
                                       double machineZero = 1e-12,
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
                    // transported first moment from input
                    double mx = flattened_state[idx++];
                    double my = flattened_state[idx++];
                    double mz = flattened_state[idx++];
                    c.firstMoment = Eigen::Vector3d(mx, my, mz);
                }
                else if (from_ith_moment == 0) {
                    // assume centroid = geometric cell center
                    double cx = (i + 0.5) - 0.5 * stencil_size;
                    double cy = (j + 0.5) - 0.5 * stencil_size;
                    double cz = (k + 0.5) - 0.5 * stencil_size;
                    double vol = c.alpha * cell_volume;
                    c.firstMoment = vol * Eigen::Vector3d(cx, cy, cz);
                }
                else if (from_ith_moment == 3) {
                    // WIP
                    // placeholder
                    throw std::invalid_argument("Invalid mode for compute2ndMoment");
                }
                else {
                    throw std::invalid_argument("Invalid mode for compute2ndMoment");
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

    // Assemble central 2nd moment tensor about c_liq
    Eigen::Matrix3d mu = Eigen::Matrix3d::Zero();
    for (auto& c : cells) {
        double vol = c.alpha * cell_volume;
        if (vol < machineZero) continue;

        Eigen::Vector3d centroid = c.firstMoment / vol;
        Eigen::Vector3d r = centroid - c_liq;

        mu += vol * (r * r.transpose());
    }

    return mu;
}

// Compute inertia tensor from stencil data
// flattened_state interpretation depends on 'from_ith_moment':
//   from_ith_moment = 0 → use geometric cell centers
//   from_ith_moment = 1 → read first moments from flattened_state (current behavior)
//   from_ith_moment = 3 → RESERVED: transported moment immediately (to be implemented)
inline Eigen::Matrix3d computeInertiaTensor(const std::vector<double>& flattened_state,
                                            int stencil_size,
                                            int from_ith_moment = 1,
                                            double machineZero = 1e-12,
                                            double cell_volume = 1.0)
{
    // First compute the central 2nd moment tensor mu about c_liq
    Eigen::Matrix3d mu = compute2ndMoment(flattened_state, stencil_size, from_ith_moment,
                                          machineZero, cell_volume);

    // Convert central 2nd moment -> inertia tensor:
    //   I = tr(mu) * Identity - mu
    Eigen::Matrix3d I = mu.trace() * Eigen::Matrix3d::Identity() - mu;

    return I;
}

} // namespace IRL
