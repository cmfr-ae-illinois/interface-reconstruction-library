#pragma once
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include "common_functions.h"

namespace IRL {

// Compute 2nd moment (central) from stencil data
// flattened_state interpretation depends on 'from_ith_moment':
//   from_ith_moment = 0 → use geometric cell centers
//   from_ith_moment = 1 → read first moments from flattened_state (current behavior)
//   from_ith_moment = 3 → RESERVED: transported moment immediately (to be implemented)
//
// Returns the central 2nd moment tensor about the neighborhood liquid centroid:
//   mu = ∑_cells vol * (r r^T), where r = (cell_centroid - c_liq)
inline Eigen::Matrix3d compute2ndMoment(const std::vector<float>& flat_stencil,
                                        int stencil_size,
                                        int from_ith_moment = 1,
                                        bool include_Surface_Area = false,
                                        double machineZero = 1e-12,
                                        double cell_volume = 1.0)
{
    // This returns the CENTRAL 2nd moment tensor
    //   mu = sum_cells V * (r r^T),
    // where r = (cell centroid - liquid centroid of stencil).

    std::vector<CellData> stencil;
    unpackStencil(flat_stencil,
                  stencil,
                  nullptr,   // do not read stored global 2nd moments here
                  nullptr,   // do not read eigenvalues here
                  stencil_size,
                  from_ith_moment,
                  include_Surface_Area,
                  false);

    Eigen::Vector3d c_liq = Eigen::Vector3d::Zero();
    double Vsum = 0.0;

    // 1) Compute liquid centroid of the stencil
    if (from_ith_moment >= 1) {
        // mx,my,mz are stored as first moments = V * centroid
        for (const auto& c : stencil) {
            const double vol = static_cast<double>(c.vfrac) * cell_volume;
            if (vol <= machineZero) continue;

            Vsum += vol;
            c_liq += Eigen::Vector3d(
                static_cast<double>(c.mx),
                static_cast<double>(c.my),
                static_cast<double>(c.mz)
            );
        }

        if (Vsum > machineZero) {
            c_liq /= Vsum;
        } else {
            return Eigen::Matrix3d::Zero();
        }
    } else {
        float cx = 0.0f, cy = 0.0f, cz = 0.0f, Vsum_f = 0.0f;
        approximateCentroidFromVfrac(stencil, stencil_size, cx, cy, cz, Vsum_f);

        c_liq = Eigen::Vector3d(
            static_cast<double>(cx),
            static_cast<double>(cy),
            static_cast<double>(cz)
        );
        Vsum = static_cast<double>(Vsum_f) * cell_volume;

        if (Vsum <= machineZero) {
            return Eigen::Matrix3d::Zero();
        }
    }

    // 2) Assemble central 2nd moment tensor
    Eigen::Matrix3d mu = Eigen::Matrix3d::Zero();
    const double c0 = 0.5 * (static_cast<double>(stencil_size) - 1.0);

    for (int i = 0; i < stencil_size; ++i) {
        for (int j = 0; j < stencil_size; ++j) {
            for (int k = 0; k < stencil_size; ++k) {
                const CellData& c = stencil[cellIndex(i, j, k, stencil_size)];
                const double vol = static_cast<double>(c.vfrac) * cell_volume;

                if (vol <= machineZero) continue;

                // assume mx,my,mz are first moments
                Eigen::Vector3d x(
                    static_cast<double>(c.mx) / vol,
                    static_cast<double>(c.my) / vol,
                    static_cast<double>(c.mz) / vol
                );

                Eigen::Vector3d r = x - c_liq;
                mu += vol * (r * r.transpose());
            }
        }
    }

    return mu;
}

// Compute inertia tensor from stencil data
// flattened_state interpretation depends on 'from_ith_moment':
//   from_ith_moment = 0 → use geometric cell centers
//   from_ith_moment = 1 → read first moments from flattened_state (current behavior)
//   from_ith_moment = 3 → RESERVED: transported moment immediately (to be implemented)
inline Eigen::Matrix3d computeInertiaTensor(const std::vector<float>& flattened_state,
                                            int stencil_size,
                                            int from_ith_moment = 1,
                                            bool include_Surface_Area = false,
                                            double machineZero = 1e-12,
                                            double cell_volume = 1.0)
{
    // First compute the central 2nd moment tensor mu about c_liq
    Eigen::Matrix3d mu = compute2ndMoment(flattened_state, stencil_size, from_ith_moment, include_Surface_Area,
                                          machineZero, cell_volume);

    // Convert central 2nd moment -> inertia tensor:
    //   I = tr(mu) * Identity - mu
    Eigen::Matrix3d I = mu.trace() * Eigen::Matrix3d::Identity() - mu;

    return I;
}

inline void appendInertiaEigenvalues(std::vector<float>& flattened_state,
                             int stencil_size,
                             int include_Moments = 1,
                             int from_ith_moment = 1,
                             bool include_Surface_Area = false,
                             double machineZero = 1e-12)
{
    Eigen::Matrix3d I = IRL::computeInertiaTensor(flattened_state, stencil_size, from_ith_moment, include_Surface_Area, machineZero);

    // Get eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(I);
    Eigen::Vector3d evals = solver.eigenvalues();

    // Sort eigenvalues descending: I1 >= I2 >= I3
    std::sort(evals.data(), evals.data() + 3, std::greater<double>());
    double I1 = evals[0], I2 = evals[1], I3 = evals[2];

    if (include_Moments <= 3) {
        flattened_state.push_back(static_cast<float>(I1));
        flattened_state.push_back(static_cast<float>(I2));
        flattened_state.push_back(static_cast<float>(I3));
    }

    if (include_Moments == 4) {
        // if include_Moments == 4, use only the three eigenvalues
        flattened_state.clear();
        flattened_state.push_back(static_cast<float>(I1));
        flattened_state.push_back(static_cast<float>(I2));
        flattened_state.push_back(static_cast<float>(I3));
    }

}

} // namespace IRL
