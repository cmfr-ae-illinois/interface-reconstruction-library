// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_MOF2M_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_MOF2M_H_

#include "examples/amrex_advector/reconstruction_hybrid.h"
#include "examples/amrex_advector/reconstruction_pu.h"
#include "examples/amrex_advector/reconstruction_vf2.h"
#include "irl/amrex/sepunion_multifab.h"

#include <Eigen/QR>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

using namespace amrex;

struct MOF2M {
  static void GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr) {
    MOF1::GetReconstruction(interface, interface_with_ghost, moments, geom,
                            nullptr);

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const Real vol = dx[0] * dx[1] * dx[2];

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);
      Array4<Real const> moments_array = moments.const_array(mfi);
      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // Compute cell volume fraction
        const double vfrac = moments_array(i, j, k, comp_vf);
        // Skip cell if vfrac ~0 or vfrac~1
        if (vfrac < IRL::global_constants::VF_LOW ||
            vfrac > IRL::global_constants::VF_HIGH) {
          const double distance = std::copysign(1.0, vfrac - 0.5);
          interface_array(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          return;
        }

        // Construct local cell as IRL object
        const double x = problo[0] + i * dx[0];
        const double y = problo[1] + j * dx[1];
        const double z = problo[2] + k * dx[2];
        const IRL::Pt lower_cell_pt(x, y, z);
        const IRL::Pt upper_cell_pt(x + dx[0], y + dx[1], z + dx[2]);
        const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
        // Compute gas and liquid moments (M0 and M1)
        const double liq_m0 = moments_array(i, j, k, comp_m0);
        const double gas_m0 = vol - liq_m0;
        const IRL::Pt liq_m1 =
            (1.0 / liq_m0) * IRL::Pt(moments_array(i, j, k, comp_m1_l),
                                     moments_array(i, j, k, comp_m1_l + 1),
                                     moments_array(i, j, k, comp_m1_l + 2));
        const IRL::Pt gas_m1 =
            (1.0 / gas_m0) * IRL::Pt(moments_array(i, j, k, comp_m1_g),
                                     moments_array(i, j, k, comp_m1_g + 1),
                                     moments_array(i, j, k, comp_m1_g + 2));
        auto centroid_line = IRL::Normal(gas_m1 - liq_m1);
        IRL::GeneralMoments3D<2> liquid_moments, gas_moments;
        liquid_moments[0] = moments_array(i, j, k, comp_m0);
        liquid_moments[1] = moments_array(i, j, k, comp_m1_l);
        liquid_moments[2] = moments_array(i, j, k, comp_m1_l + 1);
        liquid_moments[3] = moments_array(i, j, k, comp_m1_l + 2);
        liquid_moments[4] = moments_array(i, j, k, comp_m2_l);
        liquid_moments[5] = moments_array(i, j, k, comp_m2_l + 1);
        liquid_moments[6] = moments_array(i, j, k, comp_m2_l + 2);
        liquid_moments[7] = moments_array(i, j, k, comp_m2_l + 3);
        liquid_moments[8] = moments_array(i, j, k, comp_m2_l + 4);
        liquid_moments[9] = moments_array(i, j, k, comp_m2_l + 5);
        gas_moments[0] = vol - moments_array(i, j, k, comp_m0);
        gas_moments[1] = moments_array(i, j, k, comp_m1_g);
        gas_moments[2] = moments_array(i, j, k, comp_m1_g + 1);
        gas_moments[3] = moments_array(i, j, k, comp_m1_g + 2);
        gas_moments[4] = moments_array(i, j, k, comp_m2_g);
        gas_moments[5] = moments_array(i, j, k, comp_m2_g + 1);
        gas_moments[6] = moments_array(i, j, k, comp_m2_g + 2);
        gas_moments[7] = moments_array(i, j, k, comp_m2_g + 3);
        gas_moments[8] = moments_array(i, j, k, comp_m2_g + 4);
        gas_moments[9] = moments_array(i, j, k, comp_m2_g + 5);

        MOF2Functor myMOF2Functor(5, 20, cell, liquid_moments, gas_moments,
                                  cell, vfrac);
        myMOF2Functor.setguess(centroid_line, interface_array(i, j, k));
        Eigen::NumericalDiff<MOF2Functor, Eigen::Forward> NDMOF2Functor(
            myMOF2Functor, 1.0e-12);
        Eigen::LevenbergMarquardt<
            Eigen::NumericalDiff<MOF2Functor, Eigen::Forward>, double>
            MOF2LM(NDMOF2Functor);
        MOF2LM.parameters.ftol =
            std::sqrt(std::numeric_limits<double>::epsilon());
        MOF2LM.parameters.xtol =
            std::sqrt(std::numeric_limits<double>::epsilon());
        MOF2LM.parameters.factor = 1.0e-3;
        MOF2LM.parameters.maxfev = 500;
        Eigen::VectorXd x_vec(5);
        for (int ii = 0; ii < 5; ii++) x_vec(ii) = 1.0;
        Eigen::VectorXd fvec(20);
        myMOF2Functor.errorvec(x_vec, fvec);
        const double init_error_li = fvec.lpNorm<Eigen::Infinity>();
        double final_error_li = init_error_li;
        int it = 0, itmax = 5;
        do {
          for (int ii = 0; ii < 5; ii++) x_vec(ii) = 1.0;
          Eigen::LevenbergMarquardtSpace::Status status =
              MOF2LM.minimizeInit(x_vec);
          do {
            status = MOF2LM.minimizeOneStep(x_vec);
          } while (status == Eigen::LevenbergMarquardtSpace::Running);
          myMOF2Functor.errorvec(x_vec, fvec);
          final_error_li = fvec.lpNorm<Eigen::Infinity>();
          MOF2LM.parameters.factor *= 10.0;
          it++;
        } while (init_error_li < final_error_li && it < itmax);

        IRL::Paraboloid paraboloid = myMOF2Functor.getparaboloid(x_vec);
        interface_array(i, j, k) = paraboloid;
      });
    }  // end mfi

    // This copies the interface, update ghosts, and shifts datums across
    // periodic boundaries
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

struct SuperMOF2M {
  static void GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr) {
    MOF1::GetReconstruction(interface, interface_with_ghost, moments, geom,
                            nullptr);

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const Real vol = dx[0] * dx[1] * dx[2];

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);
      Array4<Real const> moments_array = moments.const_array(mfi);
      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // Compute cell volume fraction
        const double vfrac = moments_array(i, j, k, comp_vf);
        // Skip cell if vfrac ~0 or vfrac~1
        if (vfrac < IRL::global_constants::VF_LOW ||
            vfrac > IRL::global_constants::VF_HIGH) {
          const double distance = std::copysign(1.0, vfrac - 0.5);
          interface_array(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          return;
        }

        // Construct local cell as IRL object
        const double x = problo[0] + i * dx[0];
        const double y = problo[1] + j * dx[1];
        const double z = problo[2] + k * dx[2];
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(x, y, z), IRL::Pt(x + dx[0], y + dx[1], z + dx[2]));
        const auto super_cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(x - dx[0], y - dx[1], z - dx[2]),
            IRL::Pt(x + 2.0 * dx[0], y + 2.0 * dx[1], z + 2.0 * dx[2]));
        auto super_l_mom = IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);
        auto super_g_mom = IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);
        for (int ii = -1; ii <= 1; ++ii) {
          for (int jj = -1; jj <= 1; ++jj) {
            for (int kk = -1; kk <= 1; ++kk) {
              super_l_mom[0] += moments_array(i + ii, j + jj, k + kk, comp_m0);
              super_l_mom[1] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m1_l);
              super_l_mom[2] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m1_l + 1);
              super_l_mom[3] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m1_l + 2);
              super_l_mom[4] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_l);
              super_l_mom[5] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_l + 1);
              super_l_mom[6] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_l + 2);
              super_l_mom[7] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_l + 3);
              super_l_mom[8] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_l + 4);
              super_l_mom[9] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_l + 5);
              super_g_mom[0] +=
                  vol - moments_array(i + ii, j + jj, k + kk, comp_m0);
              super_g_mom[1] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m1_g);
              super_g_mom[2] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m1_g + 1);
              super_g_mom[3] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m1_g + 2);
              super_g_mom[4] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_g);
              super_g_mom[5] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_g + 1);
              super_g_mom[6] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_g + 2);
              super_g_mom[7] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_g + 3);
              super_g_mom[8] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_g + 4);
              super_g_mom[9] +=
                  moments_array(i + ii, j + jj, k + kk, comp_m2_g + 5);
            }
          }
        }

        IRL::Pt liquid_centroid =
            IRL::Pt(super_l_mom[1], super_l_mom[2], super_l_mom[3]);
        IRL::Pt gas_centroid =
            IRL::Pt(super_g_mom[1], super_g_mom[2], super_g_mom[3]);
        liquid_centroid *= 1.0 / IRL::safelyEpsilon(super_l_mom[0]);
        gas_centroid *= 1.0 / IRL::safelyEpsilon(super_g_mom[0]);
        auto centroid_line = IRL::Normal(gas_centroid - liquid_centroid);
        MOF2Functor myMOF2Functor(5, 20, super_cell, super_l_mom, super_g_mom,
                                  cell, vfrac);
        myMOF2Functor.setguess(centroid_line, interface_array(i, j, k));
        Eigen::NumericalDiff<MOF2Functor, Eigen::Forward> NDMOF2Functor(
            myMOF2Functor, 1.0e-12);
        Eigen::LevenbergMarquardt<
            Eigen::NumericalDiff<MOF2Functor, Eigen::Forward>, double>
            MOF2LM(NDMOF2Functor);
        MOF2LM.parameters.ftol =
            std::sqrt(std::numeric_limits<double>::epsilon());
        MOF2LM.parameters.xtol =
            std::sqrt(std::numeric_limits<double>::epsilon());
        MOF2LM.parameters.factor = 1.0e-3;
        MOF2LM.parameters.maxfev = 500;
        Eigen::VectorXd x_vec(5);
        for (int ii = 0; ii < 5; ii++) x_vec(ii) = 1.0;
        // MOF2LM.minimize(x);
        Eigen::VectorXd fvec(20);
        myMOF2Functor.errorvec(x_vec, fvec);
        const double init_error_li = fvec.lpNorm<Eigen::Infinity>();
        double final_error_li = init_error_li;
        int it = 0, itmax = 5;
        do {
          for (int ii = 0; ii < 5; ii++) x_vec(ii) = 1.0;
          Eigen::LevenbergMarquardtSpace::Status status =
              MOF2LM.minimizeInit(x_vec);
          do {
            status = MOF2LM.minimizeOneStep(x_vec);
          } while (status == Eigen::LevenbergMarquardtSpace::Running);
          myMOF2Functor.errorvec(x_vec, fvec);
          final_error_li = fvec.lpNorm<Eigen::Infinity>();
          MOF2LM.parameters.factor *= 10.0;
          it++;
        } while (init_error_li < final_error_li && it < itmax);

        IRL::Paraboloid paraboloid = myMOF2Functor.getparaboloid(x_vec);
        interface_array(i, j, k) = paraboloid;
      });
    }  // end mfi

    // This copies the interface, update ghosts, and shifts datums across
    // periodic boundaries
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_MOF2M_H_
