// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_MOF2_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_MOF2_H_

#include "examples/amrex_advector/reconstruction_hybrid.h"
#include "examples/amrex_advector/reconstruction_pu.h"
#include "examples/amrex_advector/reconstruction_vf2.h"
#include "irl/amrex/sepunion_multifab.h"

#include <Eigen/QR>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

using namespace amrex;

void RecenterMoments(IRL::GeneralMoments3D<2>* moments, const IRL::Pt& center) {
  const double M0 = (*moments)[0];
  const Eigen::Matrix<double, 3, 1> M1{(*moments)[1], (*moments)[2],
                                       (*moments)[3]};
  const Eigen::Matrix<double, 3, 3> M2{
      {(*moments)[4], (*moments)[5], (*moments)[6]},
      {(*moments)[5], (*moments)[7], (*moments)[8]},
      {(*moments)[6], (*moments)[8], (*moments)[9]}};
  const Eigen::Matrix<double, 3, 1> C{center[0], center[1], center[2]};
  const auto M1_final = M1 - M0 * C;
  const auto M2_final =
      M2 - M1 * C.transpose() - C * M1.transpose() + M0 * C * C.transpose();
  (*moments)[1] = M1_final(0);
  (*moments)[2] = M1_final(1);
  (*moments)[3] = M1_final(2);
  (*moments)[4] = M2_final(0, 0);
  (*moments)[5] = M2_final(0, 1);
  (*moments)[6] = M2_final(0, 2);
  (*moments)[7] = M2_final(1, 1);
  (*moments)[8] = M2_final(1, 2);
  (*moments)[9] = M2_final(2, 2);
}

struct MOF2Functor {
  typedef Eigen::VectorXd InputType;
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  typedef double Scalar;
  enum {
    InputsAtCompileTime = Eigen::Dynamic,
    ValuesAtCompileTime = Eigen::Dynamic
  };

  // Variables
  int m_inputs, m_values;
  double m_cell_volume, m_vfrac, m_a, m_b;
  IRL::RectangularCuboid m_cell, m_cell_constraint;
  IRL::GeneralMoments3D<2> m_liquid_moments, m_gas_moments;
  IRL::ReferenceFrame m_frame;
  IRL::Pt m_datum, m_liquid_centroid, m_gas_centroid, m_cell_centroid;
  double m_length_scale, m_m0_scale, m_m1_scale_liquid, m_m1_scale_gas,
      m_m2_scale_liquid, m_m2_scale_gas, m_vfrac_constraint;

  // Constructor
  MOF2Functor(int inputs, int values, const IRL::RectangularCuboid& cell,
              const IRL::GeneralMoments3D<2>& liquid_moments,
              const IRL::GeneralMoments3D<2>& gas_moments,
              const IRL::RectangularCuboid& constraint_cell,
              const double constraint_vfrac)
      : m_inputs(inputs),
        m_values(values),
        m_cell(cell),
        m_cell_volume(IRL::safelyTiny(cell.calculateVolume())),
        m_cell_centroid(cell.calculateCentroid()),
        m_liquid_moments(liquid_moments),
        m_gas_moments(gas_moments),
        m_vfrac_constraint(constraint_vfrac),
        m_cell_constraint(constraint_cell),
        m_liquid_centroid(
            IRL::Pt(liquid_moments[1], liquid_moments[2], liquid_moments[3]) *
            (1.0 / IRL::safelyTiny(liquid_moments[0]))),
        m_gas_centroid(IRL::Pt(gas_moments[1], gas_moments[2], gas_moments[3]) *
                       (1.0 / IRL::safelyTiny(gas_moments[0]))) {
    m_vfrac = liquid_moments[0] / m_cell_volume;
    m_length_scale = std::cbrt(m_cell_volume);
    m_m0_scale = std::pow(IRL::safelyTiny(m_liquid_moments[0]), -1.0);
    m_m1_scale_liquid =
        std::pow(IRL::safelyTiny(m_liquid_moments[0]), -4.0 / 3.0);
    m_m1_scale_gas = std::pow(IRL::safelyTiny(m_gas_moments[0]), -4.0 / 3.0);
    m_m2_scale_liquid =
        std::pow(IRL::safelyTiny(m_liquid_moments[0]), -5.0 / 3.0);
    m_m2_scale_gas = std::pow(IRL::safelyTiny(m_gas_moments[0]), -5.0 / 3.0);
    m_datum =
        IRL::Pt((1.0 - m_vfrac) * m_liquid_centroid + m_vfrac * m_gas_centroid);
    RecenterMoments(&m_liquid_moments, m_liquid_centroid);
    RecenterMoments(&m_gas_moments, m_gas_centroid);
  }

  void setguess(const IRL::Normal& guess_direction,
                const IRL::SeparatorUnion& interface) {
    if (interface.type() == IRL::SeparatorUnion::SeparatorType::OnePlane) {
      m_frame = IRL::ReferenceFrame::fromNormal(interface.getPlane().normal());
      m_a = 0.0;
      m_b = 0.0;
    } else if (interface.type() ==
               IRL::SeparatorUnion::SeparatorType::Paraboloid) {
      m_frame = interface.getParaboloid().getReferenceFrame();
      m_a = interface.getParaboloid().getAlignedParaboloid().a();
      m_b = interface.getParaboloid().getAlignedParaboloid().b();
    } else {
      m_frame = IRL::ReferenceFrame::fromNormal(guess_direction);
      m_a = 0.0;
      m_b = 0.0;
    }
  }

  const IRL::Paraboloid getparaboloid(const Eigen::VectorXd& x) const {
    const IRL::Pt datum = m_datum;
    IRL::ReferenceFrame newframe;
    IRL::Paraboloid paraboloid;
    // constexpr double vfrac_tol = 0.01;
    constexpr double vfrac_tol = 0.0;
    // For almost empty of full cells (<1% andd > 99%), we only allow rotation
    // around x,y only and same a and b coefficients
    if (m_vfrac < vfrac_tol || m_vfrac > 1.0 - vfrac_tol) {
      const auto rotation = IRL::UnitQuaternion(2.0 * M_PI * x(1), m_frame[1]) *
                            IRL::UnitQuaternion(2.0 * M_PI * x(0), m_frame[0]);
      newframe = rotation * m_frame;
      double delta = (x(2) - 0.999) / m_length_scale;
      delta = (std::abs(delta * m_length_scale) > 3.0)
                  ? std::copysign(3.0 / m_length_scale, delta)
                  : delta;
      paraboloid = IRL::Paraboloid(datum, newframe, delta, delta);
    }
    // For any other cells, we allow rotations around x,y,z and free
    // coefficients
    else {
      const auto rotation = IRL::UnitQuaternion(2.0 * M_PI * x(4), m_frame[2]) *
                            IRL::UnitQuaternion(2.0 * M_PI * x(1), m_frame[1]) *
                            IRL::UnitQuaternion(2.0 * M_PI * x(0), m_frame[0]);
      newframe = rotation * m_frame;
      const double delta_a = (x(2) - 0.999) / m_length_scale;
      const double delta_b = (x(3) - 1.001) / m_length_scale;
      double coeff_a = m_a + delta_a;
      double coeff_b = m_b + delta_b;
      coeff_a = (std::abs(coeff_a * m_length_scale) > 3.0)
                    ? std::copysign(3.0 / m_length_scale, coeff_a)
                    : coeff_a;
      coeff_b = (std::abs(coeff_b * m_length_scale) > 3.0)
                    ? std::copysign(3.0 / m_length_scale, coeff_b)
                    : coeff_b;
      paraboloid = IRL::Paraboloid(datum, newframe, coeff_a, coeff_b);
    }
    IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
        solver_distance(m_cell_constraint, m_vfrac_constraint, 1.0e-14,
                        paraboloid);
    auto new_datum =
        IRL::Pt(datum + solver_distance.getDistance() * newframe[2]);
    paraboloid.setDatum(new_datum);
    return paraboloid;
  }

  void errorvec(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    const auto paraboloid = this->getparaboloid(x);
    const auto sep_moments =
        IRL::getVolumeMoments<IRL::SeparatedMoments<IRL::GeneralMoments3D<2>>>(
            m_cell, paraboloid);
    auto liquid_moments = sep_moments[0];
    auto gas_moments = sep_moments[1];
    RecenterMoments(&liquid_moments, m_liquid_centroid);
    RecenterMoments(&gas_moments, m_gas_centroid);
    fvec.setZero();
    const double k1dx =
        2.0 * (m_a + (x(2) - 0.999) / m_length_scale) * m_length_scale;
    const double k2dx =
        2.0 * (m_b + (x(3) - 1.001) / m_length_scale) * m_length_scale;
    constexpr double maxkdx = 6.0;
    constexpr double mu = 50.0;
    // Penalty to prevent kappa * dx > maxkdx
    fvec(0) = mu * ((std::max(0.0, std::abs(k1dx) - maxkdx)) +
                    (std::max(0.0, std::abs(k2dx) - maxkdx)) +
                    (std::max(0.0, std::abs(1.0 - x(0)) - 0.05)) +
                    (std::max(0.0, std::abs(1.0 - x(1)) - 0.05)));
    // Next line is commented out because M0 constraint is enforced by
    // getparaboloid() call
    // fvec(1) = m_m0_scale * (liquid_moments[0] - m_liquid_moments[0]);
    for (int i = 0; i < 3; ++i) {
      fvec(2 + i) =
          m_m1_scale_liquid * (liquid_moments[1 + i] - m_liquid_moments[1 + i]);
      fvec(5 + i) =
          m_m1_scale_gas * (gas_moments[1 + i] - m_gas_moments[1 + i]);
    }
    for (int i = 0; i < 6; ++i) {
      fvec(8 + i) =
          m_m2_scale_liquid * (liquid_moments[4 + i] - m_liquid_moments[4 + i]);
      fvec(14 + i) =
          m_m2_scale_gas * (gas_moments[4 + i] - m_gas_moments[4 + i]);
    }
    fvec(9) *= 2.0;
    fvec(10) *= 2.0;
    fvec(12) *= 2.0;
    fvec(15) *= 2.0;
    fvec(16) *= 2.0;
    fvec(18) *= 2.0;
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    this->errorvec(x, fvec);
    return 0;
  }

  // Get methods for users to determine function input and output dimensions
  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

struct MOF2 {
  static int GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr,
      Real* reconstruction_loop_time = nullptr) {
    MOF1::GetReconstruction(interface, interface_with_ghost, moments, geom,
                            nullptr, reconstruction_loop_time);

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const Real vol = dx[0] * dx[1] * dx[2];

    ReconstructionLoopTimer loop_timer(reconstruction_loop_time);
    loop_timer.start();
    amrex::Gpu::DeviceScalar<int> counter(0);
    int* niter = counter.dataPtr();
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
            amrex::Gpu::Atomic::Add(niter, 1);
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
    loop_timer.stop();

    // This copies the interface, update ghosts, and shifts datums across
    // periodic boundaries
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
    return counter.dataValue();
  }
};

struct SuperMOF2 {
  static int GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr,
      Real* reconstruction_loop_time = nullptr) {
    MOF1::GetReconstruction(interface, interface_with_ghost, moments, geom,
                            nullptr, reconstruction_loop_time);

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const Real vol = dx[0] * dx[1] * dx[2];

    ReconstructionLoopTimer loop_timer(reconstruction_loop_time);
    loop_timer.start();
    amrex::Gpu::DeviceScalar<int> counter(0);
    int* niter = counter.dataPtr();
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
            amrex::Gpu::Atomic::Add(niter, 1);
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
    loop_timer.stop();

    // This copies the interface, update ghosts, and shifts datums across
    // periodic boundaries
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
    return counter.dataValue();
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_MOF2_H_
