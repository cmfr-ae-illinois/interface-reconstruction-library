// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry
// operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>
#include <random>

#include "irl/moments/general_moments.h"
#include "irl/moments/separated_volume_moments.h"

#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/rotations.h"
#include "irl/interface_reconstruction_methods/constrained_optimization_behavior.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"
#include "irl/optimization/constrained_levenberg_marquardt.h"
#include "irl/paraboloid_reconstruction/gradient_paraboloid.h"
#include "irl/paraboloid_reconstruction/hessian_paraboloid.h"

#include "gtest/gtest.h"

#include "irl/data_structures/small_vector.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron_paraboloid.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron_paraboloid.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/paraboloid_reconstruction/parametrized_surface.h"
#include "irl/planar_reconstruction/planar_separator.h"

#include <Eigen/Dense>  // Eigen header
#include "tests/src/basic_mesh.h"
#include "tests/src/data.h"
#include "tests/src/vtk.h"

namespace IRL {

// / Requirements for OptimizingClass:
// / - `double calculateScalarError(void)` : A method to calculate a scalar
// error / that we are trying to minimize. / -
// `Eigen::Matrix<double,kRows,1>calculateVectorError(void)` : A method that /
// returns the vector return (correct_values - guess_values) by value / - `void
// updateGuess(Eigen::Matrix<double,kColumns,1>)` : A method that takes / in the
// delta change and computes a new guess vector (which it is storing / itself)
// / - `void updateBestGuess(void)` : A method that updates the best guess
// / and all other things necessary before a new Jacobian
// / is calculated and a new step is taken.
// / - `Eigen::Matrix<double,kRows,1>calculateChangeInGuess(void)` : A method
// / that calculates the difference between Guess variables and bestGuess (for
// / use in calculating derivative in Jacobian).
// / - void increaseLambda(double*) : A method to increase the value of
// / lambda (for failed attempts at finding a new minimum).
// / - void decreaseLambda(double*) : A method to decrease the value of
// / lambda (for successful attempts at finding a new minimum).
// / - `bool errorTooHigh(const double)` : A method that takes a scalar error
// / and returns a boolean whether the error is low
// / enough to stop optimization and return.
// / -`bool iterationTooHigh(const int)` : A method that takes
// / the number of iterations and returns a bool
// / whether the maximum number of allowable iterations
// / has been exceeded.
// / - `minimumReached(const Eigen::Matrix<double,kColumns,1>)` : A method that
// / takes delta and determines if the optimization has reached a minimum,
// return / a bool `true` if optimization should exit.
// `shouldComputeJacobian(const int, / const int)` : A method that returns a
// bool for whether or not a jacobian / should be computed when given the
// current iteration and the last iteration / the Jacobian was computed for.
// /
// / kRows is the number of rows involved in the error vector of the
// /  Levenberg-Marquardt system [y-f].
// /
// / kColumns is the number of columns involved in the Jacobian for
// / the Levenberg-Marquardt system, equal to the number of parameters
// / being fit.

template <UnsignedIndex_t kRows = 20, UnsignedIndex_t kColumns = 8>
class Paraboloid_MOF_Object {
 public:
  Paraboloid_MOF_Object() = default;

  Paraboloid_MOF_Object(const RectangularCuboid a_cell, const Pt a_lower_corner,
                        const Pt a_upper_corner,
                        SeparatedMoments<GeneralMoments3D<2>> a_moments,
                        const Paraboloid a_paraboloid) {
    cell_m = a_cell;
    best_reconstruction_m = a_paraboloid;
    GeneralMoments3D<2>& moments_a = a_moments[0];
    GeneralMoments3D<2>& moments_b = a_moments[1];
    m0_ref_a_m = 1.0 * moments_a[0];
    // std::min(moments_a[0], moments_b[0]);
    //  moments_a[0];
    m0_ref_b_m = 1.0 * moments_b[0];
    // std::min(moments_a[0], moments_b[0]);
    //  moments_b[0];
    m1_ref_a_m = pow(m0_ref_a_m, 4.0 / 3.0);
    // m0_ref_a_m;
    // pow(m0_ref_a_m, 4.0 / 3.0);
    m1_ref_b_m = pow(m0_ref_b_m, 4.0 / 3.0);
    // m0_ref_b_m;
    // pow(m0_ref_b_m, 4.0 / 3.0);
    m2_ref_a_m = pow(m0_ref_a_m, 5.0 / 3.0);
    // m0_ref_a_m;
    // pow(m0_ref_a_m, 5.0 / 3.0);
    m2_ref_b_m = pow(m0_ref_b_m, 5.0 / 3.0);
    // m0_ref_b_m;
    // pow(m0_ref_b_m, 5.0 / 3.0);
    ref_length_m = pow(cell_m.calculateVolume(), 1.0 / 3.0);
    lower_corner_m = a_lower_corner;
    upper_corner_m = a_upper_corner;
    // a_moments.normalizeAsInvariant();
    correct_values_m.Zero();
    guess_values_m.Zero();
    // for (UnsignedIndex_t i = 0; i < kRows; i++) {
    //   correct_values_m(i) = a_moments[i];
    // }
    correct_values_m(0) = moments_a[0] / m0_ref_a_m;
    correct_values_m(1) = moments_a[1] / m1_ref_a_m;
    correct_values_m(2) = moments_a[2] / m1_ref_a_m;
    correct_values_m(3) = moments_a[3] / m1_ref_a_m;
    correct_values_m(4) = moments_a[4] / m2_ref_a_m;
    correct_values_m(5) = moments_a[5] / m2_ref_a_m;
    correct_values_m(6) = moments_a[6] / m2_ref_a_m;
    correct_values_m(7) = moments_a[7] / m2_ref_a_m;
    correct_values_m(8) = moments_a[8] / m2_ref_a_m;
    correct_values_m(9) = moments_a[9] / m2_ref_a_m;
    correct_values_m(10) = moments_b[0] / m0_ref_b_m;
    correct_values_m(11) = moments_b[1] / m1_ref_b_m;
    correct_values_m(12) = moments_b[2] / m1_ref_b_m;
    correct_values_m(13) = moments_b[3] / m1_ref_b_m;
    correct_values_m(14) = moments_b[4] / m2_ref_b_m;
    correct_values_m(15) = moments_b[5] / m2_ref_b_m;
    correct_values_m(16) = moments_b[6] / m2_ref_b_m;
    correct_values_m(17) = moments_b[7] / m2_ref_b_m;
    correct_values_m(18) = moments_b[8] / m2_ref_b_m;
    correct_values_m(19) = moments_b[9] / m2_ref_b_m;
  }

  double calculateScalarError(void) {
    return (correct_values_m - guess_values_m).squaredNorm();
  }

  Eigen::Matrix<double, kRows, 1> calculateVectorError(void) {
    return correct_values_m - guess_values_m;
  }

  void updateGuess(Eigen::Matrix<double, kColumns, 1>* a_delta) {
    // Get old datum, frame, and curvatures
    const auto& datum = best_reconstruction_m.getDatum();
    const auto& frame = best_reconstruction_m.getReferenceFrame();
    const auto& algnd_para = best_reconstruction_m.getAlignedParaboloid();
    // Create quaternion for rotating the frame
    UnitQuaternion x_rotation((*a_delta)(5), frame[0]);
    UnitQuaternion y_rotation((*a_delta)(6), frame[1]);
    UnitQuaternion z_rotation((*a_delta)(7), frame[2]);
    // auto total_rotation = x_rotation * y_rotation * z_rotation;
    auto total_rotation = z_rotation;
    total_rotation.normalize();
    // Shift datum
    // const Pt new_datum = datum + (*a_delta)(2) * frame[0] +
    //                      (*a_delta)(3) * frame[1] + (*a_delta)(4) * frame[2];
    const Pt new_datum = datum + (*a_delta)(4) * frame[2];
    // const Pt new_datum = datum;
    // Rotate frame
    // const ReferenceFrame new_frame = total_rotation * frame;
    const ReferenceFrame new_frame = frame;
    // Update curvatures
    const double new_a = algnd_para.a() + (*a_delta)(0);
    const double new_b = algnd_para.b() + (*a_delta)(1);
    // Return new paraboloid guess
    guess_reconstruction_m = Paraboloid(new_datum, new_frame, new_a, new_b);
    guess_values_m.setZero();
    // auto moments =
    //     getVolumeMoments<VolumeMoments>(cell_m, guess_reconstruction_m);
    // guess_values_m(0) = moments.volume() / m0_ref_m;
    // guess_values_m(1) = moments.centroid()[0] / m1_ref_m;
    // guess_values_m(2) = moments.centroid()[1] / m1_ref_m;
    // guess_values_m(3) = moments.centroid()[2] / m1_ref_m;

    auto moments = getVolumeMoments<SeparatedMoments<GeneralMoments3D<2>>>(
        cell_m, guess_reconstruction_m);
    GeneralMoments3D<2>& moments_a = moments[0];
    GeneralMoments3D<2>& moments_b = moments[1];
    // moments.normalizeAsInvariant();
    // for (UnsignedIndex_t i = 0; i < kRows; ++i) {
    //   guess_values_m(i) = moments[i];
    // }
    guess_values_m(0) = moments_a[0] / m0_ref_a_m;
    guess_values_m(1) = moments_a[1] / m1_ref_a_m;
    guess_values_m(2) = moments_a[2] / m1_ref_a_m;
    guess_values_m(3) = moments_a[3] / m1_ref_a_m;
    guess_values_m(4) = moments_a[4] / m2_ref_a_m;
    guess_values_m(5) = moments_a[5] / m2_ref_a_m;
    guess_values_m(6) = moments_a[6] / m2_ref_a_m;
    guess_values_m(7) = moments_a[7] / m2_ref_a_m;
    guess_values_m(8) = moments_a[8] / m2_ref_a_m;
    guess_values_m(9) = moments_a[9] / m2_ref_a_m;
    guess_values_m(10) = moments_b[0] / m0_ref_b_m;
    guess_values_m(11) = moments_b[1] / m1_ref_b_m;
    guess_values_m(12) = moments_b[2] / m1_ref_b_m;
    guess_values_m(13) = moments_b[3] / m1_ref_b_m;
    guess_values_m(14) = moments_b[4] / m2_ref_b_m;
    guess_values_m(15) = moments_b[5] / m2_ref_b_m;
    guess_values_m(16) = moments_b[6] / m2_ref_b_m;
    guess_values_m(17) = moments_b[7] / m2_ref_b_m;
    guess_values_m(18) = moments_b[8] / m2_ref_b_m;
    guess_values_m(19) = moments_b[9] / m2_ref_b_m;
  }

  void calculateJacobian(
      Eigen::Matrix<double, kColumns, kRows>* a_jacobian_transpose,
      Eigen::Matrix<double, kColumns, kColumns>* a_jacTjac) {
    // // Set up temporary delta
    // Eigen::Matrix<double, kColumns, 1> solo_delta;
    // // Calculate tranpose of Jacobian
    // using MyGradientType = ParaboloidGradientLocalBase<double>;
    // using MyScalarType = ScalarWithGradientBase<double, MyGradientType>;
    // using MyPtType = PtBase<MyScalarType>;
    // using MyMomentType = VolumeMomentsBase<MyScalarType>;
    // auto dp_frame = best_reconstruction_m.getReferenceFrame();
    // auto dp_datum = best_reconstruction_m.getDatum();
    // auto dp_aligned = best_reconstruction_m.getAlignedParaboloid();
    // auto frame = ReferenceFrameBase<MyScalarType>(
    //     NormalBase<MyScalarType>(MyScalarType(dp_frame[0][0]),
    //                              MyScalarType(dp_frame[0][1]),
    //                              MyScalarType(dp_frame[0][2])),
    //     NormalBase<MyScalarType>(MyScalarType(dp_frame[1][0]),
    //                              MyScalarType(dp_frame[1][1]),
    //                              MyScalarType(dp_frame[1][2])),
    //     NormalBase<MyScalarType>(MyScalarType(dp_frame[2][0]),
    //                              MyScalarType(dp_frame[2][1]),
    //                              MyScalarType(dp_frame[2][2])));
    // auto datum = PtBase<MyScalarType>(MyScalarType(dp_datum[0]),
    //                                   MyScalarType(dp_datum[1]),
    //                                   MyScalarType(dp_datum[2]));
    // auto paraboloid =
    //     ParaboloidBase<MyScalarType>(datum, frame,
    //     MyScalarType(dp_aligned.a()),
    //                                  MyScalarType(dp_aligned.b()));
    // // auto lower_corner = MyPtType(MyScalarType(lower_corner_m[0]),
    // //                              MyScalarType(lower_corner_m[1]),
    // //                              MyScalarType(lower_corner_m[2]));
    // // auto upper_corner = MyPtType(MyScalarType(upper_corner_m[0]),
    // //                              MyScalarType(upper_corner_m[1]),
    // //                              MyScalarType(upper_corner_m[2]));
    // // auto cube = StoredRectangularCuboid<MyPtType>::fromBoundingPts(
    // //     lower_corner, upper_corner);
    // auto cube = StoredRectangularCuboid<MyPtType>::fromOtherPolytope(cell_m);
    // auto moments = getVolumeMoments<MyMomentType>(cube, paraboloid);
    // const MyScalarType& volume = moments.volume();

    // // for (int parameter = 0; parameter < kColumns; ++parameter) {
    // //   (*a_jacobian_transpose)(parameter, 0) =
    // //       volume.gradient().getGrad()[parameter];
    // //   (*a_jacobian_transpose)(parameter, 1) =
    // //       moments.centroid()[0].gradient().getGrad()(parameter);
    // //   (*a_jacobian_transpose)(parameter, 2) =
    // //       moments.centroid()[1].gradient().getGrad()(parameter);
    // //   (*a_jacobian_transpose)(parameter, 3) =
    // //       moments.centroid()[2].gradient().getGrad()(parameter);
    // // }
    // (*a_jacobian_transpose).Zero();

    // (*a_jacobian_transpose)(0, 0) = volume.gradient().getGradA() / m0_ref_m;
    // (*a_jacobian_transpose)(1, 0) = volume.gradient().getGradB() / m0_ref_m;
    // (*a_jacobian_transpose)(4, 0) = volume.gradient().getGradTz() / m0_ref_m;
    // (*a_jacobian_transpose)(7, 0) = volume.gradient().getGradRz() / m0_ref_m;
    // for (int i = 0; i < 3; i++) {
    //   (*a_jacobian_transpose)(0, i + 1) =
    //       moments.centroid()[i].gradient().getGradA() / m1_ref_m;
    //   (*a_jacobian_transpose)(1, i + 1) =
    //       moments.centroid()[i].gradient().getGradB() / m1_ref_m;
    //   (*a_jacobian_transpose)(4, i + 1) =
    //       moments.centroid()[i].gradient().getGradTz() / m1_ref_m;
    //   (*a_jacobian_transpose)(7, i + 1) =
    //       moments.centroid()[i].gradient().getGradRz() / m1_ref_m;
    // }

    // // Calculate Tranpose(Jacobian) * Jacobian
    // *a_jacTjac = (*a_jacobian_transpose) *
    // (a_jacobian_transpose->transpose());
  }

  void updateBestGuess(void) {
    best_values_m = guess_values_m;
    best_reconstruction_m = guess_reconstruction_m;
  }

  Eigen::Matrix<double, kRows, 1> calculateChangeInGuess(void) {
    return guess_values_m - best_values_m;
  }

  void increaseLambda(double* a_lambda) {
    (*a_lambda) *= optimization_behavior_m.lambda_increase;
  }

  void decreaseLambda(double* a_lambda) {
    (*a_lambda) *= optimization_behavior_m.lambda_decrease;
  }

  bool errorTooHigh(const double a_error) {
    return a_error > optimization_behavior_m.acceptable_error;
  }

  bool iterationTooHigh(const UnsignedIndex_t a_iteration) {
    return a_iteration > optimization_behavior_m.maximum_iterations;
  }

  bool shouldComputeJacobian(const UnsignedIndex_t a_iteration,
                             const UnsignedIndex_t a_last_jacobian) {
    return a_iteration - a_last_jacobian >
           optimization_behavior_m.delay_jacobian_amount;
  }

  bool minimumReached(const Eigen::Matrix<double, kColumns, 1> a_delta) {
    return a_delta.squaredNorm() < 1.0e-30;  // 1.0e-16;
  }

  Paraboloid& getBestReconstruction(void) { return best_reconstruction_m; }

  void clipChange(Eigen::Matrix<double, kColumns, 1>* a_delta) {
    double scale = 1.0;
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      if (std::fabs((*a_delta)(2 + d)) > 0.5 * ref_length_m) {
        scale =
            std::min(scale, 0.5 * ref_length_m / std::fabs((*a_delta)(2 + d)));
      }
      if (std::fabs((*a_delta)(5 + d)) > 0.1 * M_PI) {
        scale = std::min(scale, 0.1 * M_PI / std::fabs((*a_delta)(5 + d)));
      }
    }
    if (std::fabs((*a_delta)(0)) > 0.1 / ref_length_m) {
      scale = std::min(scale, 0.1 / ref_length_m / std::fabs((*a_delta)(0)));
    }
    if (std::fabs((*a_delta)(1)) > 0.1 / ref_length_m) {
      scale = std::min(scale, 0.1 / ref_length_m / std::fabs((*a_delta)(1)));
    }
    if (scale != 1.0) {
      for (UnsignedIndex_t d = 0; d < kColumns; ++d) {
        (*a_delta)(d) *= scale;
      }
    }
  }

  Eigen::Matrix<double, kColumns, 1>& getDefaultInitialDelta(void) {
    for (UnsignedIndex_t i = 0; i < kColumns; ++i) {
      initial_delta_m(i) = sqrt(DBL_EPSILON);
    }
    return initial_delta_m;
  }

  ~Paraboloid_MOF_Object() = default;

 private:
  OptimizationBehavior optimization_behavior_m;
  RectangularCuboid cell_m;
  double m0_ref_a_m, m1_ref_a_m, m2_ref_a_m;
  double m0_ref_b_m, m1_ref_b_m, m2_ref_b_m;
  double ref_length_m;
  Eigen::Matrix<double, kRows, 1> correct_values_m;
  Eigen::Matrix<double, kRows, 1> guess_values_m;
  Eigen::Matrix<double, kRows, 1> best_values_m;
  Eigen::Matrix<double, kColumns, 1> initial_delta_m;
  Paraboloid guess_reconstruction_m, best_reconstruction_m;
  Pt lower_corner_m, upper_corner_m;
};

// PolynomialLM poly;
// poly.setup(0.0);
// Eigen::Matrix<double, 1, 1> initial_delta;
// initial_delta(0) = 0.1;
// LevenbergMarquardt<PolynomialLM, 1, 1> lm_solver;
// lm_solver.solve(&poly, initial_delta);
// EXPECT_NEAR(poly.getSolution(), 4.0, 1.0e-4)
//     << "Reason for exit: " << lm_solver.getReason() << '\n';

TEST(ParaboloidMOF, MOF) {
  // using ScalarType = Quad_t;
  // using Pt = PtBase<ScalarType>;
  // using Normal = NormalBase<ScalarType>;
  // using ReferenceFrame = ReferenceFrameBase<ScalarType>;
  // using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;
  // using Paraboloid = ParaboloidBase<ScalarType>;

  // Construct NxNxN mesh
  int n_cells_diam = 10;
  double radius = static_cast<double>(n_cells_diam) / 2.0;
  int n_cells = 2 + n_cells_diam;
  constexpr const int number_of_ghost_cells = 0;
  const double dx = 2.0 * radius / static_cast<double>(n_cells_diam);
  const double lx = dx * static_cast<double>(n_cells);
  Pt lower_domain(-lx / 2.0, -lx / 2.0, -lx / 2.0);
  Pt upper_domain(lx / 2.0, lx / 2.0, lx / 2.0);
  BasicMesh mesh(n_cells, n_cells, n_cells, number_of_ghost_cells);
  mesh.setCellBoundaries(lower_domain, upper_domain);

  // Initialize VF from sphere
  std::random_device rd;
  std::mt19937_64 eng(rd());  // rd());
  std::uniform_real_distribution<double> random_tr(-0.5, 0.5);

  auto center = Pt(1.0e-6, 1.0e-6, 1.0e-6);
  // auto center = Pt(random_tr(eng), random_tr(eng), random_tr(eng));
  double curvature = 1.0 / radius;

  Data<double> liquid_vf(mesh);
  Data<SeparatedMoments<GeneralMoments3D<2>>> liquid_moments(mesh);
  Data<Paraboloid> liquid_gas_interface(mesh);

  int sub_div = 10;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        double sub_dx =
            (mesh.x(i + 1) - mesh.x(i)) / static_cast<double>(sub_div);
        double sub_dy =
            (mesh.y(j + 1) - mesh.y(j)) / static_cast<double>(sub_div);
        double sub_dz =
            (mesh.z(k + 1) - mesh.z(k)) / static_cast<double>(sub_div);
        liquid_moments(i, j, k) =
            SeparatedMoments<GeneralMoments3D<2>>::fromScalarConstant(0.0);
        const Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
        const auto cell_centroid = Pt(0, 0, 0);
        // Pt(0.5 * (lower_cell_pt + upper_cell_pt));
        for (int kk = 0; kk < sub_div; ++kk) {
          for (int jj = 0; jj < sub_div; ++jj) {
            for (int ii = 0; ii < sub_div; ++ii) {
              const Pt lower_subcell_pt(
                  mesh.x(i) + static_cast<double>(ii) * sub_dx,
                  mesh.y(j) + static_cast<double>(jj) * sub_dy,
                  mesh.z(k) + static_cast<double>(kk) * sub_dz);
              const Pt upper_subcell_pt(
                  mesh.x(i) + static_cast<double>(ii + 1) * sub_dx,
                  mesh.y(j) + static_cast<double>(jj + 1) * sub_dy,
                  mesh.z(k) + static_cast<double>(kk + 1) * sub_dz);
              const auto sub_cell = RectangularCuboid::fromBoundingPts(
                  Pt(lower_subcell_pt - cell_centroid),
                  Pt(upper_subcell_pt - cell_centroid));
              Normal normal =
                  0.5 * Pt(lower_subcell_pt + upper_subcell_pt) - center;
              normal.normalize();
              int largest_dir = 0;
              if (std::fabs(normal[largest_dir]) < std::fabs(normal[1]))
                largest_dir = 1;
              if (std::fabs(normal[largest_dir]) < std::fabs(normal[2]))
                largest_dir = 2;
              ReferenceFrame frame;
              if (largest_dir == 0)
                frame[0] = crossProduct(normal, Normal(0.0, 1.0, 0.0));
              else if (largest_dir == 0)
                frame[0] = crossProduct(normal, Normal(0.0, 0.0, 1.0));
              else
                frame[0] = crossProduct(normal, Normal(1.0, 0.0, 0.0));
              frame[0].normalize();
              frame[1] = crossProduct(normal, frame[0]);
              frame[2] = normal;
              Paraboloid paraboloid(
                  Pt(center + radius * normal - cell_centroid), frame,
                  0.5 * curvature, 0.5 * curvature);
              liquid_moments(i, j, k) +=
                  getVolumeMoments<SeparatedMoments<GeneralMoments3D<2>>>(
                      sub_cell, paraboloid);
              // auto moments =
              //     getVolumeMoments<VolumeMoments>(sub_cell, paraboloid);
              // liquid_moments(i, j, k)[0] += moments.volume();
              // liquid_moments(i, j, k)[1] += moments.centroid()[0];
              // liquid_moments(i, j, k)[2] += moments.centroid()[1];
              // liquid_moments(i, j, k)[3] += moments.centroid()[2];
            }
          }
        }
        {
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              Pt(lower_cell_pt - cell_centroid),
              Pt(upper_cell_pt - cell_centroid));
          liquid_vf(i, j, k) =
              liquid_moments(i, j, k)[0].volume() / cell.calculateVolume();
          liquid_vf(i, j, k) = std::max(liquid_vf(i, j, k), 0.0);
          liquid_vf(i, j, k) = std::min(liquid_vf(i, j, k), 1.0);

          const auto centroid = Pt(
              cell_centroid + Pt(liquid_moments(i, j, k)[0][1],
                                 liquid_moments(i, j, k)[0][2],
                                 liquid_moments(i, j, k)[0][3]) /
                                  safelyEpsilon(liquid_moments(i, j, k)[0][0]));

          Normal normal = centroid - center;
          normal.normalize();
          if (i == 1 && j == 1 & k == 1) normal[0] += 0.0;
          normal.normalize();
          int largest_dir = 0;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[1]))
            largest_dir = 1;
          if (std::fabs(normal[largest_dir]) < std::fabs(normal[2]))
            largest_dir = 2;
          ReferenceFrame frame;
          if (largest_dir == 0)
            frame[0] = crossProduct(normal, Normal(0.0, 1.0, 0.0));
          else if (largest_dir == 0)
            frame[0] = crossProduct(normal, Normal(0.0, 0.0, 1.0));
          else
            frame[0] = crossProduct(normal, Normal(1.0, 0.0, 0.0));
          frame[0].normalize();
          frame[1] = crossProduct(normal, frame[0]);
          frame[2] = normal;
          double a = 0.5 * curvature;
          double b = 0.5 * curvature;
          if (i == 1 && j == 1 & k == 1) {
            a = 0.5 * curvature;
            b = 0.5 * curvature;
          }
          // a += 0.4 * curvature;
          // b -= 0.4 * curvature;
          a = 1.0e-6;
          b = 1.0e-6;
          // Paraboloid paraboloid(center + radius * normal, frame, a, b);
          Paraboloid paraboloid(center + radius * normal - cell_centroid, frame,
                                a, b);
          ProgressiveDistanceSolverParaboloid<RectangularCuboid> solver(
              cell, liquid_vf(i, j, k), 1.0e-14, paraboloid);
          // Paraboloid new_paraboloid(
          //     Pt(center + (radius + solver.getDistance()) * normal), frame,
          //     a, b);
          Paraboloid new_paraboloid(
              Pt(center + radius * normal - cell_centroid +
                 (solver.getDistance()) * normal),
              frame, a, b);
          liquid_gas_interface(i, j, k) = new_paraboloid;
        }
      }
    }
  }

  // Initialize folders/mesh for very simple I/O
  int viz_output = 0;
  double time = 0.0;
  VTKOutput vtk_io("viz_out", "viz", mesh);
  vtk_io.addData("VOF", liquid_vf);
  vtk_io.writeVTKFile(time);
  std::vector<ParametrizedSurfaceOutput> surfaces;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (liquid_vf(i, j, k) >= global_constants::VF_LOW &&
            liquid_vf(i, j, k) <= global_constants::VF_HIGH) {
          const Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
          const auto cell_centroid = Pt(0, 0, 0);
          // Pt(0.5 * (lower_cell_pt + upper_cell_pt));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          const auto aligned_paraboloid =
              liquid_gas_interface(i, j, k).getAlignedParaboloid();
          const auto datum =
              Pt(cell_centroid + liquid_gas_interface(i, j, k).getDatum());
          const auto frame = liquid_gas_interface(i, j, k).getReferenceFrame();
          const auto paraboloid = Paraboloid(
              datum, frame, aligned_paraboloid.a(), aligned_paraboloid.b());
          auto volume_and_surface = getVolumeMoments<
              AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(cell,
                                                                   paraboloid);
          auto& surface = volume_and_surface.getSurface();
          surface.setLengthScale(std::pow(cell.calculateVolume(), 1.0 / 3.0) /
                                 10.0);
          surfaces.push_back(surface);
        }
      }
    }
  }
  vtk_io.writeVTKInterface(time, surfaces, true);
  ++viz_output;
  std::cout << "Initial surface printed" << std::endl;
  surfaces.clear();

  global_constants::VF_LOW = 1.0e-8;
  global_constants::VF_HIGH = 1.0 - global_constants::VF_LOW;

  double total_surface = 0.0;
  double avg_mean_curv = 0.0;
  double rms_curv_error = 0.0, max_curv_error = 0.0;
  double rms_counter = 0.0;
  double exact_curv = 2.0 / radius;
  for (int i = mesh.imin() + 1; i <= mesh.imax() - 1; ++i) {
    for (int j = mesh.jmin() + 1; j <= mesh.jmax() - 1; ++j) {
      for (int k = mesh.kmin() + 1; k <= mesh.kmax() - 1; ++k) {
        if (liquid_vf(i, j, k) >= global_constants::VF_LOW &&
            liquid_vf(i, j, k) <= global_constants::VF_HIGH) {
          Paraboloid paraboloid_guess = liquid_gas_interface(i, j, k);
          // // SUPERCELL
          // const Pt lower_cell_pt(mesh.x(i - 1), mesh.y(j - 1), mesh.z(k -
          // 1)); const Pt upper_cell_pt(mesh.x(i + 2), mesh.y(j + 2), mesh.z(k
          // + 2)); const auto cell_centroid = Pt(0, 0, 0); auto moments =
          //     SeparatedMoments<GeneralMoments3D<2>>::fromScalarConstant(0.0);
          // for (int ii = i - 1; ii < i + 2; ++ii) {
          //   for (int jj = j - 1; jj < j + 2; ++jj) {
          //     for (int kk = k - 1; kk < k + 2; ++kk) {
          //       moments += liquid_moments(ii, jj, kk);
          //     }
          //   }
          // }

          // LOCAL CELL
          const Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
          const auto cell_centroid = Pt(0, 0, 0);
          auto moments = liquid_moments(i, j, k);

          const auto cell = RectangularCuboid::fromBoundingPts(
              Pt(lower_cell_pt - cell_centroid),
              Pt(upper_cell_pt - cell_centroid));

          Paraboloid_MOF_Object<20, 8> mof_object(
              cell, Pt(lower_cell_pt - cell_centroid),
              Pt(upper_cell_pt - cell_centroid), moments, paraboloid_guess);
          LevenbergMarquardtNew<Paraboloid_MOF_Object<20, 8>, 20, 8> solver;
          // std::cout << "Solving... " << std::endl;
          solver.solve(&mof_object, mof_object.getDefaultInitialDelta());
          // std::cout << "Solved in   : " << solver.getIterationCount()
          //           << " iterations with error = "
          //           << sqrt(mof_object.calculateScalarError()) <<
          // std::endl;

          Paraboloid paraboloid = mof_object.getBestReconstruction();
          auto aligned_paraboloid = paraboloid.getAlignedParaboloid();
          auto datum = paraboloid.getDatum();
          auto frame = paraboloid.getReferenceFrame();
          // auto normal111 = Normal(1, 1, 1);
          // normal111.normalize();
          // if (magnitude(frame[2] - normal111) < 1.0e-6) {
          //   std::cout << "Normal  = " << frame[2] << std::endl;
          //   std::cout << "Moments = " << liquid_moments(i, j, k) <<
          //   std::endl;
          // }

          const Pt lower_cell_pt_new(mesh.x(i), mesh.y(j), mesh.z(k));
          const Pt upper_cell_pt_new(mesh.x(i + 1), mesh.y(j + 1),
                                     mesh.z(k + 1));
          const auto new_cell = RectangularCuboid::fromBoundingPts(
              Pt(lower_cell_pt_new - cell_centroid),
              Pt(lower_cell_pt_new - cell_centroid));
          ProgressiveDistanceSolverParaboloid<RectangularCuboid>
              solver_distance(new_cell, liquid_vf(i, j, k), 1.0e-14,
                              paraboloid);
          Paraboloid new_paraboloid(
              Pt(cell_centroid + datum +
                 solver_distance.getDistance() * frame[2]),
              frame, aligned_paraboloid.a(), aligned_paraboloid.b());

          const auto new_new_cell = RectangularCuboid::fromBoundingPts(
              Pt(lower_cell_pt_new), Pt(upper_cell_pt_new));
          auto volume_and_surface = getVolumeMoments<
              AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(
              new_new_cell, new_paraboloid);
          auto& surface = volume_and_surface.getSurface();
          auto area = surface.getSurfaceArea();
          auto curv_int = surface.getMeanCurvatureIntegral();
          total_surface += area;
          avg_mean_curv += curv_int;
          auto curv = curv_int / safelyEpsilon(area);
          auto curv_error = std::fabs(curv - exact_curv) / exact_curv;
          // std::cout << liquid_vf(i, j, k) << " " << curv_error << std::endl;
          max_curv_error = std::max(max_curv_error, curv_error);
          rms_curv_error += curv_error * curv_error;
          rms_counter += 1.0;
          surface.setLengthScale(std::pow(cell.calculateVolume(), 1.0 / 3.0) /
                                 15.0);
          surfaces.push_back(surface);
        }
      }
    }
  }

  rms_curv_error = std::sqrt(rms_curv_error / rms_counter);

  avg_mean_curv /= total_surface;
  total_surface = std::sqrt(total_surface);
  std::cout << "Cells / diameter = " << n_cells_diam << std::endl;
  std::cout << "Surface area = " << total_surface << " instead of "
            << std::sqrt(4.0 * M_PI * radius * radius) << " (rel error = "
            << std::fabs(total_surface -
                         std::sqrt(4.0 * M_PI * radius * radius)) /
                   std::sqrt(4.0 * M_PI * radius * radius)
            << ")" << std::endl;
  std::cout << "Avg mean curvature = " << avg_mean_curv << " instead of "
            << 2.0 / radius << " (rel error = "
            << std::fabs(avg_mean_curv - 2.0 / radius) / (2.0 / radius) << ")"
            << std::endl;

  std::cout << std::scientific << std::setprecision(6)
            << "CURV ERROR:\nRMS = " << rms_curv_error
            << "\nMAX = " << max_curv_error << std::endl;

  std::cout << "Printing final solution " << std::endl;
  vtk_io.writeVTKInterface(time, surfaces, true);
}

}  // namespace IRL
