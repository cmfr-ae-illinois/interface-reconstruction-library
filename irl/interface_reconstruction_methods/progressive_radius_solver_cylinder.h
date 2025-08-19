// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_RADIUS_SOLVER_CYLINDER_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_RADIUS_SOLVER_CYLINDER_H_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/generic_cutting/quadratic_intersection/moment_contributions.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/helper.h"
#include "irl/helpers/mymath.h"
#include "irl/parameters/constants.h"

namespace IRL {

/// \brief Volume conserving distance-finding routine for
/// single paraboloid reconstructions.
///
/// \param[in] volume_fraction_m Volume fraction to recreate
/// \param[in] volume_fraction_tolerance_m Tolerance to recreate
/// `a_volume_fraction` within
/// \param[in] reconstruction_m A copy of the reconstruction
/// to find distances for.
/// \param[out] radius_m Correct distance to plane is stored
/// in distances_m after construction of the object.
template <class CellType>
class ProgressiveRadiusSolverCylinder {
  /// \brief Max number of iterations for the secant Solver
  static constexpr UnsignedIndex_t max_iter_m = {80};
  static constexpr UnsignedIndex_t max_bisection_iter = {80};

 public:
  /// \brief Constructor that initializes the class for optimization
  /// and solves for the radius.

  /// \brief Default constructor
  ProgressiveRadiusSolverCylinder(void) = default;

  /// \brief Constructor that loads in the neccessary information to
  /// the class and then solves for volume-conserving radius of a cylinder
  /// intersection.
  /// \param[in] a_volume_fraction Volume fraction to match.
  /// \param[in] a_volume_fraction_tolerance Tolerance to recreate volume
  /// fraction within.
  /// \param[in] a_reconstruction Reconstruction to find radius for.
  ProgressiveRadiusSolverCylinder(const CellType& a_cell,
                                  const double a_volume_fraction,
                                  const double a_volume_fraction_tolerance,
                                  const Cylinder& a_reconstruction);

  /// \brief Reinitialize solver and solve for radius
  /// \param[in] a_volume_fraction Volume fraction to match.
  /// \param[in] a_volume_fraction_tolerance Tolerance to recreate volume
  /// fraction within.
  /// \param[in] a_reconstruction Reconstruction to find distances for.
  void solve(const CellType& a_cell, const double a_volume_fraction,
             const double a_volume_fraction_tolerance,
             const Cylinder& a_reconstruction);

  /// \brief Return pointer to radius_m to be used for changing a
  /// reconstruction.
  Cylinder getCylinder(void);

  /// \brief Default destructor
  ~ProgressiveRadiusSolverCylinder(void) = default;

 private:
  /// \brief Main driving function that solves for distances
  inline void solveForRadius(const CellType& a_cell);

  /// \brief Starting volume of the cell
  double initial_cell_volume_m;
  /// \brief Target volume fraction to match
  double target_volume_fraction_m;
  /// \brief Tolerance to match volume fraction within
  double volume_fraction_tolerance_m;
  /// \brief The reconstruction the distance is being found for
  Cylinder reconstruction_m;
  /// \brief Radius that satisfy tolerance.
  double r_m;
  /// \brief Stretching that satisfy tolerance.
  double b_m;
};

}  // namespace IRL

#include "irl/interface_reconstruction_methods/progressive_radius_solver_cylinder.tpp"

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_RADIUS_SOLVER_CYLINDER_H_
