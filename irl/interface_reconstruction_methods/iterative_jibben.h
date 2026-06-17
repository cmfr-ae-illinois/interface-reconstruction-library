// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_JIBBEN_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_JIBBEN_H_

#include <Eigen/Dense>
#include <cassert>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/interface_reconstruction_methods/jibben_neighborhood.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

class iJibben_3D {
 public:
  /// \brief Default constructor.
  iJibben_3D(void) = default;

  /// \brief constructor with neighborhood
  iJibben_3D(const JibbenNeighborhood* a_neighborhood_pointer);

  /// \brief Get the paraboloid coefficients by solving the least-squares
  /// system.
  void getParaboloidCoefficients(void);

  /// \brief computing mean curvature of the paraboloid
  double computeMeanCurvature(const double& xi, const double& eta) const;

  /// \brief computing normal vector
  IRL::Normal computeNormal(const double& xi, const double& eta) const;

  /// \brief Make a paraboloid in the global principal frame.
  Paraboloid makeParaboloid(const IRL::Pt& a_local_datum,
                            const IRL::ReferenceFrame& a_local_frame) const;
  Paraboloid makeParaboloid2(
      const JibbenNeighborhood* a_neighborhood_pointer) const;

  std::pair<double, Normal> computeAveragedCurvatureAndNormal() const;

  /// \brief Default destructor.
  ~iJibben_3D(void) = default;

 private:
  // /// \brief Solve the system for the reconstruction
  // Paraboloid solve(void);

  /// \brief Storage of the stencil information
  const JibbenNeighborhood* neighborhood_m;
  /// \brief Paraboloid coefficients
  std::array<double, 6> coefficients_m;
};

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_JIBBEN_H_
