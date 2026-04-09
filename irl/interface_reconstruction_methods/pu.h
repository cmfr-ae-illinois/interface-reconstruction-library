// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_H_

#include <Eigen/Dense>
#include <cassert>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

class PU {
 public:
  /// \brief Default constructor.
  PU(void) = default;

  /// \brief Solve for the paraboloid using the provided neighborhood and delta.
  Paraboloid solve(const PUNeighborhood* a_neighborhood_pointer,
                   const double a_delta);

  /// \brief default destructor
  ~PU(void) = default;

 private:
  /// \brief solve the system for reconstruction
  Paraboloid solve(void);

  /// \brief PU surface and gradient
  std::pair<double, Eigen::Vector3d> getPUAndGrad(const Pt& a_pt);

  /// \brief PU surface and gradient and Hessian
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> getPUAndGradAndHessian(
      const Pt& a_pt);

  /// \brief Projecting point on PU surface
  Pt projectPointonPU(const Pt& a_pt);

  /// \brief storing stencil information
  const PUNeighborhood* neighborhood_m;
  /// \brief kernal radius for PU
  double delta_m;
};

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_H_
