// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_JIBBEN_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_JIBBEN_H_

#include <Eigen/Dense>
#include <cassert>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/interface_reconstruction_methods/jibben_neighborhood.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

class Jibben_3D {
 public:
  /// \brief Default constructor.
  Jibben_3D(void) = default;

  /// \brief Solve the system for the reconstruction
  Paraboloid solve(const JibbenNeighborhood* a_neighborhood_pointer,
                   const double a_delta = -1.0);

  /// \brief Getting fitting error
  const double getFittingError(void) const;

  /// \brief Get volume error
  const double getVolumeError(const double& dx) const;

  /// \brief Get volume error squared
  const double getVolumeErrorSquared(const double& dx) const;

  const double getVolumeErrorSquared_w1(const double& dx) const;

  const double getVolumeErrorSquared_w2(const double& dx) const;

  /// \brief Default destructor.
  ~Jibben_3D(void) = default;

 private:
  /// \brief Solve the system for the reconstruction
  Paraboloid solve(void);

  /// \brief Storage of the stencil information
  const JibbenNeighborhood* neighborhood_m;
  /// \brief Error resulting from the fit
  double error_m;
  /// \brief Weighting function radius
  double delta_m;
  /// \brief Paraboloid coefficients
  std::array<double, 6> coefficients_m;
};

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_JIBBEN_H_
