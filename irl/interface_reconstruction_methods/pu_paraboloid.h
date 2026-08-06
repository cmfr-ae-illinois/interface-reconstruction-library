// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_PARABOLOID_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_PARABOLOID_H_

#include <limits>

#include "irl/interface_reconstruction_methods/pu.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace IRL {

// \brief Partition-of-unity paraboloid reconstruction
//
// This class uses the implicit PU surface supplied by the PU base class.
// The centroid of the center interface is first projected onto the PU
// surface. The gradient and Hessian at the projected point are then used
// to construct a local paraboloid
template <class CellType>
class PUParaboloid : public PU<CellType> {
 public:
  // Default constructor
  PUParaboloid(void) = default;

  // Construct from a neighborhood, kernel size, and cell spacing
  PUParaboloid(const PUNeighborhood<CellType>& a_neighborhood,
               const double a_kernel_size, const double a_dx);

  // Construct from a neighborhood and kernel size
  ///
  // The grid spacing must be supplied later through solve() or setDx()
  PUParaboloid(const PUNeighborhood<CellType>& a_neighborhood,
               const double a_kernel_size);

  // Reconstruct a paraboloid using the currently stored neighborhood
  // kernel size, and grid spacing
  Paraboloid solve(void);

  // Replace the neighborhood and reconstruct
  Paraboloid solve(const PUNeighborhood<CellType>& a_neighborhood,
                   const double a_kernel_size, const double a_dx);

  // Set the grid spacing
  void setDx(const double a_dx);

  // Return the currently stored grid spacing
  double getDx(void) const;

 private:
  // Construct a clearly invalid paraboloid used to indicate failure
  Paraboloid makeInvalidParaboloid(void) const;

  // cell spacing
  double dx_m = 0.0;
};

}  // namespace IRL

#include "irl/interface_reconstruction_methods/pu_paraboloid.tpp"

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_PARABOLOID_H_