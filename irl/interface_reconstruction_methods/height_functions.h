// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_HEIGHT_FUNCTIONS_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_HEIGHT_FUNCTIONS_H_

#include <Eigen/Dense>
#include <cassert>
#include <cstdint>

#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/defined_types.h"
#include "irl/variant_reconstruction/separator_variant.h"

namespace IRL {

/// \brief Outcome of a height-function reconstruction.
enum class HeightFunctionStatus : std::uint8_t {
  NOT_SOLVED,
  RETURNED_PLANE,
  RETURNED_PARABOLOID_FROM_PARTIAL_STENCIL,
  RETURNED_PARABOLOID_FROM_FULL_STENCIL
};

template <int Nx, int Ny, int Nz, class SeparatorVariantType = SeparatorVariant>
class HeightFunctions_3D {
 public:
  /// \brief Default constructor.
  HeightFunctions_3D(void) = default;
  HeightFunctions_3D(const double a_tol) : tol_m(a_tol) {}

  Normal lvira_normal(
      const std::array<std::array<std::array<double, Nz>, Ny>, Nx>& a_vfrac);

  /// \brief Solve the system for the reconstruction
  SeparatorVariantType solve(
      const std::array<std::array<std::array<double, Nz>, Ny>, Nx>& a_vfrac,
      const Pt a_stencil_center, const double a_dx,
      const Normal& a_normal = Normal(0, 0, 0));

  /// \brief Return the outcome of the last call to solve.
  HeightFunctionStatus status(void) const { return status_m; }

  const int Nhx(void) const { return Nhx_m; }
  const int Nhy(void) const { return Nhy_m; }
  const int nvalid(void) const {
    int count = 0;
    for (int i = 0; i < Nhx_m; ++i) {
      for (int j = 0; j < Nhy_m; ++j) {
        if (valid_m(i, j)) {
          ++count;
        }
      }
    }
    return count;
  }

  /// \brief Default destructor.
  ~HeightFunctions_3D(void) = default;

 private:
  int compute_coeffs(const double dx);
  int compute_heights(
      const std::array<std::array<std::array<double, Nz>, Ny>, Nx>& a_vfrac,
      const int a_dir, const int a_sign);

  int Nhx_m, Nhy_m;
  double tol_m = 1.0e-12;
  Eigen::Matrix<double, 6, 1> coefficients_m;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> heights_m;
  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> valid_m;
  /// \brief Outcome of the last call to solve
  HeightFunctionStatus status_m = HeightFunctionStatus::NOT_SOLVED;
};

}  // namespace IRL

#include "irl/interface_reconstruction_methods/height_functions.tpp"

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_HEIGHT_FUNCTIONS_H_
