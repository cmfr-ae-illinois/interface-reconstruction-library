// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_PARABOLOID_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_PARABOLOID_TPP_

#include <cmath>
#include <limits>
#include <tuple>

#include <Eigen/Dense>

namespace IRL {

template <class CellType>
PUParaboloid<CellType>::PUParaboloid(
    const PUNeighborhood<CellType>& a_neighborhood, const double a_kernel_size,
    const double a_dx)
    : PU<CellType>(a_neighborhood, a_kernel_size), dx_m(a_dx) {}

template <class CellType>
PUParaboloid<CellType>::PUParaboloid(
    const PUNeighborhood<CellType>& a_neighborhood, const double a_kernel_size)
    : PU<CellType>(a_neighborhood, a_kernel_size) {}

template <class CellType>
void PUParaboloid<CellType>::setDx(const double a_dx) {
  dx_m = a_dx;
}

template <class CellType>
double PUParaboloid<CellType>::getDx(void) const {
  return dx_m;
}

template <class CellType>
Paraboloid PUParaboloid<CellType>::solve(
    const PUNeighborhood<CellType>& a_neighborhood, const double a_kernel_size,
    const double a_dx) {
  this->setNeighborhood(a_neighborhood);
  this->setKernelSize(a_kernel_size);
  dx_m = a_dx;

  return this->solve();
}

template <class CellType>
Paraboloid PUParaboloid<CellType>::makeInvalidParaboloid(void) const {
  Paraboloid paraboloid;

  const double infinity = std::numeric_limits<double>::infinity();
  paraboloid.setDatum(Pt(infinity, infinity, infinity));

  auto coefficients = paraboloid.getAlignedParaboloid();
  coefficients.a() = 20.0 / dx_m;
  coefficients.b() = 20.0 / dx_m;
  paraboloid.setAlignedParaboloid(coefficients);

  return paraboloid;
}

template <class CellType>
Paraboloid PUParaboloid<CellType>::solve(void) {
  const auto& neighborhood = this->getNeighborhood();

  if (neighborhood.size() == 0) {
    return this->makeInvalidParaboloid();
  }

  const auto center_index = neighborhood.getCenterOfStencil();

  if (center_index >= neighborhood.size()) {
    return this->makeInvalidParaboloid();
  }

  // centroid of target interface
  const Pt center_centroid = neighborhood.getCentroid(center_index);

  if (!std::isfinite(center_centroid[0]) ||
      !std::isfinite(center_centroid[1]) ||
      !std::isfinite(center_centroid[2])) {
    return this->makeInvalidParaboloid();
  }

  // prtojecting centroid onto zero-level set of PU surface
  bool projection_success = false;
  const Pt projected_centroid =
      this->projectOntoPU(center_centroid, dx_m, projection_success);

  if (!projection_success || !std::isfinite(projected_centroid[0]) ||
      !std::isfinite(projected_centroid[1]) ||
      !std::isfinite(projected_centroid[2])) {
    return this->makeInvalidParaboloid();
  }

  // Evaluate PU value, gradient, and Hessian at the projected point.
  const auto pu_derivatives = this->getPUGradAndHess(projected_centroid);

  const double pu_value = std::get<0>(pu_derivatives);
  const Eigen::Vector3d& gradient = std::get<1>(pu_derivatives);
  const Eigen::Matrix3d& hessian = std::get<2>(pu_derivatives);

  if (!std::isfinite(pu_value) || !gradient.allFinite() ||
      !hessian.allFinite()) {
    return this->makeInvalidParaboloid();
  }

  const double gradient_squared_norm = gradient.squaredNorm();

  if (!std::isfinite(gradient_squared_norm) ||
      gradient_squared_norm <= std::numeric_limits<double>::epsilon()) {
    return this->makeInvalidParaboloid();
  }

  // construct paraboloid from derivatives at the projected point
  Paraboloid paraboloid =
      Paraboloid::fromDerivatives(projected_centroid, gradient, hessian);

  // Validate the generated paraboloid before returning it.
  const Pt& datum = paraboloid.getDatum();
  const auto& aligned_paraboloid = paraboloid.getAlignedParaboloid();

  if (!std::isfinite(datum[0]) || !std::isfinite(datum[1]) ||
      !std::isfinite(datum[2]) || !std::isfinite(aligned_paraboloid.a()) ||
      !std::isfinite(aligned_paraboloid.b())) {
    return this->makeInvalidParaboloid();
  }

  return paraboloid;
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_PARABOLOID_TPP_