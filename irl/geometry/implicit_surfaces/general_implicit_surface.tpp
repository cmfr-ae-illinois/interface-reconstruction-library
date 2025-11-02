// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_IMPLICIT_SURFACES_GENERAL_IMPLICIT_SURFACE_TPP_
#define IRL_GEOMETRY_IMPLICIT_SURFACES_GENERAL_IMPLICIT_SURFACE_TPP_

#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"

namespace IRL {

template <class ScalarType, size_t MAX_REFINE_LEVEL>
ScalarType GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL>::meanCurvature(
    const ScalarType& x, const ScalarType& y, const ScalarType& z) const {
  const Vec3 g = this->gradF(x, y, z);
  const ScalarType G = g.squaredNorm();
  if (G < 1e-30) return std::numeric_limits<ScalarType>::quiet_NaN();

  const Mat3 Hf = this->hessF(x, y, z);
  const ScalarType trH = Hf.trace();
  const ScalarType gtHg = g.transpose() * (Hf * g);

  return (-gtHg + G * trH) / std::pow(G, 1.5);
}

template <class ScalarType, size_t MAX_REFINE_LEVEL>
ScalarType
GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL>::gaussianCurvature(
    const ScalarType& x, const ScalarType& y, const ScalarType& z) const {
  const Vec3 g = this->gradF(x, y, z);
  const ScalarType norm_g = g.norm();
  const ScalarType eps = ScalarType(1e-30);
  if (norm_g <= eps) return std::numeric_limits<ScalarType>::quiet_NaN();

  const Mat3 Hf = this->hessF(x, y, z);

  // cofactor matrix
  Mat3 cof;
  cof(0, 0) = Hf(1, 1) * Hf(2, 2) - Hf(1, 2) * Hf(2, 1);
  cof(0, 1) = Hf(0, 2) * Hf(2, 1) - Hf(0, 1) * Hf(2, 2);
  cof(0, 2) = Hf(0, 1) * Hf(1, 2) - Hf(0, 2) * Hf(1, 1);

  cof(1, 0) = Hf(1, 2) * Hf(2, 0) - Hf(1, 0) * Hf(2, 2);
  cof(1, 1) = Hf(0, 0) * Hf(2, 2) - Hf(0, 2) * Hf(2, 0);
  cof(1, 2) = Hf(0, 2) * Hf(1, 0) - Hf(0, 0) * Hf(1, 2);

  cof(2, 0) = Hf(1, 0) * Hf(2, 1) - Hf(1, 1) * Hf(2, 0);
  cof(2, 1) = Hf(0, 1) * Hf(2, 0) - Hf(0, 0) * Hf(2, 1);
  cof(2, 2) = Hf(0, 0) * Hf(1, 1) - Hf(0, 1) * Hf(1, 0);

  const ScalarType num = g.transpose() * (cof * g);
  const ScalarType den = std::pow(norm_g, 4);

  return num / den;
}

template <class ScalarType, size_t MAX_REFINE_LEVEL>
ScalarType GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL>::curvedness(
    const ScalarType& x, const ScalarType& y, const ScalarType& z) const {
  const ScalarType H = this->meanCurvature(x, y, z);
  const ScalarType K = this->gaussianCurvature(x, y, z);
  if (!std::isfinite(H) || !std::isfinite(K))
    return std::numeric_limits<ScalarType>::quiet_NaN();

  const ScalarType C2 = 2.0 * H * H - K;
  return std::sqrt(std::max(ScalarType(0), C2));
}

template <class ScalarType, size_t MAX_REFINE_LEVEL>
std::pair<ScalarType, ScalarType>
GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL>::principalCurvatures(
    const ScalarType& x, const ScalarType& y, const ScalarType& z) const {
  const ScalarType H = this->meanCurvature(x, y, z);
  const ScalarType K = this->gaussianCurvature(x, y, z);

  if (!std::isfinite(H) || !std::isfinite(K)) {
    const ScalarType nan = std::numeric_limits<ScalarType>::quiet_NaN();
    return {nan, nan};
  }

  ScalarType disc = H * H - K;
  if (disc < ScalarType(0)) disc = ScalarType(0);

  const ScalarType r = std::sqrt(disc);
  ScalarType k1 = H + r;
  ScalarType k2 = H - r;

  if (k1 < k2) std::swap(k1, k2);
  return {k1, k2};
}

template <class ScalarType, size_t MAX_REFINE_LEVEL>
typename GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL>::Vec3
GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL>::projectPointOnSurface(
    const Vec3& p, int max_iter, ScalarType tol) const {
  Vec3 x_proj = p;
  for (int i = 0; i < max_iter; ++i) {
    const ScalarType f = this->F(x_proj[0], x_proj[1], x_proj[2]);
    const Vec3 g = this->gradF(x_proj[0], x_proj[1], x_proj[2]);
    const ScalarType g_norm2 = g.squaredNorm();
    if (g_norm2 < ScalarType(1e-14)) break;
    x_proj -= (f / g_norm2) * g;
    if (std::abs(f) < tol) break;
    if (i == max_iter - 1) {
      std::cout << "Max iterations reached. Projection incomplete. "
                << "f = " << std::abs(f) << std::endl;
    }
  }
  return x_proj;
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_IMPLICIT_SURFACES_GENERAL_IMPLICIT_SURFACE_TPP_