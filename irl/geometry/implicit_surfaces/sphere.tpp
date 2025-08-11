// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_IMPLICIT_SURFACES_SPHERE_TPP_
#define IRL_GEOMETRY_IMPLICIT_SURFACES_SPHERE_TPP_

#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"

namespace IRL {

template <class ScalarType, size_t MAX_REFINE_LEVEL>
class Sphere : public GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL> {
 public:
  using Vec3 = Eigen::Matrix<ScalarType, 3, 1>;
  using Mat3 = Eigen::Matrix<ScalarType, 3, 3>;

  // constructors
  Sphere()
      : xc(ScalarType(0.)),
        yc(ScalarType(0.)),
        zc(ScalarType(0.)),
        R(ScalarType(0.15)) {}
  Sphere(const ScalarType& x_center, const ScalarType& y_center,
         const ScalarType& z_center, const ScalarType& radius)
      : xc(x_center), yc(y_center), zc(z_center), R(radius) {}

  ScalarType F(const ScalarType& x, const ScalarType& y,
               const ScalarType& z) const {
    return (x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc) -
           R * R;
  }

  Vec3 gradF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const {
    ScalarType Fx = ScalarType(2.0) * (x - xc);
    ScalarType Fy = ScalarType(2.0) * (y - yc);
    ScalarType Fz = ScalarType(2.0) * (z - zc);
    return Vec3(Fx, Fy, Fz);
  }

  Mat3 hessF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const override {
    return Mat3::Identity() * ScalarType(2.);
  }

  bool hasVolume() const override { return true; }

  ScalarType volume() const override { return (4.0 / 3.0) * M_PI * R * R * R; }

  bool hasSurfaceArea() const override { return true; }

  ScalarType surfaceArea() const override { return 4.0 * M_PI * R * R; }

 private:
  ScalarType xc, yc, zc, R;
};

}  // namespace IRL

#endif  // IRL_GEOMETRY_IMPLICIT_SURFACES_SPHERE_TPP_
