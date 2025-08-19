// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_IMPLICIT_SURFACES_ELLIPSOID_TPP_
#define IRL_GEOMETRY_IMPLICIT_SURFACES_ELLIPSOID_TPP_

#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"

namespace IRL {

template <class ScalarType, size_t MAX_REFINE_LEVEL>
class Ellipsoid : public GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL> {
 public:
  using Vec3 = Eigen::Matrix<ScalarType, 3, 1>;
  using Mat3 = Eigen::Matrix<ScalarType, 3, 3>;

  // constructors
  Ellipsoid()
      : xc(ScalarType(0)),
        yc(ScalarType(0)),
        zc(ScalarType(0)),
        a(ScalarType(0.3)),
        b(ScalarType(0.15)),
        c(ScalarType(0.1)) {}

  Ellipsoid(const ScalarType& x_center, const ScalarType& y_center,
            const ScalarType& z_center, const ScalarType& axis_a,
            const ScalarType& axis_b, const ScalarType& axis_c)
      : xc(x_center),
        yc(y_center),
        zc(z_center),
        a(axis_a),
        b(axis_b),
        c(axis_c) {}

  ScalarType F(const ScalarType& x, const ScalarType& y,
               const ScalarType& z) const override {
    return ((x - xc) / a) * ((x - xc) / a) + ((y - yc) / b) * ((y - yc) / b) +
           ((z - zc) / c) * ((z - zc) / c) - ScalarType(1.);
  }

  Vec3 gradF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const override {
    return Vec3(ScalarType(2.) * (x - xc) / (a * a),
                ScalarType(2.) * (y - yc) / (b * b),
                ScalarType(2.) * (z - zc) / (c * c));
  }

  Mat3 hessF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const override {
    Mat3 H = Mat3::Zero();
    H(0, 0) = ScalarType(2.) / (a * a);
    H(1, 1) = ScalarType(2.) / (b * b);
    H(2, 2) = ScalarType(2.) / (c * c);
    return H;
  }

  bool hasVolume() const override { return true; }
  ScalarType volume() const override {
    return (ScalarType(4) / ScalarType(3)) * M_PI * a * b * c;
  }

 private:
  ScalarType xc, yc, zc;  // center
  ScalarType a, b, c;     // semi-axes
};

}  // namespace IRL

#endif  // IRL_GEOMETRY_IMPLICIT_SURFACES_ELLIPSOID_TPP_