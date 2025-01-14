// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_CYLINDER_RECONSTRUCTIONS_ALIGNED_PARABOLOID_H_
#define IRL_CYLINDER_RECONSTRUCTIONS_ALIGNED_PARABOLOID_H_

#include <math.h>
#include <cassert>
#include <ostream>

#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/polynomial.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/quadratic_reconstruction/ellipse.h"

namespace IRL {

// Infinit cylinder in the form z^2 + b*y^2 = r
template <class ScalarType>
class AlignedCylinderBase {
 public:
  // using PolynomialBase<2, ScalarType>::PolynomialBase;
  using value_type = ScalarType;

  constexpr AlignedCylinderBase(void) {
    std::fill(coefficients_m.begin(), coefficients_m.end(), ScalarType(0));
  }

  constexpr AlignedCylinderBase(
      const std::array<ScalarType, 2>& a_coefficients) {
    coefficients_m = a_coefficients;
  }

  template <class OtherScalarType>
  constexpr AlignedCylinderBase(
      const AlignedCylinderBase<OtherScalarType>& a_aligned_cylinder) {
    coefficients_m[0] = ScalarType(a_aligned_cylinder.b());
    coefficients_m[1] = ScalarType(a_aligned_cylinder.r());
  };

  ScalarType& b(void) { return coefficients_m[0]; }
  ScalarType b(void) const { return coefficients_m[0]; }
  ScalarType& r(void) { return coefficients_m[1]; }
  ScalarType r(void) const { return coefficients_m[1]; }

  // Ellipse in the form Ax^2 + By^2 + Cxy + Dx + Ey + F = 0
  EllipseBase<ScalarType> intersectWithPlane(
      const PlaneBase<ScalarType>& a_plane) const {
    EllipseBase<ScalarType> ellipse_to_return;
    if constexpr (std::is_same<ScalarType, Quad_t>::value) {
      if (fabsq(a_plane.normal()[2]) < FLT128_EPSILON) {
        // no intersection
        return ellipse_to_return;
      }

    } else {
      if (std::fabs(a_plane.normal()[2]) < DBL_EPSILON) {
        // no intersection
        return ellipse_to_return;
      }
    }
    // IRL Plane defined as n_x x + n_y y + n_z z - d = 0
    auto normal = a_plane.normal();
    auto nx = normal[0];
    auto ny = normal[1];
    auto nz = normal[2];
    auto d = a_plane.distance();

    ellipse_to_return.a() = (nx * nx) / (nz * nz);
    ellipse_to_return.b() = this->b() + (ny * ny) / (nz * nz);
    ellipse_to_return.c() =
        static_cast<ScalarType>(2) * (nx * ny) / (nz * nz);  // AlignedParaboloid is aligned along x-y
    ellipse_to_return.d() = static_cast<ScalarType>(2) * d * nx / nz;
    ellipse_to_return.e() = static_cast<ScalarType>(2) * d * ny / nz;
    ellipse_to_return.f() = d * d - this->r();
    return ellipse_to_return;
  }

  void serialize(ByteBuffer* a_buffer) const {
    a_buffer->pack(coefficients_m.data(), 2);
  }
  void unpackSerialized(ByteBuffer* a_buffer) {
    a_buffer->unpack(coefficients_m.data(), 2);
  }

 private:
  std::array<ScalarType, 2> coefficients_m;
};

using AlignedCylinder = AlignedCylinderBase<double>;

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out,
    const AlignedCylinderBase<ScalarType>& a_aligned_cylinder) {
  out << a_aligned_cylinder.r() << " = z^2 + " << a_aligned_cylinder.b() << " * y^2";
  return out;
}

}  // namespace IRL

#endif  // IRL_CYLINDER_RECONSTRUCTIONS_ALIGNED_PARABOLOID_H_
