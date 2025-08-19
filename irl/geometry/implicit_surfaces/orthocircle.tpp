// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_IMPLICIT_SURFACES_ORTHOCIRCLE_TPP_
#define IRL_GEOMETRY_IMPLICIT_SURFACES_ORTHOCIRCLE_TPP_

#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"

namespace IRL {

template <class ScalarType, size_t MAX_REFINE_LEVEL>
class Orthocircle
    : public GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL> {
 public:
  using Vec3 = Eigen::Matrix<ScalarType, 3, 1>;
  using Mat3 = Eigen::Matrix<ScalarType, 3, 3>;

  // Default constructor
  Orthocircle() : c1(ScalarType(0.075)), c2(ScalarType(3)) {}

  ScalarType F(const ScalarType& x, const ScalarType& y,
               const ScalarType& z) const override {
    ScalarType x2 = x * x, y2 = y * y, z2 = z * z;
    ScalarType A = x2 + y2 - ScalarType(1);
    ScalarType B = y2 + z2 - ScalarType(1);
    ScalarType C = z2 + x2 - ScalarType(1);
    return (A * A + z2) * (B * B + x2) * (C * C + y2) -
           c1 * c1 * (ScalarType(1) + c2 * (x2 + y2 + z2));
  }

  Vec3 gradF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const override {
    ScalarType x2 = x * x, y2 = y * y, z2 = z * z;
    ScalarType A = x2 + y2 - ScalarType(1);
    ScalarType B = y2 + z2 - ScalarType(1);
    ScalarType C = z2 + x2 - ScalarType(1);
    ScalarType c1sq = c1 * c1;

    ScalarType Fx = -ScalarType(2) * c1sq * c2 * x +
                    ScalarType(2) * x * (A * A + z2) * (C * C + y2) +
                    ScalarType(4) * x * C * (A * A + z2) * (B * B + x2) +
                    ScalarType(4) * x * A * (C * C + y2) * (B * B + x2);

    ScalarType Fy = -ScalarType(2) * c1sq * c2 * y +
                    ScalarType(4) * y * B * (A * A + z2) * (C * C + y2) +
                    ScalarType(2) * y * (A * A + z2) * (B * B + x2) +
                    ScalarType(4) * y * A * (B * B + x2) * (C * C + y2);

    ScalarType Fz = -ScalarType(2) * c1sq * c2 * z +
                    ScalarType(4) * z * B * (A * A + z2) * (C * C + y2) +
                    ScalarType(4) * z * C * (A * A + z2) * (B * B + x2) +
                    ScalarType(2) * z * (B * B + x2) * (C * C + y2);

    return Vec3(Fx, Fy, Fz);
  }

  Mat3 hessF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const override {
    Mat3 H = Mat3::Zero();

    // Precompute powers and terms
    ScalarType x2 = x * x, y2 = y * y, z2 = z * z;
    ScalarType A = x2 + y2 - ScalarType(1);
    ScalarType B = y2 + z2 - ScalarType(1);
    ScalarType C = z2 + x2 - ScalarType(1);
    ScalarType c1sq = c1 * c1;

    H(0, 0) = -ScalarType(2) * c1sq * c2 +
              ScalarType(2) * (A * A + z2) * (y2 + C * C) +
              ScalarType(4) * x *
                  (ScalarType(4) * x * C * (A * A + z2) +
                   ScalarType(4) * x * (A) * (y2 + C * C)) +
              (x2 + B * B) *
                  (ScalarType(32) * x2 * (A) * (C) +
                   (A * A + z2) * (ScalarType(8) * x2 + ScalarType(4) * C) +
                   (ScalarType(8) * x2 + ScalarType(4) * A) * (y2 + C * C));

    H(1, 1) =
        -ScalarType(2) * c1sq * c2 +
        (A * A + z2) * (y2 + C * C) * (ScalarType(8) * y2 + ScalarType(4) * B) +
        ScalarType(8) * y * B *
            (ScalarType(2) * y * (A * A + z2) +
             ScalarType(4) * y * A * (y2 + C * C)) +
        (x2 + B * B) *
            (ScalarType(16) * y2 * (A) + ScalarType(2) * (A * A + z2) +
             (ScalarType(8) * y2 + ScalarType(4) * A) * (y2 + C * C));

    H(2, 2) =
        -ScalarType(2) * c1sq * c2 +
        (A * A + z2) * (y2 + C * C) * (ScalarType(8) * z2 + ScalarType(4) * B) +
        (x2 + B * B) *
            (ScalarType(16) * z2 * (C) +
             (A * A + z2) * (ScalarType(8) * z2 + ScalarType(4) * C) +
             ScalarType(2) * (y2 + C * C)) +
        ScalarType(8) * z * B *
            (ScalarType(4) * z * C * (A * A + z2) +
             ScalarType(2) * z * (y2 + C * C));

    H(0, 1) = ScalarType(4) * x * y * (A * A + z2) +
              ScalarType(16) * x * y * C * B * (A * A + z2) +
              ScalarType(8) * x * y * A * (y2 + C * C) +
              ScalarType(16) * x * y * A * B * (y2 + C * C) +
              ScalarType(8) * x * y * A * (x2 + B * B) +
              ScalarType(16) * x * y * A * C * (x2 + B * B) +
              ScalarType(8) * x * y * (y2 + C * C) * (x2 + B * B);
    H(1, 0) = H(0, 1);

    H(0, 2) = ScalarType(8) * x * z * C * (A * A + z2) +
              ScalarType(16) * x * z * C * B * (A * A + z2) +
              ScalarType(4) * x * z * (y2 + C * C) +
              ScalarType(16) * x * A * z * B * (y2 + C * C) +
              ScalarType(8) * x * z * C * (x2 + B * B) +
              ScalarType(16) * x * A * z * C * (x2 + B * B) +
              ScalarType(8) * x * z * (A * A + z2) * (x2 + B * B);
    H(2, 0) = H(0, 2);

    H(1, 2) = ScalarType(8) * y * z * B * (A * A + z2) +
              ScalarType(16) * y * z * C * B * (A * A + z2) +
              ScalarType(8) * y * z * B * (y2 + C * C) +
              ScalarType(16) * y * A * z * B * (y2 + C * C) +
              ScalarType(8) * y * z * (A * A + z2) * (y2 + C * C) +
              ScalarType(4) * y * z * (x2 + B * B) +
              ScalarType(16) * y * A * z * C * (x2 + B * B);
    H(2, 1) = H(1, 2);

    return H;
  }

 private:
  ScalarType c1, c2;
};

}  // namespace IRL

#endif  // IRL_GEOMETRY_IMPLICIT_SURFACES_ORTHOCIRCLE_TPP_
