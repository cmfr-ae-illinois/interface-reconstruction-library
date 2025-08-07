// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_IMPLICIT_SURFACES_GENUS_TPP_
#define IRL_GEOMETRY_IMPLICIT_SURFACES_GENUS_TPP_

#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"

namespace IRL {

template <class ScalarType, size_t MAX_REFINE_LEVEL>
class Genus : public GeneralImplicitSurface<ScalarType, MAX_REFINE_LEVEL> {
 public:
  using Vec3 = Eigen::Matrix<ScalarType, 3, 1>;
  using Mat3 = Eigen::Matrix<ScalarType, 3, 3>;

  // Default constructor
  Genus() = default;

  ScalarType F(const ScalarType& x, const ScalarType& y,
               const ScalarType& z) const override {
    return ScalarType(2) * y * (y * y - ScalarType(3) * x * x) *
               (ScalarType(1) - z * z) +
           (x * x + y * y) * (x * x + y * y) -
           (ScalarType(9) * z * z - ScalarType(1)) * (ScalarType(1) - z * z);
  }

  Vec3 gradF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const override {
    ScalarType Fx = ScalarType(4) * x * (x * x + y * y) -
                    ScalarType(12) * x * y * (ScalarType(1) - z * z);

    ScalarType Fy = ScalarType(4) * y * (x * x + y * y) +
                    ScalarType(4) * y * y * (ScalarType(1) - z * z) +
                    ScalarType(2) * (-ScalarType(3) * x * x + y * y) *
                        (ScalarType(1) - z * z);

    ScalarType Fz =
        -ScalarType(4) * y * (-ScalarType(3) * x * x + y * y) * z -
        ScalarType(18) * z * (ScalarType(1) - z * z) +
        ScalarType(2) * z * (-ScalarType(1) + ScalarType(9) * z * z);

    return Vec3(Fx, Fy, Fz);
  }

  Mat3 hessF(const ScalarType& x, const ScalarType& y,
             const ScalarType& z) const override {
    Mat3 H = Mat3::Zero();

    H(0, 0) = ScalarType(8) * x * x + ScalarType(4) * (x * x + y * y) -
              ScalarType(12) * y * (ScalarType(1) - z * z);

    H(1, 1) = ScalarType(8) * y * y + ScalarType(4) * (x * x + y * y) +
              ScalarType(12) * y * (ScalarType(1) - z * z);

    H(2, 2) = -ScalarType(4) * y * (-ScalarType(3) * x * x + y * y) +
              ScalarType(72) * z * z -
              ScalarType(18) * (ScalarType(1) - z * z) +
              ScalarType(2) * (-ScalarType(1) + ScalarType(9) * z * z);

    H(0, 1) =
        ScalarType(8) * x * y - ScalarType(12) * x * (ScalarType(1) - z * z);
    H(1, 0) = H(0, 1);

    H(0, 2) = ScalarType(24) * x * y * z;
    H(2, 0) = H(0, 2);

    H(1, 2) = -ScalarType(8) * y * y * z -
              ScalarType(4) * (-ScalarType(3) * x * x + y * y) * z;
    H(2, 1) = H(1, 2);

    return H;
  }
};

}  // namespace IRL

#endif  // IRL_GEOMETRY_IMPLICIT_SURFACES_GENUS_TPP_
