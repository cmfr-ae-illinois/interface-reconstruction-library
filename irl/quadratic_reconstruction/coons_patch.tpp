// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_QUADRATIC_RECONSTRUCTION_COONS_PATCH_TPP_
#define IRL_QUADRATIC_RECONSTRUCTION_COONS_PATCH_TPP_

#include "irl/quadratic_reconstruction/coons_patch.h"

namespace IRL {

template <class ScalarType>
inline CoonsPatchBase<ScalarType>::CoonsPatchBase(void)
    : c0(), c1(), c2(), c3() {}

template <class ScalarType>
inline CoonsPatchBase<ScalarType>::CoonsPatchBase(
    const RationalBezierArcBase<ScalarType>& arc_0,
    const RationalBezierArcBase<ScalarType>& arc_1,
    const RationalBezierArcBase<ScalarType>& arc_2,
    const RationalBezierArcBase<ScalarType>& arc_3)
    : c0(arc_0), c1(arc_1), c2(arc_2), c3(arc_3) {}

// function definitions

template <class ScalarType>
inline std::array<RationalBezierArcBase<ScalarType>, 4>
CoonsPatchBase<ScalarType>::getBoundaryArcs() const {
  return {{c0, c1, c2, c3}};
}

template <class ScalarType>
inline RationalBezierArcBase<ScalarType> CoonsPatchBase<ScalarType>::getArc(
    const int edgeIndex) const {
  switch (edgeIndex) {
    case 0:
      return c0;
    case 1:
      return c1;
    case 2:
      return c2;
    case 3:
      return c3;
    default:
      throw std::out_of_range(
          "CoonsPatchBase::getArc: edgeIndex must be in [0,3]");
  }
}

template <class ScalarType>
inline PtBase<ScalarType> CoonsPatchBase<ScalarType>::evaluate(
    const ScalarType u, const ScalarType v) const {
  using Pt = PtBase<ScalarType>;

  Pt X = (ScalarType(1) - v) * c0.point(u) + v * c2.point(u) +
         (ScalarType(1) - u) * c3.point(v) + u * c1.point(v) -
         (ScalarType(1) - u) * (ScalarType(1) - v) * c0.start_point() -
         (ScalarType(1) - u) * v * c2.start_point() -
         u * (ScalarType(1) - v) * c1.start_point() - u * v * c2.end_point();

  return X;
}

template <class ScalarType>
inline Eigen::Matrix<ScalarType, 2, 2> CoonsPatchBase<ScalarType>::jacobian(
    const ScalarType u, const ScalarType v) const {
  // boundary‐curve derivatives
  auto dc0_du = c0.derivative(u);
  auto dc2_du = c2.derivative(u);
  auto dc1_dv = c1.derivative(v);
  auto dc3_dv = c3.derivative(v);

  // corner‐points of the bilinear blend
  auto P00 = c0.start_point();
  auto P01 = c2.start_point();
  auto P10 = c1.start_point();
  auto P11 = c2.end_point();

  // points on the curves
  auto c0_u = c0.point(u);
  auto c2_u = c2.point(u);
  auto c3_v = c3.point(v);
  auto c1_v = c1.point(v);

  // partials
  ScalarType dxdu = (ScalarType(1) - v) * dc0_du[0] + v * dc2_du[0] - c3_v[0] +
                    c1_v[0] + (ScalarType(1) - v) * P00[0] + v * P01[0] -
                    (ScalarType(1) - v) * P10[0] - v * P11[0];

  ScalarType dydu = (ScalarType(1) - v) * dc0_du[1] + v * dc2_du[1] - c3_v[1] +
                    c1_v[1] + (ScalarType(1) - v) * P00[1] + v * P01[1] -
                    (ScalarType(1) - v) * P10[1] - v * P11[1];

  ScalarType dxdv = -c0_u[0] + c2_u[0] + (ScalarType(1) - u) * dc3_dv[0] +
                    u * dc1_dv[0] + (ScalarType(1) - u) * P00[0] -
                    (ScalarType(1) - u) * P01[0] + u * P10[0] - u * P11[0];

  ScalarType dydv = -c0_u[1] + c2_u[1] + (ScalarType(1) - u) * dc3_dv[1] +
                    u * dc1_dv[1] + (ScalarType(1) - u) * P00[1] -
                    (ScalarType(1) - u) * P01[1] + u * P10[1] - u * P11[1];

  Eigen::Matrix<ScalarType, 2, 2> J;
  J(0, 0) = dxdu;
  J(0, 1) = dxdv;
  J(1, 0) = dydu;
  J(1, 1) = dydv;
  return J;
}

// — detJacobian: just take the determinant of the above
template <class ScalarType>
inline ScalarType CoonsPatchBase<ScalarType>::detJacobian(
    const ScalarType u, const ScalarType v) const {
  auto J = this->jacobian(u, v);
  return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
}

}  // namespace IRL

#endif  // IRL_QUADRATIC_RECONSTRUCTION_COONS_PATCH_TPP_
