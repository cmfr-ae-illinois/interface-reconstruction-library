// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_

#include <float.h>
#include <cassert>
#include <cmath>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/scalar_with_gradient.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/mymath.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace IRL {

/* This returns the algebraic signed distance to a paraboloid
 * returns < 0 number if a_pt is below the paraboloid
 * returns > 0 number if a_pt is above the paraboloid */
template <class ScalarType>
inline ScalarType signedDistance(
    const PtBase<ScalarType>& a_pt,
    const AlignedParaboloidBase<ScalarType>& a_paraboloid) {
  return a_paraboloid.a() * a_pt[0] * a_pt[0] +
         a_paraboloid.b() * a_pt[1] * a_pt[1] + a_pt[2];
}

/******************************************************************************/
/*********************** First moment contribution ****************************/
/******************************************************************************/
/* This compute the first contribution to the moments (arising from the
 * integration of the face plane primitives on the poligonized clipped faces) */
template <class ReturnType, class ScalarType>
ReturnType computeType1Contribution(const PtBase<ScalarType>& a_ref_pt,
                                    const PtBase<ScalarType>& a_pt_0,
                                    const PtBase<ScalarType>& a_pt_1) {
  if constexpr (is_moments_volume<ReturnType>::value) {
    return ReturnType::fromScalarConstant(
        (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) / ScalarType(6.0) *
        ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
         (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1])));
  } else {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType ONETWELVTH = ONE / ScalarType(12);

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ZERO);
    const ScalarType triangle_area =
        ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
         (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1])) /
        TWO;
    moments.volume() =
        triangle_area * (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) / THREE;
    moments.centroid()[0] =
        triangle_area *
        (a_pt_0[0] * (TWO * a_pt_0[2] + a_pt_1[2] + a_ref_pt[2]) +
         a_pt_1[0] * (a_pt_0[2] + TWO * a_pt_1[2] + a_ref_pt[2]) +
         a_ref_pt[0] * (a_pt_0[2] + a_pt_1[2] + TWO * a_ref_pt[2])) *
        ONETWELVTH;
    moments.centroid()[1] =
        triangle_area *
        (a_pt_0[1] * (TWO * a_pt_0[2] + a_pt_1[2] + a_ref_pt[2]) +
         a_pt_1[1] * (a_pt_0[2] + TWO * a_pt_1[2] + a_ref_pt[2]) +
         a_ref_pt[1] * (a_pt_0[2] + a_pt_1[2] + TWO * a_ref_pt[2])) *
        ONETWELVTH;
    moments.centroid()[2] =
        triangle_area *
        (a_pt_0[2] * a_pt_0[2] + a_pt_1[2] * a_pt_1[2] +
         a_ref_pt[2] * a_ref_pt[2] + a_pt_1[2] * a_ref_pt[2] +
         a_pt_0[2] * a_pt_1[2] + a_pt_0[2] * a_ref_pt[2]) *
        ONETWELVTH;
    return moments;
  }
}

template <class ReturnType, class ScalarType>
ReturnType computeType2Contribution(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1) {
  if constexpr (is_moments_volume<ReturnType>::value) {
    const ScalarType ONETWELVTH = ScalarType(1) / ScalarType(12);
    return ReturnType::fromScalarConstant(
        ONETWELVTH * (a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) *
        (-a_pt_0[2] - a_pt_1[2] +
         a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] +
         a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1]));
  } else {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType ONETWELVTH = ONE / ScalarType(12);
    const ScalarType ONE60TH = ONE / ScalarType(60);
    const ScalarType ONE180TH = ONE / ScalarType(180);

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ZERO);
    moments.volume() = (a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) *
                       ONETWELVTH *
                       (-a_pt_0[2] - a_pt_1[2] +
                        a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] +
                        a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1]);
    moments.centroid()[0] =
        (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
        (TWO * a_aligned_paraboloid.b() * (a_pt_0[1] - a_pt_1[1]) *
             (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) +
         THREE * (a_pt_0[0] + a_pt_1[0]) * (a_pt_0[2] + a_pt_1[2])) *
        ONE60TH;
    moments.centroid()[1] =
        (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
        (TWO * a_aligned_paraboloid.a() * (a_pt_0[0] - a_pt_1[0]) *
             (a_pt_1[1] * a_pt_0[0] - a_pt_0[1] * a_pt_1[0]) +
         THREE * (a_pt_0[1] + a_pt_1[1]) * (a_pt_0[2] + a_pt_1[2])) *
        ONE60TH;
    moments.centroid()[2] =
        ((a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) *
         (TWO * a_aligned_paraboloid.a() * a_aligned_paraboloid.b() *
              ((a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
               (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1])) +
          THREE * a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] *
              (a_pt_0[2] + a_pt_1[2]) +
          THREE * a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1] *
              (a_pt_0[2] + a_pt_1[2]) -
          THREE * (a_pt_0[2] * a_pt_0[2] + a_pt_0[2] * a_pt_1[2] +
                   a_pt_1[2] * a_pt_1[2]))) *
        ONE180TH;
    return moments;
  }
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsV3SeriesOne(
    const ContainerType& a_weight) {
  std::array<ContainerType, 3> coeffs;
  coeffs.fill(ContainerType(ScalarType(0)));
  ContainerType x(ScalarType(1));
  UnsignedIndex_t i = 0;
  ScalarType max_diff;
  while (i <= 40) {
    max_diff = ScalarType(0);
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      ContainerType add_to_coeff = ScalarType(v3Series[i][j]) * x;
      coeffs[j] += add_to_coeff;
      max_diff = maximum(max_diff, fabs(add_to_coeff));
    }
    if (max_diff < ScalarType(DBL_EPSILON)) {
      break;
    }
    x *= a_weight - ContainerType(ScalarType(1));
    i++;
  }
  return coeffs;
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 12> coeffsV3andC3SeriesOne(
    const ContainerType& a_weight) {
  std::array<ContainerType, 12> coeffs;
  coeffs.fill(ContainerType(ScalarType(0)));
  ContainerType x(ScalarType(1));
  UnsignedIndex_t i = 0;
  ScalarType max_diff;
  while (i <= 40) {
    max_diff = ScalarType(0);
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      ContainerType add_to_coeff = ScalarType(v3Series[i][j]) * x;
      coeffs[j] += add_to_coeff;
      max_diff = maximum(max_diff, fabs(add_to_coeff));
    }
    if (max_diff < ScalarType(DBL_EPSILON)) {
      break;
    }
    x *= a_weight - ContainerType(ScalarType(1));
    i++;
  }
  i = 0;
  x = ContainerType(ScalarType(1));
  while (i <= 40) {
    max_diff = ScalarType(0);
    for (UnsignedIndex_t j = 0; j < 4; ++j) {
      ContainerType add_to_coeff = ScalarType(cx3Series[i][j]) * x;
      coeffs[3 + j] += add_to_coeff;
      max_diff = maximum(max_diff, fabs(add_to_coeff));
    }
    for (UnsignedIndex_t j = 0; j < 5; ++j) {
      ContainerType add_to_coeff = ScalarType(cz3Series[i][j]) * x;
      coeffs[7 + j] += add_to_coeff;
      max_diff = maximum(max_diff, fabs(add_to_coeff));
    }
    if (max_diff < ScalarType(DBL_EPSILON)) {
      break;
    }
    x *= a_weight - ContainerType(ScalarType(1));
    i++;
  }

  return coeffs;
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsV3Exact(
    const ContainerType& a_weight) {
  /* Defining constants and types */
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType THREE = ScalarType(3);
  const ScalarType FOUR = ScalarType(4);
  const ScalarType SIX = ScalarType(6);

  /* Function */
  const auto w2 = a_weight * a_weight;
  const auto w3 = w2 * a_weight;
  const auto w4 = w2 * w2;
  const auto w5 = w2 * w3;
  const auto w6 = w3 * w3;
  const auto L = ContainerType(ONE) / ((a_weight - ContainerType(ONE)) *
                                       (a_weight + ContainerType(ONE)));
  const auto L3 = L * L * L;
  const auto S = (a_weight < ContainerType(ONE))
                     ? sqrt(ContainerType(ONE) - a_weight * a_weight)
                     : sqrt(a_weight * a_weight - ContainerType(ONE));
  const auto T = (a_weight < ContainerType(ONE))
                     ? atan((ContainerType(ONE) - a_weight) / S) / S
                     : atanh((a_weight - ContainerType(ONE)) / S) / S;
  return std::array<ContainerType, 3>(
      {(TWO * w6 - THREE * w4 + ScalarType(31) * w2 -
        (ScalarType(42) * w3 + ScalarType(18) * a_weight) * T) *
           L3 / ScalarType(12),
       (TWO * w6 - ScalarType(9) * w4 - ScalarType(8) * w2 +
        ScalarType(30) * w3 * T) *
           L3 / THREE,
       (ScalarType(11) * w4 + FOUR * w2 -
        (ScalarType(12) * w5 + ScalarType(18) * w3) * T) *
           L3 / SIX});
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 12> coeffsV3andC3Exact(
    const ContainerType& a_weight) {
  /* Defining constants and types */
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType THREE = ScalarType(3);
  const ScalarType FOUR = ScalarType(4);
  const ScalarType FIVE = ScalarType(5);
  const ScalarType SIX = ScalarType(6);

  /* Function */
  const auto w2 = a_weight * a_weight;
  const auto w3 = w2 * a_weight;
  const auto w4 = w2 * w2;
  const auto w5 = w2 * w3;
  const auto w6 = w3 * w3;
  const auto w7 = w4 * w3;
  const auto w8 = w4 * w4;
  const auto w9 = w4 * w5;
  const auto w10 = w5 * w5;
  const auto L = ContainerType(ONE) / ((a_weight - ContainerType(ONE)) *
                                       (a_weight + ContainerType(ONE)));
  const auto L3 = L * L * L;
  const auto L4 = L3 * L;
  const auto L5 = L4 * L;
  const auto S = (a_weight < ContainerType(ONE))
                     ? sqrt(ContainerType(ONE) - a_weight * a_weight)
                     : sqrt(a_weight * a_weight - ContainerType(ONE));
  const auto T = (a_weight < ContainerType(ONE))
                     ? atan((ContainerType(ONE) - a_weight) / S) / S
                     : atanh((a_weight - ContainerType(ONE)) / S) / S;
  return std::array<ContainerType, 12>(
      {L3 *
           (TWO * w6 - THREE * w4 + ScalarType(31) * w2 -
            (ScalarType(42) * w3 + ScalarType(18) * a_weight) * T) /
           ScalarType(12),
       L3 *
           (TWO * w6 - ScalarType(9) * w4 - ScalarType(8) * w2 +
            ScalarType(30) * w3 * T) /
           THREE,
       L3 *
           (ScalarType(11) * w4 + FOUR * w2 -
            (ScalarType(12) * w5 + ScalarType(18) * w3) * T) /
           SIX,
       L4 * ((-T * a_weight) / ScalarType(32) +
             (ScalarType(93) * (w2)) / ScalarType(2240) -
             (ScalarType(163) * (w4)) / ScalarType(3360) +
             (FIVE * (w6)) / ScalarType(168) - (w8) / ScalarType(140)),
       L4 * ((w2) / ScalarType(70) + (-T * (w3)) / ScalarType(16) +
             (ScalarType(29) * (w4)) / ScalarType(1120) -
             (ScalarType(19) * (w6)) / ScalarType(1680) +
             (w8) / ScalarType(420)),
       -L4 *
           ((w2) / ScalarType(210) - (w4) / ScalarType(21) -
            (-T * (w5)) / ScalarType(8) -
            (ScalarType(13) * (w6)) / ScalarType(560) + (w8) / ScalarType(280)),
       L4 * ((w2) / ScalarType(35) - (ScalarType(16) * (w4)) / ScalarType(105) +
             (ScalarType(58) * (w6)) / ScalarType(105) - T * (w7) +
             (w8) / ScalarType(14)),
       L5 * ((-T * a_weight) / ScalarType(128) +
             (ScalarType(193) * (w2)) / ScalarType(16128) -
             (ScalarType(149) * (w4)) / ScalarType(8064) +
             (ScalarType(19) * (w6)) / ScalarType(1120) -
             (ScalarType(41) * (w8)) / ScalarType(5040) +
             (w10) / ScalarType(630)),
       L5 * ((FOUR * (w2)) / ScalarType(945) + (-T * (w3)) / ScalarType(48) +
             (ScalarType(65) * (w4)) / ScalarType(6048) -
             (w6) / ScalarType(144) +
             (ScalarType(11) * (w8)) / ScalarType(3780) -
             (w10) / ScalarType(1890)),
       -L5 * ((w2) / ScalarType(1890) -
              (ScalarType(13) * (w4)) / ScalarType(1890) -
              (-T * (w5)) / ScalarType(48) -
              (ScalarType(11) * (w6)) / ScalarType(2016) +
              (ScalarType(5) * (w8)) / ScalarType(3024) -
              (w10) / ScalarType(3780)),
       L5 * ((w2) / ScalarType(315) - (w4) / ScalarType(45) +
             (ScalarType(4) * (w6)) / ScalarType(35) +
             (-T * (w7)) / ScalarType(4) +
             (ScalarType(17) * (w8)) / ScalarType(504) -
             (w10) / ScalarType(252)),
       -L5 *
           ((w2) / ScalarType(63) - (ScalarType(29) * (w4)) / ScalarType(315) +
            (ScalarType(26) * (w6)) / ScalarType(105) -
            (ScalarType(194) * (w8)) / ScalarType(315) + T * (w9) -
            (w10) / ScalarType(18))});
}

template <class ReturnType, class ScalarType>
ReturnType computeType3Contribution(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const RationalBezierArcBase<ScalarType>& a_arc) {
  if constexpr (is_moments_volume<ReturnType>::value) {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType FOUR = ScalarType(4);
    const ScalarType FIVE = ScalarType(5);
    const ScalarType SIX = ScalarType(6);
    const ScalarType HALF = ONE / TWO;
    const ScalarType QUARTER = ONE / FOUR;

    const auto& pt_0 = a_arc.start_point();
    const auto& cp = a_arc.control_point();
    const auto& pt_1 = a_arc.end_point();
    const auto& weight = a_arc.weight();
    const ScalarType area_proj_triangle =
        HALF * (pt_0[0] * (pt_1[1] - cp[1]) + pt_1[0] * (cp[1] - pt_0[1]) +
                cp[0] * (pt_0[1] - pt_1[1]));
    assert(weight >= ZERO);
    std::array<ScalarType, 3> coeffs;
    if (weight < ScalarType(0.35))  // We use the exact expressions
      coeffs = coeffsV3Exact<ScalarType, ScalarType>(weight);
    else if (weight <
             ScalarType(1.7))  // We use the 40th order Taylor series (w -> 1)
      coeffs = coeffsV3SeriesOne<ScalarType, ScalarType>(weight);
    else if (weight <
             ScalarType(1.0e9))  // We use the series expansion (w -> infty)
      coeffs = coeffsV3Exact<ScalarType, ScalarType>(weight);
    else  // This is within EPSILON of the actual value
      coeffs = std::array<ScalarType, 3>({ONE / SIX, TWO / THREE, ZERO});

    return ReturnType::fromScalarConstant(
        area_proj_triangle *
        (coeffs[0] *
             signedDistance<ScalarType>(HALF * (pt_0 + pt_1), a_paraboloid) +
         coeffs[1] * signedDistance<ScalarType>(
                         QUARTER * (pt_0 + pt_1) + HALF * cp, a_paraboloid) +
         coeffs[2] * signedDistance<ScalarType>(cp, a_paraboloid)));
  } else {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType FOUR = ScalarType(4);
    const ScalarType FIVE = ScalarType(5);
    const ScalarType SIX = ScalarType(6);
    const ScalarType HALF = ONE / TWO;
    const ScalarType QUARTER = ONE / FOUR;

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ZERO);
    const auto& pt_0 = a_arc.start_point();
    const auto& cp = a_arc.control_point();
    const auto& pt_1 = a_arc.end_point();
    const auto weight = a_arc.weight();
    const auto A = a_paraboloid.a(), B = a_paraboloid.b();
    const auto X0 = pt_0[0], X1 = cp[0], X2 = pt_1[0];
    const auto Y0 = pt_0[1], Y1 = cp[1], Y2 = pt_1[1];
    const auto Z0 = pt_0[2], Z1p = cp[2], Z2 = pt_1[2];
    const ScalarType Z1P = -A * X1 * X1 - B * Y1 * Y1;
    const ScalarType AA = A * A, BB = B * B, AB = A * B;
    const ScalarType X00 = X0 * X0, X11 = X1 * X1, X22 = X2 * X2;
    const ScalarType X000 = X00 * X0, X111 = X11 * X1, X222 = X22 * X2;
    const ScalarType X0000 = X00 * X00, X1111 = X11 * X11, X2222 = X22 * X22;
    const ScalarType Y00 = Y0 * Y0, Y11 = Y1 * Y1, Y22 = Y2 * Y2;
    const ScalarType Y000 = Y00 * Y0, Y111 = Y11 * Y1, Y222 = Y22 * Y2;
    const ScalarType Y0000 = Y00 * Y00, Y1111 = Y11 * Y11, Y2222 = Y22 * Y22;
    const ScalarType Z00 = Z0 * Z0, Z22 = Z2 * Z2;
    const ScalarType Z1p1p = Z1p * Z1p, Z1P1P = Z1P * Z1P;
    const ScalarType X02 = X0 * X2, X12 = X1 * X2, X01 = X0 * X1;
    const ScalarType Y02 = Y0 * Y2, Y12 = Y1 * Y2, Y01 = Y0 * Y1;
    const ScalarType Z02 = Z0 * Z2, Z01p = Z0 * Z1p, Z01P = Z0 * Z1P;
    const ScalarType Z1p1P = Z1p * Z1P, Z1p2 = Z1p * Z2, Z1P2 = Z1P * Z2;
    const ScalarType X0Z0 = X0 * Z0, X0Z1p = X0 * Z1p, X0Z1P = X0 * Z1P,
                     X0Z2 = X0 * Z2;
    const ScalarType X1Z0 = X1 * Z0, X1Z1p = X1 * Z1p, X1Z1P = X1 * Z1P,
                     X1Z2 = X1 * Z2;
    const ScalarType X2Z0 = X2 * Z0, X2Z1p = X2 * Z1p, X2Z1P = X2 * Z1P,
                     X2Z2 = X2 * Z2;
    const ScalarType Y0Z0 = Y0 * Z0, Y0Z1p = Y0 * Z1p, Y0Z1P = Y0 * Z1P,
                     Y0Z2 = Y0 * Z2;
    const ScalarType Y1Z0 = Y1 * Z0, Y1Z1p = Y1 * Z1p, Y1Z1P = Y1 * Z1P,
                     Y1Z2 = Y1 * Z2;
    const ScalarType Y2Z0 = Y2 * Z0, Y2Z1p = Y2 * Z1p, Y2Z1P = Y2 * Z1P,
                     Y2Z2 = Y2 * Z2;
    const ScalarType area_proj_triangle =
        HALF * (X0 * (Y2 - Y1) + X2 * (Y1 - Y0) + X1 * (Y0 - Y2));
    assert(weight >= ZERO);
    // Compute coefficients (functions of the weight)
    std::array<ScalarType, 12> coeffs;
    if (weight < ScalarType(0.35))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
    } else if (weight <
               ScalarType(1.7))  // We use the 40th order Taylor series (w -> 1)
    {
      coeffs = coeffsV3andC3SeriesOne<ScalarType, ScalarType>(weight);
    } else if (weight < ScalarType(1.0e9))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
    }
    // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
    //   coeffs = coeffsV3andC3SeriesInfinity(weight);
    else  // This is within EPSILON of the actual value
    {
      coeffs = std::array<ScalarType, 12>(
          {ScalarType(ONE / SIX), ScalarType(TWO / THREE), ScalarType(ZERO),
           ScalarType(-ONE / ScalarType(140)),
           ScalarType(ONE / ScalarType(420)),
           ScalarType(-ONE / ScalarType(280)), ScalarType(ONE / ScalarType(14)),
           ScalarType(ONE / ScalarType(630)),
           ScalarType(-ONE / ScalarType(1890)),
           ScalarType(ONE / ScalarType(3780)),
           ScalarType(-ONE / ScalarType(252)),
           ScalarType(ONE / ScalarType(18))});
    }
    auto m0_basis = std::array<ScalarType, 3>(
        {signedDistance<ScalarType>(HALF * (pt_0 + pt_1), a_paraboloid),
         signedDistance<ScalarType>(QUARTER * (pt_0 + pt_1) + HALF * cp,
                                    a_paraboloid),
         signedDistance<ScalarType>(cp, a_paraboloid)});
    auto m1x_basis = std::array<ScalarType, 4>(
        {-SIX * (X0Z0 - X0Z2 - X2Z0 + X2Z2 - TWO * B * X2 * Y00 +
                 TWO * B * X0 * Y02 + TWO * B * X2 * Y02 - TWO * B * X0 * Y22),
         TWO * (FIVE * X0Z0 + ScalarType(10) * X0Z1p + SIX * X0Z1P +
                ScalarType(7) * X0Z2 + ScalarType(30) * A * X02 * X1 -
                ScalarType(11) * X1Z0 - FOUR * X1Z1p - ScalarType(11) * X1Z2 +
                ScalarType(7) * X2Z0 + ScalarType(10) * X2Z1p + SIX * X2Z1P +
                FIVE * X2Z2 - ScalarType(14) * B * X1 * Y00 +
                FOUR * B * X2 * Y00 + ScalarType(14) * B * X0 * Y01 -
                FOUR * B * X1 * Y01 + ScalarType(10) * B * X2 * Y01 -
                FOUR * B * X0 * Y02 + ScalarType(10) * B * X1 * Y02 -
                FOUR * B * X2 * Y02 + FOUR * B * X0 * Y11 +
                FOUR * B * X2 * Y11 + ScalarType(10) * B * X0 * Y12 -
                FOUR * B * X1 * Y12 + ScalarType(14) * B * X2 * Y12 +
                FOUR * B * X0 * Y22 - ScalarType(14) * B * X1 * Y22),
         TWO * (-FIVE * X0Z1p + ScalarType(18) * X0Z1P + X0Z2 +
                SIX * A * X02 * X1 - FIVE * X1Z0 - SIX * X1Z1p - SIX * X1Z1P -
                FIVE * X1Z2 + X2Z0 - FIVE * X2Z1p + ScalarType(18) * X2Z1P -
                ScalarType(12) * B * X1 * Y01 + TWO * B * X2 * Y01 +
                TWO * B * X1 * Y02 + ScalarType(12) * B * X0 * Y11 +
                ScalarType(12) * B * X2 * Y11 + TWO * B * X0 * Y12 -
                ScalarType(12) * B * X1 * Y12),
         TWO * (X1Z1p - X1Z1P)});
    auto m1y_basis = std::array<ScalarType, 4>(
        {SIX *
             (-Y0Z0 + Y0Z2 + TWO * A * (X22 * Y0 + X00 * Y2 - X02 * (Y0 + Y2)) +
              Y2Z0 - Y2Z2),
         TWO *
             (FIVE * Y0Z0 + ScalarType(10) * Y0Z1p + SIX * Y0Z1P +
              ScalarType(7) * Y0Z2 + ScalarType(30) * B * Y02 * Y1 -
              ScalarType(11) * Y1Z0 - FOUR * Y1Z1p - ScalarType(11) * Y1Z2 +
              TWO * A *
                  (-TWO * X02 * Y0 + TWO * X11 * Y0 + FIVE * X12 * Y0 +
                   TWO * X22 * Y0 - ScalarType(7) * X00 * Y1 + FIVE * X02 * Y1 -
                   TWO * X12 * Y1 - ScalarType(7) * X22 * Y1 + TWO * X00 * Y2 -
                   TWO * X02 * Y2 + TWO * X11 * Y2 + ScalarType(7) * X12 * Y2 +
                   X01 * (ScalarType(7) * Y0 - TWO * Y1 + FIVE * Y2)) +
              ScalarType(7) * Y2Z0 + ScalarType(10) * Y2Z1p + SIX * Y2Z1P +
              FIVE * Y2Z2),
         -TWO * (FIVE * Y0Z1p - ScalarType(18) * Y0Z1P - Y0Z2 -
                 SIX * B * Y02 * Y1 + FIVE * Y1Z0 + SIX * Y1Z1p + SIX * Y1Z1P +
                 FIVE * Y1Z2 -
                 TWO * A *
                     (X12 * Y0 - SIX * X01 * Y1 + X02 * Y1 - SIX * X12 * Y1 +
                      X01 * Y2 + SIX * X11 * (Y0 + Y2)) -
                 Y2Z0 + FIVE * Y2Z1p - ScalarType(18) * Y2Z1P),
         TWO * (Y1Z1p - Y1Z1P)});
    auto m1z_basis = std::array<ScalarType, 5>(
        {-(AA * (ScalarType(21) * X0000 + ScalarType(28) * X000 * X2 +
                 ScalarType(30) * X00 * X22 + ScalarType(28) * X0 * X222 +
                 ScalarType(21) * X2222)) -
             ScalarType(21) * BB * Y0000 - ScalarType(28) * BB * Y000 * Y2 -
             ScalarType(30) * BB * Y00 * Y22 -
             TWO * AB *
                 (X00 * (ScalarType(21) * Y00 + ScalarType(14) * Y02 +
                         FIVE * Y22) +
                  TWO * X02 *
                      (ScalarType(7) * Y00 + ScalarType(10) * Y02 +
                       ScalarType(7) * Y22) +
                  X22 * (ScalarType(5) * Y00 + ScalarType(14) * Y02 +
                         ScalarType(21) * Y22)) -
             ScalarType(28) * BB * Y0 * Y222 - ScalarType(21) * BB * Y2222 +
             ScalarType(40) * Z00 + ScalarType(48) * Z02 + ScalarType(40) * Z22,
         THREE * AA *
                 (ScalarType(21) * X000 * X1 - ScalarType(7) * X00 * X11 -
                  ScalarType(10) * X02 * X11 + ScalarType(35) * X00 * X12 -
                  ScalarType(7) * X000 * X2 - ScalarType(10) * X00 * X22 +
                  ScalarType(35) * X01 * X22 - ScalarType(7) * X11 * X22 -
                  ScalarType(7) * X0 * X222 + ScalarType(21) * X1 * X222) -
             AB * (ScalarType(7) * X11 * Y00 - ScalarType(35) * X12 * Y00 +
                   ScalarType(10) * X22 * Y00 - ScalarType(63) * X00 * Y01 +
                   ScalarType(20) * X12 * Y01 - ScalarType(35) * X22 * Y01 +
                   ScalarType(21) * X00 * Y02 + ScalarType(10) * X11 * Y02 -
                   ScalarType(70) * X12 * Y02 + ScalarType(21) * X22 * Y02 +
                   ScalarType(7) * X00 * Y11 + ScalarType(7) * X22 * Y11 -
                   ScalarType(35) * X00 * Y12 + ScalarType(28) * X12 * Y12 -
                   ScalarType(63) * X22 * Y12 +
                   X01 * (-ScalarType(63) * Y00 + ScalarType(28) * Y01 -
                          ScalarType(70) * Y02 + ScalarType(20) * Y12 -
                          ScalarType(35) * Y22) +
                   ScalarType(10) * X00 * Y22 + ScalarType(7) * X11 * Y22 -
                   ScalarType(63) * X12 * Y22 +
                   X02 * (ScalarType(21) * Y00 - ScalarType(70) * Y01 +
                          ScalarType(40) * Y02 + ScalarType(10) * Y11 -
                          ScalarType(70) * Y12 + ScalarType(21) * Y22)) -
             THREE *
                 (BB *
                      (ScalarType(7) * Y00 * Y11 + ScalarType(10) * Y02 * Y11 -
                       ScalarType(35) * Y00 * Y12 +
                       ScalarType(7) * Y000 * (-THREE * Y1 + Y2) +
                       ScalarType(10) * Y00 * Y22 - ScalarType(35) * Y01 * Y22 +
                       ScalarType(7) * Y11 * Y22 + ScalarType(7) * Y0 * Y222 -
                       ScalarType(21) * Y1 * Y222) +
                  TWO * (ScalarType(5) * Z00 + ScalarType(10) * Z01p +
                         FOUR * Z02 - TWO * Z1p1p + ScalarType(10) * Z1p2 +
                         FIVE * Z22)),
         -SIX * AA *
                 (ScalarType(46) * X02 * X11 - ScalarType(14) * X0 * X111 +
                  X1111 - ScalarType(14) * X111 * X2 -
                  ScalarType(14) * X01 * X22 + ScalarType(28) * X11 * X22 +
                  X00 * (ScalarType(28) * X11 - ScalarType(14) * X12 + X22)) -
             TWO * AB *
                 (X22 * Y00 + ScalarType(112) * X01 * Y01 -
                  ScalarType(28) * X02 * Y01 - ScalarType(14) * X22 * Y01 -
                  ScalarType(28) * X01 * Y02 + FOUR * X02 * Y02 +
                  ScalarType(28) * X00 * Y11 - ScalarType(42) * X01 * Y11 +
                  ScalarType(46) * X02 * Y11 + ScalarType(28) * X22 * Y11 -
                  TWO * X12 *
                      (ScalarType(7) * Y00 - ScalarType(46) * Y01 +
                       ScalarType(14) * Y02 + ScalarType(21) * Y11 -
                       ScalarType(56) * Y12) -
                  ScalarType(14) * X00 * Y12 + ScalarType(92) * X01 * Y12 -
                  ScalarType(28) * X02 * Y12 + X00 * Y22 -
                  ScalarType(14) * X01 * Y22 +
                  X11 * (ScalarType(28) * Y00 - ScalarType(42) * Y01 +
                         ScalarType(46) * Y02 + SIX * Y11 -
                         ScalarType(42) * Y12 + ScalarType(28) * Y22)) +
             THREE *
                 (-TWO * BB *
                      (ScalarType(46) * Y02 * Y11 - ScalarType(14) * Y0 * Y111 +
                       Y1111 - ScalarType(14) * Y111 * Y2 -
                       ScalarType(14) * Y01 * Y22 + ScalarType(28) * Y11 * Y22 +
                       Y00 * (ScalarType(28) * Y11 - ScalarType(14) * Y12 +
                              Y22)) +
                  FIVE * Z00 + ScalarType(40) * Z01p - TWO * Z02 +
                  ScalarType(8) * Z1p1p + ScalarType(40) * Z1p2 + FIVE * Z22),
         -TWO * AA *
                 (ScalarType(3) * X02 * X11 - ScalarType(7) * X0 * X111 +
                  THREE * X1111 - ScalarType(7) * X111 * X2) -
             SIX * BB * Y02 * Y11 + ScalarType(14) * BB * Y0 * Y111 -
             SIX * BB * Y1111 -
             TWO * AB *
                 (ScalarType(2) * X12 * Y01 - ScalarType(7) * X01 * Y11 +
                  X02 * Y11 - ScalarType(7) * X12 * Y11 +
                  X11 * (-ScalarType(7) * Y01 + Y02 + SIX * Y11 -
                         ScalarType(7) * Y12) +
                  TWO * X01 * Y12) +
             ScalarType(14) * BB * Y111 * Y2 - FIVE * Z01p + Z02 -
             ScalarType(7) * Z1p1p - FIVE * Z1p2,
         -(AA * X1111) - TWO * AB * X11 * Y11 - BB * Y1111 + Z1p1p});
    for (size_t i = 0; i < 3; ++i) {
      moments.volume() += coeffs[i] * m0_basis[i];
    }
    for (size_t i = 0; i < 4; ++i) {
      moments.centroid()[0] += coeffs[3 + i] * m1x_basis[i];
      moments.centroid()[1] += coeffs[3 + i] * m1y_basis[i];
    }
    for (size_t i = 0; i < 5; ++i) {
      moments.centroid()[2] += coeffs[7 + i] * m1z_basis[i];
    }
    moments.volume() *= area_proj_triangle;
    moments.centroid()[0] *= area_proj_triangle;
    moments.centroid()[1] *= area_proj_triangle;
    moments.centroid()[2] *= area_proj_triangle;
    return moments;
  }
}

template <class ReturnType, class ScalarType>
ReturnType computeFaceOnlyContribution(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PlaneBase<ScalarType>& a_face_plane,
    const PtBase<ScalarType>& a_pt_ref) {
  if constexpr (is_moments_volume<ReturnType>::value) {
    /* Defining constants and types */
    const ScalarType EPSILON = machine_epsilon<ScalarType>();
    const ScalarType ZERO = ScalarType(0);
    const ScalarType FOUR = ScalarType(4);

    /* Function */
    assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
    assert(fabs(a_face_plane.normal()[2]) > EPSILON);
    const ScalarType a =
        -a_face_plane.normal()[0] / safelyEpsilon(a_face_plane.normal()[2]);
    const ScalarType b =
        -a_face_plane.normal()[1] / safelyEpsilon(a_face_plane.normal()[2]);
    const ScalarType c =
        a_face_plane.distance() / safelyEpsilon(a_face_plane.normal()[2]);
    const ScalarType factor = FOUR * a_paraboloid.a() * a_paraboloid.b() * c -
                              a_paraboloid.a() * b * b -
                              a_paraboloid.b() * a * a;
    return ReturnType::fromScalarConstant(copysign(
        machine_pi<ScalarType>() * factor * factor /
            (ScalarType(32) *
             pow(a_paraboloid.a() * a_paraboloid.b(), ScalarType(2.5))),
        -a_face_plane.normal()[2]));
  } else {
    /* Defining constants and types */
    const ScalarType EPSILON = machine_epsilon<ScalarType>();
    const ScalarType ZERO = ScalarType(0);
    const ScalarType FOUR = ScalarType(4);
    const ScalarType FIVE = ScalarType(5);

    /* Function */
    assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
    assert(fabs(a_face_plane.normal()[2]) > EPSILON);
    const ScalarType a = -a_face_plane.normal()[0] / a_face_plane.normal()[2];
    const ScalarType b = -a_face_plane.normal()[1] / a_face_plane.normal()[2];
    const ScalarType c = a_face_plane.distance() / a_face_plane.normal()[2];
    const auto A = a_paraboloid.a(), B = a_paraboloid.b();
    auto moments = ReturnType::fromScalarConstant(ZERO);
    const ScalarType factor = (a * a * B + A * (b * b - FOUR * B * c)) *
                              (a * a * B + A * (b * b - FOUR * B * c)) *
                              machine_pi<ScalarType>();
    moments.volume() =
        copysign(factor / (ScalarType(32) * pow(A * B, ScalarType(2.5))),
                 -a_face_plane.normal()[2]);
    moments.centroid()[0] =
        a * B *
        copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]);
    moments.centroid()[1] =
        b * A *
        copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]);
    moments.centroid()[2] =
        (FIVE * A * (b * b) + FIVE * (a * a) * B - ScalarType(8) * A * B * c) *
        copysign(factor / (ScalarType(384) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]);
    return moments;
  }
}

template <class ReturnType, class ScalarType>
ReturnType computeTriangleCorrection(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1,
    const PtBase<ScalarType>& a_pt_2) {
  if constexpr (is_moments_volume<ReturnType>::value) {
    return ReturnType::fromScalarConstant(
        (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) * (a_pt_1[0] + a_pt_2[0]) +
         -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) * (a_pt_1[1] + a_pt_2[1]) -
         a_pt_0[2] - ScalarType(2) * a_pt_1[2] - a_pt_2[2]) /
        ScalarType(12) *
        ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
         (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
         (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0]));
  } else {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType FOUR = ScalarType(4);
    const ScalarType FIVE = ScalarType(5);
    const ScalarType SIX = ScalarType(6);
    const ScalarType HALF = ONE / TWO;

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ZERO);
    const ScalarType A = a_paraboloid.a(), B = a_paraboloid.b();
    const ScalarType X0 = a_pt_0[0], X1 = a_pt_1[0], X2 = a_pt_2[0];
    const ScalarType Y0 = a_pt_0[1], Y1 = a_pt_1[1], Y2 = a_pt_2[1];
    const ScalarType Z0 = a_pt_0[2], Z1 = a_pt_1[2], Z2 = a_pt_2[2];
    const ScalarType triangle_area =
        HALF * ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
                (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
                (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0]);
    moments.volume() =
        (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) * (a_pt_1[0] + a_pt_2[0]) +
         -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) * (a_pt_1[1] + a_pt_2[1]) -
         a_pt_0[2] - TWO * a_pt_1[2] - a_pt_2[2]) *
        triangle_area / SIX;
    moments.centroid()[0] =
        triangle_area *
        ((A * (-(X0 * X0 * X0) - X1 * X1 * X1 - X1 * X1 * X2 - X1 * (X2 * X2) -
               X2 * X2 * X2 - X0 * X0 * (X1 + X2) -
               X0 * (X1 * X1 + X1 * X2 + X2 * X2))) /
             ScalarType(10) +
         (B * (-(X1 * (Y0 * Y0 + TWO * Y0 * Y1 + THREE * (Y1 * Y1) + Y0 * Y2 +
                       TWO * Y1 * Y2 + Y2 * Y2)) -
               X2 * (Y0 * Y0 + Y0 * Y1 + Y1 * Y1 + TWO * Y0 * Y2 +
                     TWO * Y1 * Y2 + THREE * (Y2 * Y2)) -
               X0 * (THREE * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                     TWO * Y0 * (Y1 + Y2)))) /
             ScalarType(30) +
         (-(X0 * (TWO * Z0 + Z1 + Z2)) - X1 * (Z0 + TWO * Z1 + Z2) -
          X2 * (Z0 + Z1 + TWO * Z2)) /
             ScalarType(12));
    moments.centroid()[1] =
        -triangle_area *
        ((B * (Y0 * Y0 * Y0 + Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) +
               Y2 * Y2 * Y2 + Y0 * Y0 * (Y1 + Y2) +
               Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2))) /
             ScalarType(10) +
         (A *
          (X0 * X0 * (THREE * Y0 + Y1 + Y2) + X1 * X1 * (Y0 + THREE * Y1 + Y2) +
           X2 * X2 * (Y0 + Y1 + THREE * Y2) + X1 * X2 * (Y0 + TWO * (Y1 + Y2)) +
           X0 * (X1 * (TWO * Y0 + TWO * Y1 + Y2) +
                 X2 * (TWO * Y0 + Y1 + TWO * Y2)))) /
             ScalarType(30) +
         (Y0 * (TWO * Z0 + Z1 + Z2) + Y1 * (Z0 + TWO * Z1 + Z2) +
          Y2 * (Z0 + Z1 + TWO * Z2)) /
             ScalarType(12));
    moments.centroid()[2] =
        triangle_area *
        ((A * A *
          (X0 * X0 * X0 * X0 + X1 * X1 * X1 * X1 + X1 * X1 * X1 * X2 +
           X1 * X1 * (X2 * X2) + X1 * (X2 * X2 * X2) + X2 * X2 * X2 * X2 +
           X0 * X0 * X0 * (X1 + X2) + X0 * X0 * (X1 * X1 + X1 * X2 + X2 * X2) +
           X0 *
               (X1 * X1 * X1 + X1 * X1 * X2 + X1 * (X2 * X2) + X2 * X2 * X2))) /
             ScalarType(30) +
         (B * B *
          (Y0 * Y0 * Y0 * Y0 + Y1 * Y1 * Y1 * Y1 + Y1 * Y1 * Y1 * Y2 +
           Y1 * Y1 * (Y2 * Y2) + Y1 * (Y2 * Y2 * Y2) + Y2 * Y2 * Y2 * Y2 +
           Y0 * Y0 * Y0 * (Y1 + Y2) + Y0 * Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2) +
           Y0 *
               (Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) + Y2 * Y2 * Y2))) /
             ScalarType(30) +
         (A * B *
          (X1 * X2 *
               (Y0 * Y0 + THREE * (Y1 * Y1) + FOUR * Y1 * Y2 +
                THREE * (Y2 * Y2) + TWO * Y0 * (Y1 + Y2)) +
           X0 * X0 *
               (SIX * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                THREE * Y0 * (Y1 + Y2)) +
           X1 * X1 *
               (Y0 * Y0 + SIX * (Y1 * Y1) + THREE * Y1 * Y2 + Y2 * Y2 +
                Y0 * (THREE * Y1 + Y2)) +
           X2 * X2 *
               (Y0 * Y0 + Y1 * Y1 + THREE * Y1 * Y2 + SIX * (Y2 * Y2) +
                Y0 * (Y1 + THREE * Y2)) +
           X0 * (X1 * (THREE * (Y0 * Y0) + FOUR * Y0 * Y1 + THREE * (Y1 * Y1) +
                       TWO * Y0 * Y2 + TWO * Y1 * Y2 + Y2 * Y2) +
                 X2 * (THREE * (Y0 * Y0) + TWO * Y0 * Y1 + Y1 * Y1 +
                       FOUR * Y0 * Y2 + TWO * Y1 * Y2 + THREE * (Y2 * Y2))))) /
             ScalarType(90) +
         (-(Z0 * Z0) - Z1 * Z1 - Z1 * Z2 - Z2 * Z2 - Z0 * (Z1 + Z2)) /
             ScalarType(12));
    return moments;
  }
}

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_
