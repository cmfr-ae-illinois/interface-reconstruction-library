// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_MOMENT_CONTRIBUTIONS_TPP_
#define IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_MOMENT_CONTRIBUTIONS_TPP_

#include <float.h>
#include <cassert>
#include <cmath>

#include "external/NumericalIntegration/NumericalIntegration.h"
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
#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/generic_cutting/quadratic_intersection/moment_contributions.h"

// this enable a lot of debug text when computing
//#define VALDEBUG

namespace IRL {

template <class ScalarType, UnsignedIndex_t OrderX, UnsignedIndex_t OrderY,
          UnsignedIndex_t OrderZ>
inline ScalarType MomentCylinderIntegrand(
    const PtBase<ScalarType>& a_position,
    const PtBase<ScalarType>& a_derivative, const ScalarType& B,
    const ScalarType& R) {
  const ScalarType &x = a_position[0], &y = a_position[1],
                   &dy = a_derivative[1];
  const ScalarType inv2 = CstHalf<ScalarType>();
  const ScalarType inv3 = CstThird<ScalarType>();
  const ScalarType inv4 = CstFourth<ScalarType>();
  const ScalarType inv5 = CstFifth<ScalarType>();
  const ScalarType inv6 = CstSixth<ScalarType>();
  const ScalarType inv7 = CstSeventh<ScalarType>();

  // M0
  if constexpr (OrderX == 0 && OrderY == 0 && OrderZ == 0) {
    return -dy * x * sqrt(R - B * y * y);
  }
  // M1x
  else if constexpr (OrderX == 1 && OrderY == 0 && OrderZ == 0) {
    const auto x2 = x * x;
    return -dy * inv2 * x2 * sqrt(R - B * y * y);
  }
  // M1y
  else if constexpr (OrderX == 0 && OrderY == 1 && OrderZ == 0) {
    return -dy * x * y * sqrt(R - B * y * y);
  }
  // M1z
  else if constexpr (OrderX == 0 && OrderY == 0 && OrderZ == 1) {
    const auto y2 = y * y;
    return dy * inv2 * x * (R - B * y2);
  }
  return ScalarType(0);
}

template <UnsignedIndex_t ProjDir, class ReturnType, class ScalarType,
          class MomentFunctorType>
inline ReturnType MomentsIntegrandCylinderArc(const ScalarType a_t,
                                      const MomentFunctorType& a_functor) {
  using ReturnScalarType = typename ReturnType::value_type;
  const auto der_t = a_functor.der_t(a_t);
  const auto pos_t = a_functor.pos_t(a_t);
  const auto& normal = a_functor.face_normal();
  const auto& dist = a_functor.face_distance();
  const auto& weight = a_functor.norm_weight();
  const auto& b = a_functor.b();
  const auto& r = a_functor.r();
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    return ReturnType::fromScalarConstant(ReturnScalarType(
        MomentCylinderIntegrand<ScalarType, 0, 0, 0>(pos_t, der_t, b, r) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight)));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
    auto integrand = ReturnType::fromScalarConstant(ReturnScalarType(0));
    integrand.volume() = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 0, 0>(pos_t, der_t, b, r) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand.centroid()[0] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 1, 0, 0>(pos_t, der_t, b, r) -
        MomentPlaneIntegrand<ScalarType, 1, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand.centroid()[1] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 1, 0>(pos_t, der_t, b, r) -
        MomentPlaneIntegrand<ScalarType, 0, 1, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand.centroid()[2] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 0, 1>(pos_t, der_t, b, r) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 1, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    return integrand;
  } else {
    std::cout << "Cylinder M>1 not available yet" << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

template <class ReturnType, class ScalarType, UnsignedIndex_t QuadRuleOrder>
class CylinderMomentArcIntegrator {
 public:
  using ReturnScalarType = typename ReturnType::value_type;
  using AlignedCylinder = AlignedCylinderBase<ScalarType>;
  using RationalBezierArc = RationalBezierArcBase<ScalarType>;
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  CylinderMomentArcIntegrator(const AlignedCylinder& a_cylinder,
                                const RationalBezierArc& a_arc,
                                const Normal& a_face_normal,
                                const UnsignedIndex_t a_proj_dir)
      : lower_limit_m(ScalarType(0)),
        upper_limit_m(ScalarType(1)),
        arc_m(a_arc),
        face_normal_m(a_face_normal),
        b_m(a_cylinder.b()),
        r_m(a_cylinder.r()),
        face_distance_m(a_face_normal * a_arc.start_point()) {
    switch (a_proj_dir) {
      case 0:
        integrand_m = MomentsIntegrandCylinderArc<0>;
        norm_weight_m = face_normal_m[2] / face_normal_m[0];
        break;
      case 1:
        integrand_m = MomentsIntegrandCylinderArc<1>;
        norm_weight_m = -face_normal_m[2] / face_normal_m[1];
        break;
      case 2:
        integrand_m = MomentsIntegrandCylinderArc<2>;
        norm_weight_m = ScalarType(1) / face_normal_m[2];
        break;
      default:
        UnsignedIndex_t max_component_index = 0;
        ScalarType max_component = fabs(face_normal_m[0]);
        for (UnsignedIndex_t d = 1; d < 3; ++d) {
          if (fabs(face_normal_m[d]) > max_component) {
            max_component_index = d;
            max_component = fabs(face_normal_m[d]);
          }
        }
        switch (max_component_index) {
          case 0:
            integrand_m = MomentsIntegrandCylinderArc<0>;
            norm_weight_m = face_normal_m[2] / face_normal_m[0];
            break;
          case 1:
            integrand_m = MomentsIntegrandCylinderArc<1>;
            norm_weight_m = -face_normal_m[2] / face_normal_m[1];
            break;
          default:
            integrand_m = MomentsIntegrandCylinderArc<2>;
            norm_weight_m = ScalarType(1) / face_normal_m[2];
        }
    }
  }
  const ReturnType operator()(const ScalarType a_t) const {
    return (*integrand_m)(a_t, (*this));
  }
  const ReturnType integrate(void) {
    const auto one = ScalarType(1), inv2 = ScalarType(0.5);
    const auto& abscissea = AbscissaeGauss<ScalarType, QuadRuleOrder>();
    const auto& weights = WeightsGauss<ScalarType, QuadRuleOrder>();
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(0));
    for (UnsignedIndex_t i = 0; i < QuadRuleOrder; ++i) {
      const auto t = lower_limit_m + (upper_limit_m - lower_limit_m) * inv2 *
                                         (one + abscissea[i]);
      moments += ReturnScalarType(inv2 * weights[i]) * (*this)(t);
    }
    return moments;
  }
  inline const Pt der_t(const ScalarType a_t) const {
    return arc_m.derivative(a_t);
  }
  inline const Pt pos_t(const ScalarType a_t) const { return arc_m.point(a_t); }
  inline const Normal& face_normal(void) const { return face_normal_m; }
  inline const ScalarType& face_distance(void) const { return face_distance_m; }
  inline const ScalarType& norm_weight(void) const { return norm_weight_m; }
  inline const ScalarType& b(void) const { return b_m; }
  inline const ScalarType& r(void) const { return r_m; }

 private:
  const RationalBezierArc arc_m;
  const Normal face_normal_m;
  const ScalarType face_distance_m, b_m, r_m;
  ScalarType norm_weight_m;
  ScalarType lower_limit_m, upper_limit_m;
  ReturnType (*integrand_m)(const ScalarType,
                            const CylinderMomentArcIntegrator&);
};

/******************************************************************************/
/*********************** First moment contribution
 * ****************************/
/******************************************************************************/
/* This compute the first contribution to the moments (arising from the
 * integration of the face plane primitives on the poligonized clipped faces)
 */

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsVC3SeriesOne(
    const ContainerType& a_weight) {
  std::array<ContainerType, 3> coeffs;
  coeffs.fill(ContainerType(ScalarType(0)));
  coeffs[0]= ScalarType(1)/ScalarType(6);
  ContainerType x(1);
  UnsignedIndex_t i = 0;
  ScalarType max_diff;
  while (i <= 40) {
    max_diff = ScalarType(0);
    for (UnsignedIndex_t j = 0; j < 2; ++j) {
      ContainerType add_to_coeff = ScalarType(vc3Series[i][j]) * x;
      coeffs[j+1] += add_to_coeff;
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

// template <class ContainerType, class ScalarType>
// inline std::array<ContainerType, 12> coeffsV3andC3SeriesOne(
//     const ContainerType& a_weight) {
//   std::array<ContainerType, 12> coeffs;
//   coeffs.fill(ContainerType(ScalarType(0)));
//   ContainerType x(ScalarType(1));
//   UnsignedIndex_t i = 0;
//   ScalarType max_diff;
//   while (i <= 40) {
//     max_diff = ScalarType(0);
//     for (UnsignedIndex_t j = 0; j < 3; ++j) {
//       ContainerType add_to_coeff = ScalarType(v3Series[i][j]) * x;
//       coeffs[j] += add_to_coeff;
//       max_diff = maximum(max_diff, fabs(add_to_coeff));
//     }
//     if (max_diff < ScalarType(DBL_EPSILON)) {
//       break;
//     }
//     x *= a_weight - ContainerType(ScalarType(1));
//     i++;
//   }
//   i = 0;
//   x = ContainerType(ScalarType(1));
//   while (i <= 40) {
//     max_diff = ScalarType(0);
//     for (UnsignedIndex_t j = 0; j < 4; ++j) {
//       ContainerType add_to_coeff = ScalarType(cx3Series[i][j]) * x;
//       coeffs[3 + j] += add_to_coeff;
//       max_diff = maximum(max_diff, fabs(add_to_coeff));
//     }
//     for (UnsignedIndex_t j = 0; j < 5; ++j) {
//       ContainerType add_to_coeff = ScalarType(cz3Series[i][j]) * x;
//       coeffs[7 + j] += add_to_coeff;
//       max_diff = maximum(max_diff, fabs(add_to_coeff));
//     }
//     if (max_diff < ScalarType(DBL_EPSILON)) {
//       break;
//     }
//     x *= a_weight - ContainerType(ScalarType(1));
//     i++;
//   }

//   return coeffs;
// }

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsVC3Exact(
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
  const auto L = ContainerType(ONE) / ((a_weight - ContainerType(ONE)) *
                                       (a_weight + ContainerType(ONE)));
  const auto L2 = L * L;
  const auto S = (a_weight < ContainerType(ONE))
                     ? sqrt(ContainerType(ONE) - a_weight * a_weight)
                     : sqrt(a_weight * a_weight - ContainerType(ONE));
  const auto T = (a_weight < ContainerType(ONE))
                     ? atan((ContainerType(ONE) - a_weight) / S) / S
                     : atanh((a_weight - ContainerType(ONE)) / S) / S;
  return std::array<ContainerType, 3>(
      { ONE / ScalarType(6),
       (TWO * w4 - ScalarType(5) * w2 + 
         SIX * a_weight * T) *
           L2 / ScalarType(12),
       (THREE * w2 -
        (FOUR * w3 + TWO * a_weight) * T) *
           L2 / FOUR});
}

// template <class ContainerType, class ScalarType>
// inline std::array<ContainerType, 12> coeffsV3andC3Exact(
//     const ContainerType& a_weight) {
//   /* Defining constants and types */
//   const ScalarType ONE = ScalarType(1);
//   const ScalarType TWO = ScalarType(2);
//   const ScalarType THREE = ScalarType(3);
//   const ScalarType FOUR = ScalarType(4);
//   const ScalarType FIVE = ScalarType(5);
//   const ScalarType SIX = ScalarType(6);

//   /* Function */
//   const auto w2 = a_weight * a_weight;
//   const auto w3 = w2 * a_weight;
//   const auto w4 = w2 * w2;
//   const auto w5 = w2 * w3;
//   const auto w6 = w3 * w3;
//   const auto w7 = w4 * w3;
//   const auto w8 = w4 * w4;
//   const auto w9 = w4 * w5;
//   const auto w10 = w5 * w5;
//   const auto L = ContainerType(ONE) / ((a_weight - ContainerType(ONE)) *
//                                        (a_weight + ContainerType(ONE)));
//   const auto L3 = L * L * L;
//   const auto L4 = L3 * L;
//   const auto L5 = L4 * L;
//   const auto S = (a_weight < ContainerType(ONE))
//                      ? sqrt(ContainerType(ONE) - a_weight * a_weight)
//                      : sqrt(a_weight * a_weight - ContainerType(ONE));
//   const auto T = (a_weight < ContainerType(ONE))
//                      ? atan((ContainerType(ONE) - a_weight) / S) / S
//                      : atanh((a_weight - ContainerType(ONE)) / S) / S;
//   return std::array<ContainerType, 12>(
//       {L3 *
//            (TWO * w6 - THREE * w4 + ScalarType(31) * w2 -
//             (ScalarType(42) * w3 + ScalarType(18) * a_weight) * T) /
//            ScalarType(12),
//        L3 *
//            (TWO * w6 - ScalarType(9) * w4 - ScalarType(8) * w2 +
//             ScalarType(30) * w3 * T) /
//            THREE,
//        L3 *
//            (ScalarType(11) * w4 + FOUR * w2 -
//             (ScalarType(12) * w5 + ScalarType(18) * w3) * T) /
//            SIX,
//        L4 * ((-T * a_weight) / ScalarType(32) +
//              (ScalarType(93) * (w2)) / ScalarType(2240) -
//              (ScalarType(163) * (w4)) / ScalarType(3360) +
//              (FIVE * (w6)) / ScalarType(168) - (w8) / ScalarType(140)),
//        L4 * ((w2) / ScalarType(70) + (-T * (w3)) / ScalarType(16) +
//              (ScalarType(29) * (w4)) / ScalarType(1120) -
//              (ScalarType(19) * (w6)) / ScalarType(1680) +
//              (w8) / ScalarType(420)),
//        -L4 *
//            ((w2) / ScalarType(210) - (w4) / ScalarType(21) -
//             (-T * (w5)) / ScalarType(8) -
//             (ScalarType(13) * (w6)) / ScalarType(560) + (w8) / ScalarType(280)),
//        L4 * ((w2) / ScalarType(35) - (ScalarType(16) * (w4)) / ScalarType(105) +
//              (ScalarType(58) * (w6)) / ScalarType(105) - T * (w7) +
//              (w8) / ScalarType(14)),
//        L5 * ((-T * a_weight) / ScalarType(128) +
//              (ScalarType(193) * (w2)) / ScalarType(16128) -
//              (ScalarType(149) * (w4)) / ScalarType(8064) +
//              (ScalarType(19) * (w6)) / ScalarType(1120) -
//              (ScalarType(41) * (w8)) / ScalarType(5040) +
//              (w10) / ScalarType(630)),
//        L5 * ((FOUR * (w2)) / ScalarType(945) + (-T * (w3)) / ScalarType(48) +
//              (ScalarType(65) * (w4)) / ScalarType(6048) -
//              (w6) / ScalarType(144) +
//              (ScalarType(11) * (w8)) / ScalarType(3780) -
//              (w10) / ScalarType(1890)),
//        -L5 * ((w2) / ScalarType(1890) -
//               (ScalarType(13) * (w4)) / ScalarType(1890) -
//               (-T * (w5)) / ScalarType(48) -
//               (ScalarType(11) * (w6)) / ScalarType(2016) +
//               (ScalarType(5) * (w8)) / ScalarType(3024) -
//               (w10) / ScalarType(3780)),
//        L5 * ((w2) / ScalarType(315) - (w4) / ScalarType(45) +
//              (ScalarType(4) * (w6)) / ScalarType(35) +
//              (-T * (w7)) / ScalarType(4) +
//              (ScalarType(17) * (w8)) / ScalarType(504) -
//              (w10) / ScalarType(252)),
//        -L5 *
//            ((w2) / ScalarType(63) - (ScalarType(29) * (w4)) / ScalarType(315) +
//             (ScalarType(26) * (w6)) / ScalarType(105) -
//             (ScalarType(194) * (w8)) / ScalarType(315) + T * (w9) -
//             (w10) / ScalarType(18))});
// }

template <class ReturnType, class ScalarType>
ReturnType computeType3Contribution(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const RationalBezierArcBase<ScalarType>& a_arc,
    const NormalBase<ScalarType>& a_face_normal) {
  using ReturnScalarType = typename ReturnType::value_type;
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
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


    #ifdef VALDEBUG
    std::cout << "M3 computation :" << std::endl;
    std::cout << "x0 : " << pt_0 << std::endl;
    std::cout << "x* : " << cp << std::endl;
    std::cout << "x1 : " << pt_1 << std::endl;
    std::cout << "w : " << weight << std::endl;
    #endif

    assert(weight >= ZERO);
    std::array<ScalarType, 3> coeffs;
    if (weight < ScalarType(0.35)) { // We use the exact expressions
    #ifdef VALDEBUG
    std::cout << "weight is low, using exact value" << std::endl;
    #endif
      coeffs = coeffsVC3Exact<ScalarType, ScalarType>(weight);
    } else if (weight <
             ScalarType(1.7)) { // We use the 40th order Taylor series (w -> 1)
    #ifdef VALDEBUG
    std::cout << "weight is close to 1, using Taylor series" << std::endl;
    #endif
      coeffs = coeffsVC3SeriesOne<ScalarType, ScalarType>(weight);
    } else if (weight <
             ScalarType(1.0e9)) { // We use the series expansion (w -> infty)
    #ifdef VALDEBUG
    std::cout << "weight is high, using exact value" << std::endl;
    #endif
      coeffs = coeffsVC3Exact<ScalarType, ScalarType>(weight);
    } else { // This is within EPSILON of the actual value
    #ifdef VALDEBUG
    std::cout << "weight is BIG, using limite" << std::endl;
    #endif
      coeffs = std::array<ScalarType, 3>({ONE / SIX, ONE / SIX, ZERO});
    }

    #ifdef VALDEBUG
    std::cout << "coeffs :" << std::endl;
    std::cout << "0 : " << coeffs[0] << std::endl;
    std::cout << "1 : " << coeffs[1] << std::endl;
    std::cout << "2 : " << coeffs[2] << std::endl;
    #endif

    const auto& C11 = (pt_1[1]-pt_0[1]) * 
        (TWO * pt_0[0] * pt_0[2] + TWO * pt_1[0] * pt_1[2] +
          pt_0[0] * pt_1[2] + pt_1[0] * pt_0[2]);
    const auto& C12 = (pt_0[0] + pt_1[0] + cp[0]) *
        (pt_0[1] * (pt_1[2] - cp[2]) + pt_1[1] * (cp[2] - pt_0[2]) + 
         cp[1] * (pt_0[2] - pt_1[2]));
    const auto& C13 = cp[0] * 
        (pt_0[1] * (pt_1[2] - cp[2]) + pt_1[1] * (cp[2] - pt_0[2]) + 
         cp[1] * (pt_0[2] - pt_1[2]));
    #ifdef VALDEBUG
    std::cout << "C vector :" << std::endl;
    std::cout << "C11 : " << C11 << std::endl;
    std::cout << "C12 : " << C12 << std::endl;
    std::cout << "C13 : " << C13 << std::endl;
    #endif
    return ReturnType::fromScalarConstant(ReturnScalarType(
        (coeffs[0] * C11 +
         coeffs[1] * C12 +
         coeffs[2] * C13)));
  // } else if constexpr (std::is_same_v<ReturnType,
  //                                     VolumeMomentsBase<ReturnScalarType>>) {
  //   /* Defining constants and types */
  //   const ScalarType ZERO = ScalarType(0);
  //   const ScalarType ONE = ScalarType(1);
  //   const ScalarType TWO = ScalarType(2);
  //   const ScalarType THREE = ScalarType(3);
  //   const ScalarType FOUR = ScalarType(4);
  //   const ScalarType FIVE = ScalarType(5);
  //   const ScalarType SIX = ScalarType(6);
  //   const ScalarType HALF = ONE / TWO;
  //   const ScalarType QUARTER = ONE / FOUR;

  //   /* Function */
  //   auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
  //   const auto& pt_0 = a_arc.start_point();
  //   const auto& cp = a_arc.control_point();
  //   const auto& pt_1 = a_arc.end_point();
  //   const auto weight = a_arc.weight();
  //   const auto A = a_paraboloid.a(), B = a_paraboloid.b();
  //   const auto X0 = pt_0[0], X1 = cp[0], X2 = pt_1[0];
  //   const auto Y0 = pt_0[1], Y1 = cp[1], Y2 = pt_1[1];
  //   const auto Z0 = pt_0[2], Z1p = cp[2], Z2 = pt_1[2];
  //   const ScalarType Z1P = -A * X1 * X1 - B * Y1 * Y1;
  //   const ScalarType AA = A * A, BB = B * B, AB = A * B;
  //   const ScalarType X00 = X0 * X0, X11 = X1 * X1, X22 = X2 * X2;
  //   const ScalarType X000 = X00 * X0, X111 = X11 * X1, X222 = X22 * X2;
  //   const ScalarType X0000 = X00 * X00, X1111 = X11 * X11, X2222 = X22 * X22;
  //   const ScalarType Y00 = Y0 * Y0, Y11 = Y1 * Y1, Y22 = Y2 * Y2;
  //   const ScalarType Y000 = Y00 * Y0, Y111 = Y11 * Y1, Y222 = Y22 * Y2;
  //   const ScalarType Y0000 = Y00 * Y00, Y1111 = Y11 * Y11, Y2222 = Y22 * Y22;
  //   const ScalarType Z00 = Z0 * Z0, Z22 = Z2 * Z2;
  //   const ScalarType Z1p1p = Z1p * Z1p, Z1P1P = Z1P * Z1P;
  //   const ScalarType X02 = X0 * X2, X12 = X1 * X2, X01 = X0 * X1;
  //   const ScalarType Y02 = Y0 * Y2, Y12 = Y1 * Y2, Y01 = Y0 * Y1;
  //   const ScalarType Z02 = Z0 * Z2, Z01p = Z0 * Z1p, Z01P = Z0 * Z1P;
  //   const ScalarType Z1p1P = Z1p * Z1P, Z1p2 = Z1p * Z2, Z1P2 = Z1P * Z2;
  //   const ScalarType X0Z0 = X0 * Z0, X0Z1p = X0 * Z1p, X0Z1P = X0 * Z1P,
  //                    X0Z2 = X0 * Z2;
  //   const ScalarType X1Z0 = X1 * Z0, X1Z1p = X1 * Z1p, X1Z1P = X1 * Z1P,
  //                    X1Z2 = X1 * Z2;
  //   const ScalarType X2Z0 = X2 * Z0, X2Z1p = X2 * Z1p, X2Z1P = X2 * Z1P,
  //                    X2Z2 = X2 * Z2;
  //   const ScalarType Y0Z0 = Y0 * Z0, Y0Z1p = Y0 * Z1p, Y0Z1P = Y0 * Z1P,
  //                    Y0Z2 = Y0 * Z2;
  //   const ScalarType Y1Z0 = Y1 * Z0, Y1Z1p = Y1 * Z1p, Y1Z1P = Y1 * Z1P,
  //                    Y1Z2 = Y1 * Z2;
  //   const ScalarType Y2Z0 = Y2 * Z0, Y2Z1p = Y2 * Z1p, Y2Z1P = Y2 * Z1P,
  //                    Y2Z2 = Y2 * Z2;
  //   const ScalarType area_proj_triangle =
  //       HALF * (X0 * (Y2 - Y1) + X2 * (Y1 - Y0) + X1 * (Y0 - Y2));
  //   assert(weight >= ZERO);
  //   // Compute coefficients (functions of the weight)
  //   std::array<ScalarType, 12> coeffs;
  //   if (weight < ScalarType(0.35))  // We use the exact expressions
  //   {
  //     coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
  //   } else if (weight <
  //              ScalarType(1.7))  // We use the 40th order Taylor series (w -> 1)
  //   {
  //     coeffs = coeffsV3andC3SeriesOne<ScalarType, ScalarType>(weight);
  //   } else if (weight < ScalarType(1.0e9))  // We use the exact expressions
  //   {
  //     coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
  //   }
  //   // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
  //   //   coeffs = coeffsV3andC3SeriesInfinity(weight);
  //   else  // This is within EPSILON of the actual value
  //   {
  //     coeffs = std::array<ScalarType, 12>(
  //         {ScalarType(ONE / SIX), ScalarType(TWO / THREE), ScalarType(ZERO),
  //          ScalarType(-ONE / ScalarType(140)),
  //          ScalarType(ONE / ScalarType(420)),
  //          ScalarType(-ONE / ScalarType(280)), ScalarType(ONE / ScalarType(14)),
  //          ScalarType(ONE / ScalarType(630)),
  //          ScalarType(-ONE / ScalarType(1890)),
  //          ScalarType(ONE / ScalarType(3780)),
  //          ScalarType(-ONE / ScalarType(252)),
  //          ScalarType(ONE / ScalarType(18))});
  //   }
  //   auto m0_basis = std::array<ScalarType, 3>(
  //       {signedDistance<ScalarType>(HALF * (pt_0 + pt_1), a_paraboloid),
  //        signedDistance<ScalarType>(QUARTER * (pt_0 + pt_1) + HALF * cp,
  //                                   a_paraboloid),
  //        signedDistance<ScalarType>(cp, a_paraboloid)});
  //   auto m1x_basis = std::array<ScalarType, 4>(
  //       {-SIX * (X0Z0 - X0Z2 - X2Z0 + X2Z2 - TWO * B * X2 * Y00 +
  //                TWO * B * X0 * Y02 + TWO * B * X2 * Y02 - TWO * B * X0 * Y22),
  //        TWO * (FIVE * X0Z0 + ScalarType(10) * X0Z1p + SIX * X0Z1P +
  //               ScalarType(7) * X0Z2 + ScalarType(30) * A * X02 * X1 -
  //               ScalarType(11) * X1Z0 - FOUR * X1Z1p - ScalarType(11) * X1Z2 +
  //               ScalarType(7) * X2Z0 + ScalarType(10) * X2Z1p + SIX * X2Z1P +
  //               FIVE * X2Z2 - ScalarType(14) * B * X1 * Y00 +
  //               FOUR * B * X2 * Y00 + ScalarType(14) * B * X0 * Y01 -
  //               FOUR * B * X1 * Y01 + ScalarType(10) * B * X2 * Y01 -
  //               FOUR * B * X0 * Y02 + ScalarType(10) * B * X1 * Y02 -
  //               FOUR * B * X2 * Y02 + FOUR * B * X0 * Y11 +
  //               FOUR * B * X2 * Y11 + ScalarType(10) * B * X0 * Y12 -
  //               FOUR * B * X1 * Y12 + ScalarType(14) * B * X2 * Y12 +
  //               FOUR * B * X0 * Y22 - ScalarType(14) * B * X1 * Y22),
  //        TWO * (-FIVE * X0Z1p + ScalarType(18) * X0Z1P + X0Z2 +
  //               SIX * A * X02 * X1 - FIVE * X1Z0 - SIX * X1Z1p - SIX * X1Z1P -
  //               FIVE * X1Z2 + X2Z0 - FIVE * X2Z1p + ScalarType(18) * X2Z1P -
  //               ScalarType(12) * B * X1 * Y01 + TWO * B * X2 * Y01 +
  //               TWO * B * X1 * Y02 + ScalarType(12) * B * X0 * Y11 +
  //               ScalarType(12) * B * X2 * Y11 + TWO * B * X0 * Y12 -
  //               ScalarType(12) * B * X1 * Y12),
  //        TWO * (X1Z1p - X1Z1P)});
  //   auto m1y_basis = std::array<ScalarType, 4>(
  //       {SIX *
  //            (-Y0Z0 + Y0Z2 + TWO * A * (X22 * Y0 + X00 * Y2 - X02 * (Y0 + Y2)) +
  //             Y2Z0 - Y2Z2),
  //        TWO *
  //            (FIVE * Y0Z0 + ScalarType(10) * Y0Z1p + SIX * Y0Z1P +
  //             ScalarType(7) * Y0Z2 + ScalarType(30) * B * Y02 * Y1 -
  //             ScalarType(11) * Y1Z0 - FOUR * Y1Z1p - ScalarType(11) * Y1Z2 +
  //             TWO * A *
  //                 (-TWO * X02 * Y0 + TWO * X11 * Y0 + FIVE * X12 * Y0 +
  //                  TWO * X22 * Y0 - ScalarType(7) * X00 * Y1 + FIVE * X02 * Y1 -
  //                  TWO * X12 * Y1 - ScalarType(7) * X22 * Y1 + TWO * X00 * Y2 -
  //                  TWO * X02 * Y2 + TWO * X11 * Y2 + ScalarType(7) * X12 * Y2 +
  //                  X01 * (ScalarType(7) * Y0 - TWO * Y1 + FIVE * Y2)) +
  //             ScalarType(7) * Y2Z0 + ScalarType(10) * Y2Z1p + SIX * Y2Z1P +
  //             FIVE * Y2Z2),
  //        -TWO * (FIVE * Y0Z1p - ScalarType(18) * Y0Z1P - Y0Z2 -
  //                SIX * B * Y02 * Y1 + FIVE * Y1Z0 + SIX * Y1Z1p + SIX * Y1Z1P +
  //                FIVE * Y1Z2 -
  //                TWO * A *
  //                    (X12 * Y0 - SIX * X01 * Y1 + X02 * Y1 - SIX * X12 * Y1 +
  //                     X01 * Y2 + SIX * X11 * (Y0 + Y2)) -
  //                Y2Z0 + FIVE * Y2Z1p - ScalarType(18) * Y2Z1P),
  //        TWO * (Y1Z1p - Y1Z1P)});
  //   auto m1z_basis = std::array<ScalarType, 5>(
  //       {-(AA * (ScalarType(21) * X0000 + ScalarType(28) * X000 * X2 +
  //                ScalarType(30) * X00 * X22 + ScalarType(28) * X0 * X222 +
  //                ScalarType(21) * X2222)) -
  //            ScalarType(21) * BB * Y0000 - ScalarType(28) * BB * Y000 * Y2 -
  //            ScalarType(30) * BB * Y00 * Y22 -
  //            TWO * AB *
  //                (X00 * (ScalarType(21) * Y00 + ScalarType(14) * Y02 +
  //                        FIVE * Y22) +
  //                 TWO * X02 *
  //                     (ScalarType(7) * Y00 + ScalarType(10) * Y02 +
  //                      ScalarType(7) * Y22) +
  //                 X22 * (ScalarType(5) * Y00 + ScalarType(14) * Y02 +
  //                        ScalarType(21) * Y22)) -
  //            ScalarType(28) * BB * Y0 * Y222 - ScalarType(21) * BB * Y2222 +
  //            ScalarType(40) * Z00 + ScalarType(48) * Z02 + ScalarType(40) * Z22,
  //        THREE * AA *
  //                (ScalarType(21) * X000 * X1 - ScalarType(7) * X00 * X11 -
  //                 ScalarType(10) * X02 * X11 + ScalarType(35) * X00 * X12 -
  //                 ScalarType(7) * X000 * X2 - ScalarType(10) * X00 * X22 +
  //                 ScalarType(35) * X01 * X22 - ScalarType(7) * X11 * X22 -
  //                 ScalarType(7) * X0 * X222 + ScalarType(21) * X1 * X222) -
  //            AB * (ScalarType(7) * X11 * Y00 - ScalarType(35) * X12 * Y00 +
  //                  ScalarType(10) * X22 * Y00 - ScalarType(63) * X00 * Y01 +
  //                  ScalarType(20) * X12 * Y01 - ScalarType(35) * X22 * Y01 +
  //                  ScalarType(21) * X00 * Y02 + ScalarType(10) * X11 * Y02 -
  //                  ScalarType(70) * X12 * Y02 + ScalarType(21) * X22 * Y02 +
  //                  ScalarType(7) * X00 * Y11 + ScalarType(7) * X22 * Y11 -
  //                  ScalarType(35) * X00 * Y12 + ScalarType(28) * X12 * Y12 -
  //                  ScalarType(63) * X22 * Y12 +
  //                  X01 * (-ScalarType(63) * Y00 + ScalarType(28) * Y01 -
  //                         ScalarType(70) * Y02 + ScalarType(20) * Y12 -
  //                         ScalarType(35) * Y22) +
  //                  ScalarType(10) * X00 * Y22 + ScalarType(7) * X11 * Y22 -
  //                  ScalarType(63) * X12 * Y22 +
  //                  X02 * (ScalarType(21) * Y00 - ScalarType(70) * Y01 +
  //                         ScalarType(40) * Y02 + ScalarType(10) * Y11 -
  //                         ScalarType(70) * Y12 + ScalarType(21) * Y22)) -
  //            THREE *
  //                (BB *
  //                     (ScalarType(7) * Y00 * Y11 + ScalarType(10) * Y02 * Y11 -
  //                      ScalarType(35) * Y00 * Y12 +
  //                      ScalarType(7) * Y000 * (-THREE * Y1 + Y2) +
  //                      ScalarType(10) * Y00 * Y22 - ScalarType(35) * Y01 * Y22 +
  //                      ScalarType(7) * Y11 * Y22 + ScalarType(7) * Y0 * Y222 -
  //                      ScalarType(21) * Y1 * Y222) +
  //                 TWO * (ScalarType(5) * Z00 + ScalarType(10) * Z01p +
  //                        FOUR * Z02 - TWO * Z1p1p + ScalarType(10) * Z1p2 +
  //                        FIVE * Z22)),
  //        -SIX * AA *
  //                (ScalarType(46) * X02 * X11 - ScalarType(14) * X0 * X111 +
  //                 X1111 - ScalarType(14) * X111 * X2 -
  //                 ScalarType(14) * X01 * X22 + ScalarType(28) * X11 * X22 +
  //                 X00 * (ScalarType(28) * X11 - ScalarType(14) * X12 + X22)) -
  //            TWO * AB *
  //                (X22 * Y00 + ScalarType(112) * X01 * Y01 -
  //                 ScalarType(28) * X02 * Y01 - ScalarType(14) * X22 * Y01 -
  //                 ScalarType(28) * X01 * Y02 + FOUR * X02 * Y02 +
  //                 ScalarType(28) * X00 * Y11 - ScalarType(42) * X01 * Y11 +
  //                 ScalarType(46) * X02 * Y11 + ScalarType(28) * X22 * Y11 -
  //                 TWO * X12 *
  //                     (ScalarType(7) * Y00 - ScalarType(46) * Y01 +
  //                      ScalarType(14) * Y02 + ScalarType(21) * Y11 -
  //                      ScalarType(56) * Y12) -
  //                 ScalarType(14) * X00 * Y12 + ScalarType(92) * X01 * Y12 -
  //                 ScalarType(28) * X02 * Y12 + X00 * Y22 -
  //                 ScalarType(14) * X01 * Y22 +
  //                 X11 * (ScalarType(28) * Y00 - ScalarType(42) * Y01 +
  //                        ScalarType(46) * Y02 + SIX * Y11 -
  //                        ScalarType(42) * Y12 + ScalarType(28) * Y22)) +
  //            THREE *
  //                (-TWO * BB *
  //                     (ScalarType(46) * Y02 * Y11 - ScalarType(14) * Y0 * Y111 +
  //                      Y1111 - ScalarType(14) * Y111 * Y2 -
  //                      ScalarType(14) * Y01 * Y22 + ScalarType(28) * Y11 * Y22 +
  //                      Y00 * (ScalarType(28) * Y11 - ScalarType(14) * Y12 +
  //                             Y22)) +
  //                 FIVE * Z00 + ScalarType(40) * Z01p - TWO * Z02 +
  //                 ScalarType(8) * Z1p1p + ScalarType(40) * Z1p2 + FIVE * Z22),
  //        -TWO * AA *
  //                (ScalarType(3) * X02 * X11 - ScalarType(7) * X0 * X111 +
  //                 THREE * X1111 - ScalarType(7) * X111 * X2) -
  //            SIX * BB * Y02 * Y11 + ScalarType(14) * BB * Y0 * Y111 -
  //            SIX * BB * Y1111 -
  //            TWO * AB *
  //                (ScalarType(2) * X12 * Y01 - ScalarType(7) * X01 * Y11 +
  //                 X02 * Y11 - ScalarType(7) * X12 * Y11 +
  //                 X11 * (-ScalarType(7) * Y01 + Y02 + SIX * Y11 -
  //                        ScalarType(7) * Y12) +
  //                 TWO * X01 * Y12) +
  //            ScalarType(14) * BB * Y111 * Y2 - FIVE * Z01p + Z02 -
  //            ScalarType(7) * Z1p1p - FIVE * Z1p2,
  //        -(AA * X1111) - TWO * AB * X11 * Y11 - BB * Y1111 + Z1p1p});
  //   for (size_t i = 0; i < 3; ++i) {
  //     moments.volume() += ReturnScalarType(coeffs[i] * m0_basis[i]);
  //   }
  //   for (size_t i = 0; i < 4; ++i) {
  //     moments.centroid()[0] += ReturnScalarType(coeffs[3 + i] * m1x_basis[i]);
  //     moments.centroid()[1] += ReturnScalarType(coeffs[3 + i] * m1y_basis[i]);
  //   }
  //   for (size_t i = 0; i < 5; ++i) {
  //     moments.centroid()[2] += ReturnScalarType(coeffs[7 + i] * m1z_basis[i]);
  //   }
  //   moments.volume() *= ReturnScalarType(area_proj_triangle);
  //   moments.centroid()[0] *= ReturnScalarType(area_proj_triangle);
  //   moments.centroid()[1] *= ReturnScalarType(area_proj_triangle);
  //   moments.centroid()[2] *= ReturnScalarType(area_proj_triangle);
  //   return moments;
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    CylinderMomentArcIntegrator<ReturnType, ScalarType, 10> integrator(
        a_cylinder, a_arc, a_face_normal, 3);
    return integrator.integrate();
  } else {
    std::cout << "Type 3 for moments with order > 2 not yet implemented"
              << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

// template <class ReturnType, class ScalarType>
// ReturnType computeFaceOnlyContribution(
//     const AlignedParaboloidBase<ScalarType>& a_paraboloid,
//     const PlaneBase<ScalarType>& a_face_plane,
//     const PtBase<ScalarType>& a_pt_ref) {
//   using ReturnScalarType = typename ReturnType::value_type;
//   if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
//     /* Defining constants and types */
//     const ScalarType EPSILON = machine_epsilon<ScalarType>();
//     const ScalarType ZERO = ScalarType(0);
//     const ScalarType FOUR = ScalarType(4);

//     /* Function */
//     assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
//     assert(fabs(a_face_plane.normal()[2]) > EPSILON);
//     const ScalarType a =
//         -a_face_plane.normal()[0] / safelyEpsilon(a_face_plane.normal()[2]);
//     const ScalarType b =
//         -a_face_plane.normal()[1] / safelyEpsilon(a_face_plane.normal()[2]);
//     const ScalarType c =
//         a_face_plane.distance() / safelyEpsilon(a_face_plane.normal()[2]);
//     const ScalarType factor = FOUR * a_paraboloid.a() * a_paraboloid.b() * c -
//                               a_paraboloid.a() * b * b -
//                               a_paraboloid.b() * a * a;
//     return ReturnType::fromScalarConstant(ReturnScalarType(copysign(
//         machine_pi<ScalarType>() * factor * factor /
//             (ScalarType(32) *
//              pow(a_paraboloid.a() * a_paraboloid.b(), ScalarType(2.5))),
//         -a_face_plane.normal()[2])));
//   } else if constexpr (std::is_same_v<ReturnType,
//                                       VolumeMomentsBase<ReturnScalarType>>) {
//     /* Defining constants and types */
//     const ScalarType EPSILON = machine_epsilon<ScalarType>();
//     const ScalarType ZERO = ScalarType(0);
//     const ScalarType FOUR = ScalarType(4);
//     const ScalarType FIVE = ScalarType(5);

//     /* Function */
//     assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
//     assert(fabs(a_face_plane.normal()[2]) > EPSILON);
//     const ScalarType a = -a_face_plane.normal()[0] / a_face_plane.normal()[2];
//     const ScalarType b = -a_face_plane.normal()[1] / a_face_plane.normal()[2];
//     const ScalarType c = a_face_plane.distance() / a_face_plane.normal()[2];
//     const auto A = a_paraboloid.a(), B = a_paraboloid.b();
//     auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
//     const ScalarType factor = (a * a * B + A * (b * b - FOUR * B * c)) *
//                               (a * a * B + A * (b * b - FOUR * B * c)) *
//                               machine_pi<ScalarType>();
//     moments.volume() = ReturnScalarType(
//         copysign(factor / (ScalarType(32) * pow(A * B, ScalarType(2.5))),
//                  -a_face_plane.normal()[2]));
//     moments.centroid()[0] = ReturnScalarType(
//         a * B *
//         copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
//                  a_face_plane.normal()[2]));
//     moments.centroid()[1] = ReturnScalarType(
//         b * A *
//         copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
//                  a_face_plane.normal()[2]));
//     moments.centroid()[2] = ReturnScalarType(
//         (FIVE * A * (b * b) + FIVE * (a * a) * B - ScalarType(8) * A * B * c) *
//         copysign(factor / (ScalarType(384) * pow(A * B, ScalarType(3.5))),
//                  a_face_plane.normal()[2]));
//     return moments;
//   } else if constexpr (std::is_same_v<
//                            ReturnType,
//                            GeneralMomentsBase<2, 3, ReturnScalarType>>) {
//     /* Defining constants and types */
//     const ScalarType EPSILON = machine_epsilon<ScalarType>();
//     const ScalarType ZERO = ScalarType(0);
//     const ScalarType TWO = ScalarType(2);
//     const ScalarType FOUR = ScalarType(4);
//     const ScalarType FIVE = ScalarType(5);
//     const ScalarType SEVEN = ScalarType(7);

//     /* Function */
//     assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
//     assert(fabs(a_face_plane.normal()[2]) > EPSILON);
//     const ScalarType a = -a_face_plane.normal()[0] / a_face_plane.normal()[2];
//     const ScalarType b = -a_face_plane.normal()[1] / a_face_plane.normal()[2];
//     const ScalarType c = a_face_plane.distance() / a_face_plane.normal()[2];
//     const auto A = a_paraboloid.a(), B = a_paraboloid.b();
//     auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
//     const ScalarType factor = (a * a * B + A * (b * b - FOUR * B * c)) *
//                               (a * a * B + A * (b * b - FOUR * B * c)) *
//                               machine_pi<ScalarType>();
//     moments[0] = ReturnScalarType(
//         copysign(factor / (ScalarType(32) * pow(A * B, ScalarType(2.5))),
//                  -a_face_plane.normal()[2]));
//     moments[1] = ReturnScalarType(
//         a * B *
//         copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
//                  a_face_plane.normal()[2]));
//     moments[2] = ReturnScalarType(
//         b * A *
//         copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
//                  a_face_plane.normal()[2]));
//     moments[3] = ReturnScalarType(
//         (FIVE * A * (b * b) + FIVE * (a * a) * B - ScalarType(8) * A * B * c) *
//         copysign(factor / (ScalarType(384) * pow(A * B, ScalarType(3.5))),
//                  a_face_plane.normal()[2]));
//     moments[4] = ReturnScalarType(
//         (A * B * b * b + SEVEN * B * B * a * a - FOUR * A * B * B * c) *
//         copysign(factor / (ScalarType(768) * pow(A * B, ScalarType(4.5))),
//                  a_face_plane.normal()[2]));
//     moments[5] = ReturnScalarType(
//         (a * b) *
//         copysign(factor / (ScalarType(128) * pow(A * B, ScalarType(3.5))),
//                  a_face_plane.normal()[2]));
//     moments[6] = ReturnScalarType(
//         (B * B * a * a * a + A * B * a * b * b - TWO * A * B * B * a * c) *
//         copysign(factor / (ScalarType(128) * pow(A * B, ScalarType(4.5))),
//                  a_face_plane.normal()[2]));
//     moments[7] = ReturnScalarType(
//         (A * B * a * a + SEVEN * A * A * b * b - FOUR * A * A * B * c) *
//         copysign(factor / (ScalarType(768) * pow(A * B, ScalarType(4.5))),
//                  a_face_plane.normal()[2]));
//     moments[8] = ReturnScalarType(
//         (A * A * b * b * b + A * B * a * a * b - TWO * A * A * B * b * c) *
//         copysign(factor / (ScalarType(128) * pow(A * B, ScalarType(4.5))),
//                  a_face_plane.normal()[2]));
//     moments[9] = ReturnScalarType(
//         (SEVEN * A * A * b * b * b * b + TWO * SEVEN * A * B * a * a * b * b +
//          SEVEN * B * B * a * a * a * a -
//          ScalarType(24) * A * A * B * b * b * c -
//          ScalarType(24) * A * B * B * a * a * c +
//          ScalarType(16) * A * A * B * B * c * c) *
//         copysign(factor / (ScalarType(1024) * pow(A * B, ScalarType(4.5))),
//                  a_face_plane.normal()[2]));
//     return moments;
//   } else {
//     std::cout << "Type 4 for moments with order > 2 not yet implemented"
//               << std::endl;
//     return ReturnType::fromScalarConstant(ReturnScalarType(0));
//   }
// }

// template <class ReturnType, class ScalarType>
// ReturnType computeTriangleCorrection(
//     const AlignedParaboloidBase<ScalarType>& a_paraboloid,
//     const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1,
//     const PtBase<ScalarType>& a_pt_2) {
//   using ReturnScalarType = typename ReturnType::value_type;
//   if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
//     return ReturnType::fromScalarConstant(ReturnScalarType(
//         (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) * (a_pt_1[0] + a_pt_2[0]) +
//          -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) * (a_pt_1[1] + a_pt_2[1]) -
//          a_pt_0[2] - ScalarType(2) * a_pt_1[2] - a_pt_2[2]) /
//         ScalarType(12) *
//         ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
//          (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
//          (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0])));
//     return ReturnType::fromScalarConstant(ReturnScalarType(0));
//   } else if constexpr (std::is_same_v<ReturnType,
//                                       VolumeMomentsBase<ReturnScalarType>>) {
//     /* Defining constants and types */
//     const ScalarType ZERO = ScalarType(0);
//     const ScalarType ONE = ScalarType(1);
//     const ScalarType TWO = ScalarType(2);
//     const ScalarType THREE = ScalarType(3);
//     const ScalarType FOUR = ScalarType(4);
//     const ScalarType FIVE = ScalarType(5);
//     const ScalarType SIX = ScalarType(6);
//     const ScalarType HALF = ONE / TWO;

//     /* Function */
//     auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
//     const ScalarType A = a_paraboloid.a(), B = a_paraboloid.b();
//     const ScalarType X0 = a_pt_0[0], X1 = a_pt_1[0], X2 = a_pt_2[0];
//     const ScalarType Y0 = a_pt_0[1], Y1 = a_pt_1[1], Y2 = a_pt_2[1];
//     const ScalarType Z0 = a_pt_0[2], Z1 = a_pt_1[2], Z2 = a_pt_2[2];
//     const ScalarType triangle_area =
//         HALF * ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
//                 (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
//                 (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0]);
//     moments.volume() = ReturnScalarType(
//         (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) * (a_pt_1[0] + a_pt_2[0]) +
//          -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) * (a_pt_1[1] + a_pt_2[1]) -
//          a_pt_0[2] - TWO * a_pt_1[2] - a_pt_2[2]) *
//         triangle_area / SIX);
//     moments.centroid()[0] = ReturnScalarType(
//         triangle_area *
//         ((A * (-(X0 * X0 * X0) - X1 * X1 * X1 - X1 * X1 * X2 - X1 * (X2 * X2) -
//                X2 * X2 * X2 - X0 * X0 * (X1 + X2) -
//                X0 * (X1 * X1 + X1 * X2 + X2 * X2))) /
//              ScalarType(10) +
//          (B * (-(X1 * (Y0 * Y0 + TWO * Y0 * Y1 + THREE * (Y1 * Y1) + Y0 * Y2 +
//                        TWO * Y1 * Y2 + Y2 * Y2)) -
//                X2 * (Y0 * Y0 + Y0 * Y1 + Y1 * Y1 + TWO * Y0 * Y2 +
//                      TWO * Y1 * Y2 + THREE * (Y2 * Y2)) -
//                X0 * (THREE * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
//                      TWO * Y0 * (Y1 + Y2)))) /
//              ScalarType(30) +
//          (-(X0 * (TWO * Z0 + Z1 + Z2)) - X1 * (Z0 + TWO * Z1 + Z2) -
//           X2 * (Z0 + Z1 + TWO * Z2)) /
//              ScalarType(12)));
//     moments.centroid()[1] = ReturnScalarType(
//         -triangle_area *
//         ((B * (Y0 * Y0 * Y0 + Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) +
//                Y2 * Y2 * Y2 + Y0 * Y0 * (Y1 + Y2) +
//                Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2))) /
//              ScalarType(10) +
//          (A *
//           (X0 * X0 * (THREE * Y0 + Y1 + Y2) + X1 * X1 * (Y0 + THREE * Y1 + Y2) +
//            X2 * X2 * (Y0 + Y1 + THREE * Y2) + X1 * X2 * (Y0 + TWO * (Y1 + Y2)) +
//            X0 * (X1 * (TWO * Y0 + TWO * Y1 + Y2) +
//                  X2 * (TWO * Y0 + Y1 + TWO * Y2)))) /
//              ScalarType(30) +
//          (Y0 * (TWO * Z0 + Z1 + Z2) + Y1 * (Z0 + TWO * Z1 + Z2) +
//           Y2 * (Z0 + Z1 + TWO * Z2)) /
//              ScalarType(12)));
//     moments.centroid()[2] = ReturnScalarType(
//         triangle_area *
//         ((A * A *
//           (X0 * X0 * X0 * X0 + X1 * X1 * X1 * X1 + X1 * X1 * X1 * X2 +
//            X1 * X1 * (X2 * X2) + X1 * (X2 * X2 * X2) + X2 * X2 * X2 * X2 +
//            X0 * X0 * X0 * (X1 + X2) + X0 * X0 * (X1 * X1 + X1 * X2 + X2 * X2) +
//            X0 *
//                (X1 * X1 * X1 + X1 * X1 * X2 + X1 * (X2 * X2) + X2 * X2 * X2))) /
//              ScalarType(30) +
//          (B * B *
//           (Y0 * Y0 * Y0 * Y0 + Y1 * Y1 * Y1 * Y1 + Y1 * Y1 * Y1 * Y2 +
//            Y1 * Y1 * (Y2 * Y2) + Y1 * (Y2 * Y2 * Y2) + Y2 * Y2 * Y2 * Y2 +
//            Y0 * Y0 * Y0 * (Y1 + Y2) + Y0 * Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2) +
//            Y0 *
//                (Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) + Y2 * Y2 * Y2))) /
//              ScalarType(30) +
//          (A * B *
//           (X1 * X2 *
//                (Y0 * Y0 + THREE * (Y1 * Y1) + FOUR * Y1 * Y2 +
//                 THREE * (Y2 * Y2) + TWO * Y0 * (Y1 + Y2)) +
//            X0 * X0 *
//                (SIX * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
//                 THREE * Y0 * (Y1 + Y2)) +
//            X1 * X1 *
//                (Y0 * Y0 + SIX * (Y1 * Y1) + THREE * Y1 * Y2 + Y2 * Y2 +
//                 Y0 * (THREE * Y1 + Y2)) +
//            X2 * X2 *
//                (Y0 * Y0 + Y1 * Y1 + THREE * Y1 * Y2 + SIX * (Y2 * Y2) +
//                 Y0 * (Y1 + THREE * Y2)) +
//            X0 * (X1 * (THREE * (Y0 * Y0) + FOUR * Y0 * Y1 + THREE * (Y1 * Y1) +
//                        TWO * Y0 * Y2 + TWO * Y1 * Y2 + Y2 * Y2) +
//                  X2 * (THREE * (Y0 * Y0) + TWO * Y0 * Y1 + Y1 * Y1 +
//                        FOUR * Y0 * Y2 + TWO * Y1 * Y2 + THREE * (Y2 * Y2))))) /
//              ScalarType(90) +
//          (-(Z0 * Z0) - Z1 * Z1 - Z1 * Z2 - Z2 * Z2 - Z0 * (Z1 + Z2)) /
//              ScalarType(12)));
//     return moments;
//   } else if constexpr (std::is_same_v<
//                            ReturnType,
//                            GeneralMomentsBase<2, 3, ReturnScalarType>>) {
//     return ReturnType::fromScalarConstant(ReturnScalarType(0));
//   } else {
//     std::cout << "Type 5 for moments with order > 2 not yet implemented"
//               << std::endl;
//     return ReturnType::fromScalarConstant(ReturnScalarType(0));
//   }
// }

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_MOMENT_CONTRIBUTIONS_TPP_
