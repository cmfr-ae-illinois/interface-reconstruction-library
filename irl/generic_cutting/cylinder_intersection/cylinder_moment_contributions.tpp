// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
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
//#define DEBUG_CYL_IRL

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
        MomentCylinderIntegrand<ScalarType, 0, 0, 0>(pos_t, der_t, b, r) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand.centroid()[0] = ReturnScalarType(
        MomentCylinderIntegrand<ScalarType, 1, 0, 0>(pos_t, der_t, b, r) -
        MomentPlaneIntegrand<ScalarType, 1, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand.centroid()[1] = ReturnScalarType(
        MomentCylinderIntegrand<ScalarType, 0, 1, 0>(pos_t, der_t, b, r) -
        MomentPlaneIntegrand<ScalarType, 0, 1, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand.centroid()[2] = ReturnScalarType(
        MomentCylinderIntegrand<ScalarType, 0, 0, 1>(pos_t, der_t, b, r) -
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
template <class ReturnType, class ScalarType>
ReturnType computeType2Contribution(
    const AlignedCylinderBase<ScalarType>& a_aligned_cylinder,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1) {
  using ReturnScalarType = typename ReturnType::value_type;
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    const ScalarType ONESIXTH = ScalarType(1) / ScalarType(6);
    const ScalarType TWO = ScalarType(2);
    return ReturnType::fromScalarConstant(ReturnScalarType(
        ONESIXTH * (a_pt_0[1] - a_pt_1[1]) *
        (a_pt_0[0] * (TWO * a_pt_0[2] + a_pt_1[2]) +
         a_pt_1[0] * (TWO * a_pt_1[2] + a_pt_0[2]))));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType ONETWELVTH = ONE / ScalarType(12);
    const ScalarType ONE60TH = ONE / ScalarType(60);
    const ScalarType ONE180TH = ONE / ScalarType(180);

    const ScalarType DELTA_Y = (a_pt_0[1] - a_pt_1[1]);
    const ScalarType SX_SZ = (a_pt_0[0] + a_pt_1[0]) * 
                              (a_pt_0[2] + a_pt_1[2]);

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    moments.volume() = ReturnScalarType(
        DELTA_Y * (a_pt_0[0] * (TWO * a_pt_0[2] + a_pt_1[2]) +
                   a_pt_1[0] * (a_pt_0[2] + TWO * a_pt_1[2])) / ScalarType(6));
    moments.centroid()[0] = ReturnScalarType(
        DELTA_Y * (SX_SZ * (a_pt_0[0] + a_pt_1[0]) + 
                    TWO * a_pt_0[0] * a_pt_0[0] * a_pt_0[2] + 
                    TWO * a_pt_1[0] * a_pt_1[0] * a_pt_1[2]) / ScalarType(24));
    moments.centroid()[1] = ReturnScalarType(
        DELTA_Y * (SX_SZ * (a_pt_0[1] + a_pt_1[1]) + 
                    TWO * a_pt_0[0] * a_pt_0[1] * a_pt_0[2] + 
                    TWO * a_pt_1[0] * a_pt_1[1] * a_pt_1[2]) / ScalarType(12));
    moments.centroid()[2] = ReturnScalarType(
        DELTA_Y * (SX_SZ * (a_pt_0[2] + a_pt_1[2]) + 
                    TWO * a_pt_0[0] * a_pt_0[2] * a_pt_0[2] + 
                    TWO * a_pt_1[0] * a_pt_1[2] * a_pt_1[2]) / ScalarType(24));
    return moments;
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  } else {
    std::cout << "Type 2 for moments with order > 2 not yet implemented"
              << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 2> coeffsVC3SeriesOne(
    const ContainerType& a_weight) {
  std::array<ContainerType, 2> coeffs;
  coeffs.fill(ContainerType(ScalarType(0)));
  ContainerType x(1);
  UnsignedIndex_t i = 0;
  ScalarType max_diff;
  while (i <= 40) {
    max_diff = ScalarType(0);
    for (UnsignedIndex_t j = 0; j < 2; ++j) {
      ContainerType add_to_coeff = ScalarType(vc3Series[i][j]) * x;
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
inline std::array<ContainerType, 6> coeffsVC3andCC3SeriesOne(
    const ContainerType& a_weight) {
  std::array<ContainerType, 6> coeffs;
  coeffs.fill(ContainerType(ScalarType(0)));
  ContainerType x(ScalarType(1));
  UnsignedIndex_t i = 0;
  ScalarType max_diff;
  while (i <= 40) {
    max_diff = ScalarType(0);
    for (UnsignedIndex_t j = 0; j < 2; ++j) {
      ContainerType add_to_coeff = ScalarType(vc3Series[i][j]) * x;
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
      ContainerType add_to_coeff = ScalarType(cc3Series[i][j]) * x;
      coeffs[2 + j] += add_to_coeff;
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
inline std::array<ContainerType, 2> coeffsVC3Exact(
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
  return std::array<ContainerType, 2>(
      {(TWO * w4 - ScalarType(5) * w2 + 
         SIX * a_weight * T) *
           L2 / SIX,
       (TWO * w2 + w4 -
        SIX * w3 * T) *
           L2 / THREE});
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 6> coeffsVC3andCC3Exact(
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
  const auto L = ContainerType(ONE) / ((a_weight - ContainerType(ONE)) *
                                       (a_weight + ContainerType(ONE)));
  const auto L2 = L * L;
  const auto L3 = L * L * L;
  const auto S = (a_weight < ContainerType(ONE))
                     ? sqrt(ContainerType(ONE) - a_weight * a_weight)
                     : sqrt(a_weight * a_weight - ContainerType(ONE));
  const auto T = (a_weight < ContainerType(ONE))
                     ? atan((ContainerType(ONE) - a_weight) / S) / S
                     : atanh((a_weight - ContainerType(ONE)) / S) / S;
  return std::array<ContainerType, 6>(
      {(TWO * w4 - ScalarType(5) * w2 + 
         SIX * a_weight * T) *
           L2 / SIX,
       (TWO * w2 + w4 -
        SIX * w3 * T) *
           L2 / THREE,
       (ScalarType(23) * w2 
          - ScalarType(12) * w4
          + FOUR * w6
          - SIX * T * (THREE * a_weight
                       + TWO * w3)) *
            L3 / ScalarType(96),
       (- FIVE * w2 +
          TWO * w4 + SIX * T * a_weight) *
            L2 / ScalarType(48),
       (- ScalarType(8) * w2 - ScalarType(9) * w4 +
          TWO * w6 + ScalarType(30) * T * w3) *
            L3 / ScalarType(24),
       (ScalarType(13) * w4 + TWO * w6 -
          T * (SIX * w3 + ScalarType(24) * w5)) *
            L3 / ScalarType(24)});
}

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
    const ScalarType area_proj_triangle = 
        HALF * (pt_0[1] * (pt_1[2] - cp[2]) +
                pt_1[1] * (cp[2] - pt_0[2]) +
                cp[1] * (pt_0[2] - pt_1[2]));


    #ifdef DEBUG_CYL_IRL
      std::cout << "M3 computation :" << std::endl;
      std::cout << "x0 : " << pt_0 << std::endl;
      std::cout << "x* : " << cp << std::endl;
      std::cout << "x1 : " << pt_1 << std::endl;
      std::cout << "w : " << static_cast<double>(weight) << std::endl;
      std::cout << "area : " << static_cast<double>(area_proj_triangle) << std::endl;
    #endif

    assert(weight >= ZERO);
    std::array<ScalarType, 2> coeffs;
    if (weight < ScalarType(0.35)) { // We use the exact expressions
      #ifdef DEBUG_CYL_IRL
        std::cout << "weight is low, using exact value" << std::endl;
      #endif
      coeffs = coeffsVC3Exact<ScalarType, ScalarType>(weight);
    } else if (weight <
             ScalarType(1.7)) { // We use the 40th order Taylor series (w -> 1)
      #ifdef DEBUG_CYL_IRL
        std::cout << "weight is close to 1, using Taylor series" << std::endl;
      #endif
      coeffs = coeffsVC3SeriesOne<ScalarType, ScalarType>(weight);
    } else if (weight <
             ScalarType(1.0e9)) { // We use the series expansion (w -> infty)
      #ifdef DEBUG_CYL_IRL
        std::cout << "weight is high, using exact value" << std::endl;
      #endif
      coeffs = coeffsVC3Exact<ScalarType, ScalarType>(weight);
    } else { // This is within EPSILON of the actual value
      #ifdef DEBUG_CYL_IRL
        std::cout << "weight is BIG, using limite" << std::endl;
      #endif
      coeffs = std::array<ScalarType, 2>({ONE / THREE, ONE / THREE});
    }

    #ifdef DEBUG_CYL_IRL
      std::cout << "coeffs :" << std::endl;
      std::cout << "0 : " << static_cast<double>(coeffs[0]) << std::endl;
      std::cout << "1 : " << static_cast<double>(coeffs[1]) << std::endl;
    #endif
    const auto& C11 = (pt_0[0] + pt_1[0]);
    const auto& C12 = cp[0];
    #ifdef DEBUG_CYL_IRL
      std::cout << "C vector :" << std::endl;
      std::cout << "C11 : " << static_cast<double>(C11) << std::endl;
      std::cout << "C12 : " << static_cast<double>(C12) << std::endl;
    #endif
    return ReturnType::fromScalarConstant(ReturnScalarType(
        area_proj_triangle * (coeffs[0] * C11 +
         coeffs[1] * C12)));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
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
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    const auto& pt_0 = a_arc.start_point();
    const auto& cp = a_arc.control_point();
    const auto& pt_1 = a_arc.end_point();
    const auto weight = a_arc.weight();
    const auto X0 = pt_0[0], X1 = pt_1[0], X2 = cp[0];
    const auto Y0 = pt_0[1], Y1 = pt_1[1], Y2 = cp[1];
    const auto Z0 = pt_0[2], Z1 = pt_1[2], Z2 = cp[2];
    const ScalarType area_proj_triangle =
        HALF * (Y0 * (Z1 - Z2) + Y1 * (Z2 - Z0) + Y2 * (Z0 - Z1));
    #ifdef DEBUG_CYL_IRL
      std::cout << "M3 computation :" << std::endl;
      std::cout << "x0 : " << pt_0 << std::endl;
      std::cout << "x* : " << cp << std::endl;
      std::cout << "x1 : " << pt_1 << std::endl;
      std::cout << "w : " << static_cast<double>(weight) << std::endl;
      std::cout << "area : " << static_cast<double>(area_proj_triangle) << std::endl;
    #endif
    assert(weight >= ZERO);
    // Compute coefficients (functions of the weight)
    std::array<ScalarType, 6> coeffs;
    if (weight < ScalarType(0.35))  // We use the exact expressions
    {
      #ifdef DEBUG_CYL_IRL
        std::cout << "weight is low, using exact value" << std::endl;
      #endif
      coeffs = coeffsVC3andCC3Exact<ScalarType, ScalarType>(weight);
    } else if (weight <
               ScalarType(1.7))  // We use the 40th order Taylor series (w -> 1)
    {
      #ifdef DEBUG_CYL_IRL
        std::cout << "weight is close to 1, using Taylor series" << std::endl;
      #endif
      coeffs = coeffsVC3andCC3SeriesOne<ScalarType, ScalarType>(weight);
    } else if (weight < ScalarType(1.0e9))  // We use the exact expressions
    {
      #ifdef DEBUG_CYL_IRL
        std::cout << "weight is high, using exact value" << std::endl;
      #endif
      coeffs = coeffsVC3andCC3Exact<ScalarType, ScalarType>(weight);
    }
    // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
    //   coeffs = coeffsV3andC3SeriesInfinity(weight);
    else  // This is within EPSILON of the actual value
    {
      #ifdef DEBUG_CYL_IRL
        std::cout << "weight is BIG, using limite" << std::endl;
      #endif
      coeffs = std::array<ScalarType, 6>(
          {ScalarType(ONE / THREE), ScalarType(TWO / THREE),
           ScalarType(ONE / ScalarType(24)),
           ScalarType(ONE / ScalarType(24)),
           ScalarType(ONE / ScalarType(12)),
           ScalarType(ONE / ScalarType(12))});
    }
    #ifdef DEBUG_CYL_IRL
      std::cout << "coeffs :" << std::endl;
      std::cout << "0 : " << static_cast<double>(coeffs[0]) << std::endl;
      std::cout << "1 : " << static_cast<double>(coeffs[1]) << std::endl;
      std::cout << "2 : " << static_cast<double>(coeffs[2]) << std::endl;
      std::cout << "3 : " << static_cast<double>(coeffs[3]) << std::endl;
      std::cout << "4 : " << static_cast<double>(coeffs[4]) << std::endl;
      std::cout << "5 : " << static_cast<double>(coeffs[5]) << std::endl;
    #endif
    auto m0_basis = std::array<ScalarType, 2>(
        {X0 + X1,
         X2});
    auto m1x_basis = std::array<ScalarType, 4>(
        {(X0 + X1) * (X0 + X1),
         (X0 * X0 + X1 * X1),
         (X0 + X1) * X2,
         X2 * X2});
    auto m1y_basis = std::array<ScalarType, 4>(
        {TWO * (X0 + X1) * (Y0 + Y1),
         TWO * (X0 * Y0 + X1 * Y1),
         (X0 + X1) * Y2 + (Y0 + Y1) * X2,
         TWO * X2 * Y2});
    auto m1z_basis = std::array<ScalarType, 4>(
        {TWO * (X0 + X1) * (Z0 + Z1),
         TWO * (X0 * Z0 + X1 * Z1),
         (X0 + X1) * Z2 + (Z0 + Z1) * X2,
         TWO * X2 * Z2});
      #ifdef DEBUG_CYL_IRL
        std::cout << "C vector :" << std::endl;
        std::cout << "C11 : " << static_cast<double>(m0_basis[0]) << std::endl;
        std::cout << "C12 : " << static_cast<double>(m0_basis[1]) << std::endl;
        std::cout << "C21 : " << static_cast<double>(m1x_basis[0]) << std::endl;
        std::cout << "C22 : " << static_cast<double>(m1x_basis[1]) << std::endl;
        std::cout << "C23 : " << static_cast<double>(m1x_basis[2]) << std::endl;
        std::cout << "C24 : " << static_cast<double>(m1x_basis[3]) << std::endl;
        std::cout << "C31 : " << static_cast<double>(m1y_basis[0]) << std::endl;
        std::cout << "C32 : " << static_cast<double>(m1y_basis[1]) << std::endl;
        std::cout << "C33 : " << static_cast<double>(m1y_basis[2]) << std::endl;
        std::cout << "C34 : " << static_cast<double>(m1y_basis[3]) << std::endl;
        std::cout << "C41 : " << static_cast<double>(m1z_basis[0]) << std::endl;
        std::cout << "C42 : " << static_cast<double>(m1z_basis[1]) << std::endl;
        std::cout << "C43 : " << static_cast<double>(m1z_basis[2]) << std::endl;
        std::cout << "C44 : " << static_cast<double>(m1z_basis[3]) << std::endl;
      #endif
    for (size_t i = 0; i < 2; ++i) {
      moments.volume() += ReturnScalarType(coeffs[i] * m0_basis[i]);
    }
    for (size_t i = 0; i < 4; ++i) {
      moments.centroid()[0] += ReturnScalarType(coeffs[2 + i] * m1x_basis[i]);
      moments.centroid()[1] += ReturnScalarType(coeffs[2 + i] * m1y_basis[i]);
      moments.centroid()[2] += ReturnScalarType(coeffs[2 + i] * m1z_basis[i]);
    }
    moments.volume() *= ReturnScalarType(area_proj_triangle);
    moments.centroid()[0] *= ReturnScalarType(area_proj_triangle);
    moments.centroid()[1] *= ReturnScalarType(area_proj_triangle);
    moments.centroid()[2] *= ReturnScalarType(area_proj_triangle);
    return moments;
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

template <class ReturnType, class ScalarType>
ReturnType computeTriangleCorrection(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1,
    const PtBase<ScalarType>& a_pt_2) {
  using ReturnScalarType = typename ReturnType::value_type;
  /* Defining constants and types */
  const ScalarType ZERO = ScalarType(0);
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType THREE = ScalarType(3);
  const ScalarType FOUR = ScalarType(4);
  const ScalarType FIVE = ScalarType(5);
  const ScalarType SIX = ScalarType(6);
  const ScalarType HALF = ONE / TWO;
  const ScalarType SIXTH = ONE / SIX;
  const ScalarType X0 = a_pt_0[0], X1 = a_pt_1[0], X2 = a_pt_2[0];
  const ScalarType Y0 = a_pt_0[1], Y1 = a_pt_1[1], Y2 = a_pt_2[1];
  const ScalarType Z0 = a_pt_0[2], Z1 = a_pt_1[2], Z2 = a_pt_2[2];
  const ScalarType triangle_area =
      ((a_pt_0[0] - a_pt_2[0]) * (a_pt_1[1] - a_pt_2[1]) -
        (a_pt_1[0] - a_pt_2[0]) * (a_pt_0[1] - a_pt_2[1]));
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    return ReturnType::fromScalarConstant(ReturnScalarType(
        ((Y1 - Y0) * (X0 * (TWO * Z0 + Z1) + X1 * (Z0 + TWO * Z1)) +
         (Y0 - Y2) * (X2 * (TWO * Z2 + Z0) + X0 * (Z2 + TWO * Z0)) +
         (Y2 - Y1) * (X1 * (TWO * Z1 + Z2) + X2 * (Z1 + TWO * Z2)) -
         triangle_area * (Z0 + Z1 + Z2)) * SIXTH));
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    moments.volume() = ReturnScalarType(
        ((Y1 - Y0) * (X0 * (TWO * Z0 + Z1) + X1 * (Z0 + TWO * Z1)) +
         (Y0 - Y2) * (X2 * (TWO * Z2 + Z0) + X0 * (Z2 + TWO * Z0)) +
         (Y2 - Y1) * (X1 * (TWO * Z1 + Z2) + X2 * (Z1 + TWO * Z2)) -
         triangle_area * (Z0 + Z1 + Z2)) * SIXTH);
    moments.centroid()[0] = ReturnScalarType(
        ((Y1 - Y0) * ((X0 + X1) * (X0 + X1) * (Z0 + Z1) + 
                          TWO * X0 * X0 * Z0 + TWO * X1 * X1 * Z1) + 
         (Y0 - Y2) * ((X2 + X0) * (X2 + X0) * (Z2 + Z0) + 
                          TWO * X2 * X2 * Z2 + TWO * X0 * X0 * Z0) + 
         (Y2 - Y1) * ((X1 + X2) * (X1 + X2) * (Z1 + Z2) + 
                          TWO * X1 * X1 * Z1 + TWO * X2 * X2 * Z2) -
         triangle_area * (X0 * Z0 + X1 * Z1 + X2 * Z2 +
                          ((X0 + X1 + X2) * (Z0 + Z1 + Z2)))) / ScalarType(24));
    moments.centroid()[1] = ReturnScalarType(
        ((Y1 - Y0) * ((X0 + X1) * (Y0 + Y1) * (Z0 + Z1) + 
                          TWO * X0 * Y0 * Z0 + TWO * X1 * Y1 * Z1) + 
         (Y0 - Y2) * ((X2 + X0) * (Y2 + Y0) * (Z2 + Z0) + 
                          TWO * X2 * Y2 * Z2 + TWO * X0 * Y0 * Z0) + 
         (Y2 - Y1) * ((X1 + X2) * (Y1 + Y2) * (Z1 + Z2) + 
                          TWO * X1 * Y1 * Z1 + TWO * X2 * Y2 * Z2) -
         triangle_area * (Y0 * Z0 + Y1 * Z1 + Y2 * Z2 +
                          ((Y0 + Y1 + Y2) * (Z0 + Z1 + Z2))) * HALF) / ScalarType(12));
    moments.centroid()[2] = ReturnScalarType(
        ((Y1 - Y0) * ((X0 + X1) * (Z0 + Z1) * (Z0 + Z1) + 
                          TWO * X0 * Z0 * Z0 + TWO * X1 * Z1 * Z1) + 
         (Y0 - Y2) * ((X2 + X0) * (Z2 + Z0) * (Z2 + Z0) + 
                          TWO * X2 * Z2 * Z2 + TWO * X0 * Z0 * Z0) + 
         (Y2 - Y1) * ((X1 + X2) * (Z1 + Z2) * (Z1 + Z2) + 
                          TWO * X1 * Z1 * Z1 + TWO * X2 * Z2 * Z2) -
         triangle_area * (Z0 * (Z0 + Z1 + Z2) + Z1 * (Z1 + Z2) +
                          Z2 * Z2)) / ScalarType(24));
    return moments;
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  } else {
    std::cout << "Type 5 for moments with order > 2 not yet implemented"
              << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_MOMENT_CONTRIBUTIONS_TPP_
