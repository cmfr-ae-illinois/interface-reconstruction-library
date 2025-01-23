// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_INTERSECTION_TPP_
#define IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_INTERSECTION_TPP_

#include <float.h>
#include <cassert>
#include <cmath>
#include <random>
#include <type_traits>
#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_initializer.tpp"
#include "irl/generic_cutting/quadratic_intersection/moment_contributions.h"
#include "irl/generic_cutting/cylinder_intersection/cylinder_moment_contributions.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/geometry/half_edge_structures/brep_to_half_edge.h"
#include "irl/helpers/mymath.h"
#include "irl/moments/general_moments.h"
#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/quadratic_reconstruction/rational_bezier_arc.h"
#include "irl/generic_cutting/quadratic_intersection/quadratic_intersection.h"

#define NUMERICAL_INTEGRATION
// this enable a lot of debug text when computing
//#define VALDEBUG

namespace IRL {

/******************** Tangent to surface at given point ********************/
template <class ScalarType>
inline NormalBase<ScalarType> computeTangentVectorAtPoint(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_plane_normal,
    const PtBase<ScalarType>& a_pt) {
  // Defining constants
  const ScalarType EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ZERO = ScalarType(0);

  // Compute tangent
  NormalBase<ScalarType> surface_normal =
      getCylinderSurfaceNormal(a_cylinder, a_pt);
  surface_normal.approximatelyNormalize();
  NormalBase<ScalarType> tangent_at_pt =
      crossProduct(a_plane_normal, surface_normal);
  if (squaredMagnitude(tangent_at_pt) < EPSILON * EPSILON) {
    return NormalBase<ScalarType>(ZERO, ZERO, ZERO);
  }
  const ScalarType normal_correction = tangent_at_pt * a_plane_normal;
  tangent_at_pt = tangent_at_pt - normal_correction * a_plane_normal;
  return tangent_at_pt;
}

/************** Tangent to surface + orient at given point **************/
template <class ScalarType>
inline NormalBase<ScalarType> computeAndCorrectTangentVectorAtPt(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_plane_normal,
    const PtBase<ScalarType>& a_origin_pt, const PtBase<ScalarType>& a_end_pt,
    const NormalBase<ScalarType>& a_end_tangent,
    const PtBase<ScalarType>& a_intersection_pt) {
  // Defining constants
  const ScalarType ZERO = ScalarType(0);

  // Compute tangent
  NormalBase<ScalarType> tangent = computeTangentVectorAtPoint<ScalarType>(
      a_cylinder, a_plane_normal, a_intersection_pt);
  tangent.normalize();
  const NormalBase<ScalarType> edge_normal =
      crossProduct(a_plane_normal, a_end_pt - a_intersection_pt);
  if ((a_end_tangent * edge_normal > ZERO) == (tangent * edge_normal > ZERO)) {
    tangent = -tangent;
  }
  return tangent;
}

/***************** Calculate arc contribution (with split)
 * *******************/
template <class ReturnType, class ScalarType, class SurfaceOutputType,
          class PtType>
ReturnType computeType3ContributionWithSplit(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_plane_normal, const PtType& a_pt_ref,
    const PtType& a_pt_0, const PtType& a_pt_1,
    const NormalBase<ScalarType>& a_tangent_0,
    const NormalBase<ScalarType>& a_tangent_1, bool* a_requires_nudge,
    UnsignedIndex_t* a_split_counter, SurfaceOutputType* a_surface) {
  // Defining constants and types
  using ReturnScalarType = typename ReturnType::value_type;
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;
  using FloatType = float_type<ScalarType>;
  const ScalarType MACHINE_EPSILON = machine_epsilon<ScalarType>();
  const ScalarType DISTANCE_EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ANGLE_EPSILON = angle_epsilon<ScalarType>();
  const ScalarType ZERO = ScalarType(0);
  const ScalarType SPLIT_TOL1 = ScalarType(0.9);
  const ScalarType SPLIT_TOL2 = ScalarType(0.999);
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType FOUR = ScalarType(4);
  const ScalarType HALF = ONE / TWO;
  const ScalarType ONEQUARTER = ONE / FOUR;
  const ScalarType THREEQUARTERS = HALF + ONEQUARTER;
  const ScalarType ONE_HUNDRED = ScalarType(100);

    #ifdef VALDEBUG
    std::cout << "trying to compute Type 3 Contri, counter : " << *a_split_counter << std::endl;
    #endif

  // Store reference point, start point and end point of arc
  const Pt& pt_ref = a_pt_ref.getPt();
  const Pt& pt_0 = a_pt_0.getPt();
  const Pt& pt_1 = a_pt_1.getPt();
    #ifdef VALDEBUG
    std::cout << "going from " << pt_0 << std::endl;
    std::cout << "to " << pt_1 << std::endl;
    std::cout << "with tangents " << a_tangent_0 << std::endl;
    std::cout << "and " << a_tangent_1 << std::endl;
    #endif

  // Calculate edge vector and its normalized version
  const Normal edge_vector = pt_1 - pt_0;
  Normal edge_vector_normalized = edge_vector;
  edge_vector_normalized.normalize();

  // Compute dot product between normalized edge and end-point tangents
  const ScalarType tgt0_dot_edge = a_tangent_0 * edge_vector_normalized;
  const ScalarType tgt1_dot_edge = a_tangent_1 * edge_vector_normalized;

  // If start and end point are very close and tangents point toward each
  // other: the arc has no contribution to the moments
  if (squaredMagnitude(edge_vector) < DISTANCE_EPSILON * DISTANCE_EPSILON &&
      fabs(ONE - tgt0_dot_edge) < ANGLE_EPSILON &&
      fabs(ONE + tgt1_dot_edge) < ANGLE_EPSILON) {
    // For completeness, we update the parametric surface boundary with a
    // straight Bezier arc
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      auto surface_arc = RationalBezierArc(
          pt_0.toDoublePt(), 0.5 * (pt_0.toDoublePt() + pt_1.toDoublePt()),
          pt_1.toDoublePt(), 0.0);
      // TODO: check that this recast works in DP
      surface_arc.reset_start_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
      surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_1));
      a_surface->addArc(surface_arc);
    }
    return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
  }

  // We split the Bezier arc if:
  // - the tangents both point toward the same infinite point (to avoid
  // control points very far away)
  // - the start tangent points in the opposite direction as the edge
  // - the end tangent points in the same direction as the edge
  // This is strictly speaking a bit of an overkill, but increases accuracy
  // and stability
  bool split = ((a_tangent_0 * a_tangent_1) > SPLIT_TOL1 ||
                tgt0_dot_edge < ZERO || tgt1_dot_edge > ZERO);

  // Slip the arc and compute the contribution of the children arcs
  if (split) {
    #ifdef VALDEBUG
    std::cout << "we are splitting" << std::endl;
    #endif
    (*a_split_counter)++;
    // Compute average point and tangent
    const Pt average_pt = HALF * (pt_0 + pt_1);
    auto average_tangent = Normal(HALF * (a_tangent_0 + a_tangent_1));
    // If the norm of the average tangent is very small (meaning that the
    // tangents are aligned), then we switch to QP and shake the polytope
    if (squaredMagnitude(average_tangent) < DISTANCE_EPSILON) {
      *a_requires_nudge = true;
      return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    }
    // Let's make sure the average tangent belongs to the plane of the face
    const ScalarType normal_correction = average_tangent * a_plane_normal;
    average_tangent = average_tangent - normal_correction * a_plane_normal;
    average_tangent.approximatelyNormalize();
    // Find the intersection between the half-line starting from the middle
    // of the arc and pointing in the direction of the average tangent. This
    // will be the end-point and start-point of the two new arcs
    Pt projected_pt = projectPtAlongHalfLineOntoCylinder<ScalarType>(
        a_cylinder, average_tangent, average_pt);
    // If this point could not be found, switch to QP and shake the polytope
    if (projected_pt[0] == ScalarType(DBL_MAX)) {
      *a_requires_nudge = true;
      return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    }
    // If we have had to split more than N times, something is wrong: switch
    // to QP and shake the polytope
    if constexpr (std::is_same_v<FloatType, double>) {
      if (*a_split_counter > 5) {
        *a_requires_nudge = true;
        return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
      }
    } else {
      if (*a_split_counter > 10) {
        *a_requires_nudge = true;
        return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
      }
    }
    // Compute the tangent at the new projected point
    Normal tangent_projected_pt =
        computeAndCorrectTangentVectorAtPt<ScalarType>(
            a_cylinder, a_plane_normal, pt_0, pt_1, a_tangent_1,
            projected_pt);

    // If the projected tangent cannot be calculated, then the arc contribution is 0
    if (tangent_projected_pt[0] == ZERO && tangent_projected_pt[1] == ZERO &&
         tangent_projected_pt[2] == ZERO) {
      if constexpr (!std::is_same_v<SurfaceOutputType, NoSurfaceOutput>) {
        auto surface_arc = RationalBezierArc(
            pt_0.toDoublePt(), 0.5 * (pt_0.toDoublePt() + pt_1.toDoublePt()),
            pt_1.toDoublePt(), 0.0);
        surface_arc.reset_start_point_id(
            reinterpret_cast<std::uintptr_t>(&pt_0));
        surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_1));
        a_surface->addArc(surface_arc);
      }
      return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    }

    // If we want to output the surface:
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      // We need to store this vertex so that its address remains
      // unique over time (for surface output purposes)
      Pt* new_point = new Pt(projected_pt);
      PtBase<double>* new_point_double =
          new PtBase<double>(static_cast<double>(projected_pt[0]),
                             static_cast<double>(projected_pt[1]),
                             static_cast<double>(projected_pt[2]));
      a_surface->addPt(new_point_double);
      return computeType3ContributionWithSplit<ReturnType, ScalarType>(
                 a_cylinder, a_plane_normal, a_pt_ref, a_pt_0, *new_point,
                 a_tangent_0, tangent_projected_pt, a_requires_nudge,
                 a_split_counter, a_surface) +
             computeType3ContributionWithSplit<ReturnType, ScalarType>(
                 a_cylinder, a_plane_normal, a_pt_ref, *new_point, a_pt_1,
                 -tangent_projected_pt, a_tangent_1, a_requires_nudge,
                 a_split_counter, a_surface);
    } else {
      return computeType3ContributionWithSplit<ReturnType, ScalarType>(
                 a_cylinder, a_plane_normal, a_pt_ref, a_pt_0,
                 PtType(projected_pt), a_tangent_0, tangent_projected_pt,
                 a_requires_nudge, a_split_counter, a_surface) +
             computeType3ContributionWithSplit<ReturnType, ScalarType>(
                 a_cylinder, a_plane_normal, a_pt_ref, PtType(projected_pt),
                 a_pt_1, -tangent_projected_pt, a_tangent_1, a_requires_nudge,
                 a_split_counter, a_surface);
    }
  }
  // We do not split and compute the moment contributions of the arc
  else {
    const auto arc = RationalBezierArcBase<ScalarType>(
        pt_0, a_tangent_0, pt_1, a_tangent_1, a_plane_normal, a_cylinder);
    // Add the arc to the surface output
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      auto surface_arc = RationalBezierArc(
          pt_0.toDoublePt(), a_tangent_0.toDoubleNormal(), pt_1.toDoublePt(),
          a_tangent_1.toDoubleNormal(), a_plane_normal.toDoubleNormal(),
          AlignedCylinderBase<double>(a_cylinder));
      surface_arc.reset_start_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
      surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_1));
      a_surface->addArc(surface_arc);
    }
    // If the rational Bezier weight is negative, we switch to QP and shake
    // the polytope
    if (arc.weight() < ZERO) {
    #ifdef VALDEBUG
    std::cout << "the arc is computed with negative weight, nudging" << std::endl;
    #endif
      *a_requires_nudge = true;
      return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    }
    #ifdef VALDEBUG
    std::cout << "ok, let's continue" << std::endl;
    #endif
    // Calculate type 3 contribution of arc
    auto moments = computeType3Contribution<ReturnType, ScalarType>(
        a_cylinder, arc, a_plane_normal);
    // // If the arc was split, then we need to add the contribution of the
    // // space between the splitted arcs and the orignial arc
    // if (!(&a_pt_ref == &a_pt_0 || &a_pt_ref == &a_pt_1)) {
    //   moments += computeTriangleCorrection<ReturnType, ScalarType>(
    //       a_cylinder, pt_0, pt_1, pt_ref);
    // }
    return moments;
  }
}

/**************** Calculate all arc contributions ******************/
template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType, class PtType, class NormalType>
ReturnType computeNewEdgeSegmentContribution(
    const AlignedCylinderBase<ScalarType>& a_aligned_cylinder,
    const PtType& a_ref_pt, const HalfEdgeType a_entry_half_edge,
    const HalfEdgeType a_exit_half_edge, bool* skip_first,
    const NormalType& a_face_normal, const UnsignedIndex_t a_proj_dir,
    const bool a_ignore_type3, bool* a_requires_nudge,
    SurfaceOutputType* a_surface) {
  using ReturnScalarType = typename ReturnType::value_type;
  ReturnType full_moments = ReturnType::fromScalarConstant(ReturnScalarType(0));
  // Handle new edge on exit->entry
  full_moments += computeType1Contribution<ReturnType, ScalarType>(
      a_ref_pt, a_exit_half_edge->getVertex()->getLocation(),
      a_entry_half_edge->getVertex()->getLocation(), skip_first, true,
      a_face_normal, a_proj_dir);
  full_moments += computeType2Contribution<ReturnType, ScalarType>(
      a_aligned_cylinder, a_exit_half_edge->getVertex()->getLocation(),
      a_entry_half_edge->getVertex()->getLocation());
  if (!a_ignore_type3) {
    full_moments += orientAndApplyType3Correction<ReturnType, ScalarType>(
        a_aligned_cylinder, a_exit_half_edge, a_entry_half_edge,
        a_requires_nudge, a_surface);
  }
  return full_moments;
}

/************* Calculate moments from non-aligned paraboloid **********/
template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class CylinderType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithCylinder(SegmentedHalfEdgePolyhedronType* a_polytope,
                                  HalfEdgePolytopeType* a_complete_polytope,
                                  const CylinderType& a_cylinder) {
  // Defining type aliases (needed to ensure precision is consistent)
  // Definining scalar container: This can be a double/__float128, or a
  // scalar with embedded derivatives
  using ScalarType = typename CylinderType::value_type;
  using FloatType = float_type<ScalarType>;
  using PtType = typename SegmentedHalfEdgePolyhedronType::pt_type;
  static_assert(std::is_same_v<typename PtType::value_type, ScalarType>);
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using ReferenceFrame = ReferenceFrameBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;

  // Defining constants
  const ScalarType DISTANCE_EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ZERO = ScalarType(0);
  const ScalarType ONE = ScalarType(1);
  const ScalarType THREE = ScalarType(3);

  // Shortcuts if we already now that the paraboloid is entirely above or
  // below the polytope
  ReturnType moments;
  #ifdef VALDEBUG
  std::cout << "going to compute the volume of this cylinder : " << a_cylinder << std::endl;
  #endif
  if (a_cylinder.isAlwaysAbove()) {
    #ifdef VALDEBUG
    std::cout << "The cylinder is always above, premature ending" << std::endl;
    #endif
    if constexpr (has_cylinder_surface<ReturnType>::value) {
      moments.getMoments() =
          ReturnType::moment_type::calculateMoments(a_polytope);
      moments.getSurface().setCylinder(a_cylinder);
    } else {
      moments = ReturnType::calculateMoments(a_polytope);
    }
    return moments;
  } else if (a_cylinder.isAlwaysBelow()) {
    #ifdef VALDEBUG
    std::cout << "The cylinder is always below, premature ending" << std::endl;
    #endif
    if constexpr (has_cylinder_surface<ReturnType>::value) {
      moments.getMoments() = ReturnType::moment_type::fromScalarConstant(ZERO);
      moments.getSurface().setCylinder(a_cylinder);
    } else {
      moments = ReturnType::fromScalarConstant(ZERO);
    }
    return moments;
  }

  auto cylinder = a_cylinder;

  // Move into reference frame of the cylinder and compute and approximate
  // length-scale of the polytope
  const UnsignedIndex_t original_number_of_vertices =
      a_polytope->getNumberOfVertices();
  const auto& datum = cylinder.getDatum();
  const auto& ref_frame = cylinder.getReferenceFrame();

  const Pt start_pt = a_polytope->getVertex(0)->getLocation().getPt() - datum;
  ScalarType max_dist_sq = ZERO;

  #ifdef VALDEBUG
  std::cout << "The original number of vertices is : " << original_number_of_vertices << std::endl;
  #endif
  for (UnsignedIndex_t v = 0; v < original_number_of_vertices; ++v) {
    const Pt original_pt =
        a_polytope->getVertex(v)->getLocation().getPt() - datum;
    if (v > 0) {
      max_dist_sq =
          maximum(max_dist_sq, squaredMagnitude(original_pt - start_pt));
    }
    PtType projected_location;
    auto& pt = projected_location.getPt();
    for (UnsignedIndex_t n = 0; n < 3; ++n) {
      pt[n] = ref_frame[n] * original_pt;
    }
    a_polytope->getVertex(v)->setLocation(projected_location);
  }

  // Define scale so that the polyhedron's volume is O(1)
  const ScalarType inv_scale =
      maximum(ScalarType(1.0e6) * DISTANCE_EPSILON, sqrt(max_dist_sq));
  // const ScalarType inv_scale = ScalarType(ONE);
  const ScalarType inv_volume_scale = inv_scale * inv_scale * inv_scale;
  const ScalarType scale = ScalarType(ONE) / inv_scale;
  const ScalarType volume_scale = scale * scale * scale;

  // Normalized polyhedron
  for (UnsignedIndex_t v = 0; v < original_number_of_vertices; ++v) {
    auto& pt = a_polytope->getVertex(v)->getLocation().getPt();
    pt *= scale;
  }

  // Normalized cylinder
  auto scaled_aligned_cylinder = AlignedCylinder(std::array<ScalarType, 2>{
      cylinder.getAlignedCylinder().b(),
      cylinder.getAlignedCylinder().r() * scale * scale});

  // Compute moments of intersection
  if constexpr (has_cylinder_surface<ReturnType>::value) {
    moments.getSurface().setCylinder(cylinder);
    moments.getMoments() = intersectPolyhedronWithAlignedCylinder<
        typename ReturnType::moment_type>(
        a_polytope, a_complete_polytope, scaled_aligned_cylinder,
        inv_volume_scale, &moments.getSurface());
  } else {
    NoSurfaceOutput* surf = nullptr;
    moments = intersectPolyhedronWithAlignedCylinder<ReturnType>(
        a_polytope, a_complete_polytope, scaled_aligned_cylinder,
        inv_volume_scale, surf);
  }

  // Un-normalized moments
  if constexpr (has_cylinder_surface<ReturnType>::value) {
    if constexpr (!is_moments_volume<typename ReturnType::moment_type>::value) {
      moments.getMoments().centroid().getPt() *= inv_volume_scale * inv_scale;
      moments.getMoments().volume() *= inv_volume_scale;
    } else {
      moments.getMoments().volume() *= inv_volume_scale;
    }
    auto& arc_list = moments.getSurface().getArcs();
    for (std::size_t i = 0; i < arc_list.size(); ++i) {
      arc_list[i].start_point() *= static_cast<double>(inv_scale);
      arc_list[i].control_point() *= static_cast<double>(inv_scale);
      arc_list[i].end_point() *= static_cast<double>(inv_scale);
    }
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ScalarType>>) {
    moments.centroid().getPt() *= inv_volume_scale * inv_scale;
    moments.volume() *= inv_volume_scale;
  } else if constexpr (std::is_same_v<ReturnType,
                                      GeneralMomentsBase<2, 3, ScalarType>>) {
    moments[0] *= inv_volume_scale;
    moments[1] *= inv_volume_scale * inv_scale;
    moments[2] *= inv_volume_scale * inv_scale;
    moments[3] *= inv_volume_scale * inv_scale;
    moments[4] *= inv_volume_scale * inv_scale * inv_scale;
    moments[5] *= inv_volume_scale * inv_scale * inv_scale;
    moments[6] *= inv_volume_scale * inv_scale * inv_scale;
    moments[7] *= inv_volume_scale * inv_scale * inv_scale;
    moments[8] *= inv_volume_scale * inv_scale * inv_scale;
    moments[9] *= inv_volume_scale * inv_scale * inv_scale;
  } else {
    moments.volume() *= inv_volume_scale;
  }

  // Move first moment back to original frame of reference
  if constexpr (has_cylinder_surface<ReturnType>::value) {
    if constexpr (!is_moments_volume<typename ReturnType::moment_type>::value) {
      auto pt = Pt(ZERO, ZERO, ZERO);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          pt[n] += ref_frame[d][n] * moments.getMoments().centroid().getPt()[d];
        }
      }
      pt += moments.getMoments().volume() * datum;
      moments.getMoments().centroid().getPt() = pt;
    }
  } else {
    if constexpr (std::is_same_v<ReturnType, VolumeMomentsBase<ScalarType>>) {
      auto pt = Pt(ZERO, ZERO, ZERO);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          pt[n] += ref_frame[d][n] * moments.centroid().getPt()[d];
        }
      }
      pt += moments.volume() * datum;
      moments.centroid().getPt() = pt;
    } else if constexpr (std::is_same_v<ReturnType,
                                        GeneralMomentsBase<2, 3, ScalarType>>) {
      const Eigen::Matrix<ScalarType, 3, 1> D{datum[0], datum[1], datum[2]};
      const Eigen::Matrix<ScalarType, 3, 3> R{
          {ref_frame[0][0], ref_frame[1][0], ref_frame[2][0]},
          {ref_frame[0][1], ref_frame[1][1], ref_frame[2][1]},
          {ref_frame[0][2], ref_frame[1][2], ref_frame[2][2]}};
      const ScalarType M0 = moments[0];
      const Eigen::Matrix<ScalarType, 3, 1> M1prime{moments[1], moments[2],
                                                    moments[3]};
      const Eigen::Matrix<ScalarType, 3, 3> M2prime{
          {moments[4], moments[5], moments[6]},
          {moments[5], moments[7], moments[8]},
          {moments[6], moments[8], moments[9]}};
      const Eigen::Matrix<ScalarType, 3, 1> M1 = R * M1prime + M0 * D;
      const Eigen::Matrix<ScalarType, 3, 3> M2 =
          R * M2prime * R.transpose() + R * (M1prime * D.transpose()) +
          (D * M1prime.transpose()) * R.transpose() + M0 * (D * D.transpose());
      moments[1] = M1(0);
      moments[2] = M1(1);
      moments[3] = M1(2);
      moments[4] = M2(0, 0);
      moments[5] = M2(0, 1);
      moments[6] = M2(0, 2);
      moments[7] = M2(1, 1);
      moments[8] = M2(1, 2);
      moments[9] = M2(2, 2);
    }
  }

  return moments;
}

/*********** Calculate moments from aligned paraboloid *************/
template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AlignedCylinderType,
          class ScalarType, class SurfaceOutputType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithAlignedCylinder(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AlignedCylinderType& a_cylinder,
    const ScalarType a_inv_volume_scale, SurfaceOutputType* a_surface) {
  // Below function computes the entire integration (nudge counter
  // initialized to 0)
  return formCylinderIntersectionBases<ReturnType>(
      a_polytope, a_complete_polytope, a_cylinder, 0, a_surface);
}

/******************* Find intersections on segment **********************/
template <class PtType, class ScalarType>
inline void checkAndFindIntercepts(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1,
    StackVector<PtBase<ScalarType>, 2>* a_intercepts,
    const ScalarType a_nudge_epsilon, const bool a_elliptic) {
  static_assert(std::is_same_v<PtType, PtBase<ScalarType>>);
  const ScalarType EPSILON_LO = -ScalarType(0.5) * a_nudge_epsilon;
  const ScalarType EPSILON_HI = ScalarType(1) - EPSILON_LO;
  const ScalarType ZERO = ScalarType(0);
  a_intercepts->resize(0);
  // Compute coefficients of quadratic equation
  const auto& pt_0 = a_pt_0.getPt();
  const auto& pt_1 = a_pt_1.getPt();
  const auto pt_diff = pt_1 - pt_0;
  const ScalarType a = pt_diff[2] * pt_diff[2] +
                       a_cylinder.b() * pt_diff[1] * pt_diff[1];
  const ScalarType b =
      ScalarType(2) * (pt_diff[2] * pt_0[2] +
                       a_cylinder.b() * pt_diff[1] * pt_0[1]);  
  const ScalarType c = pt_0[2] * pt_0[2] +
                       a_cylinder.b() * pt_0[1] * pt_0[1] - a_cylinder.r();
  // Solve quadratic equation
  const StackVector<ScalarType, 2> solutions =
      solveQuadratic<ScalarType>(a, b, c);
  // Convert solution into point
  for (auto& solution : solutions) {
    auto pt_z = pt_0[2] + solution * pt_diff[2];
    if (solution > EPSILON_LO && solution < EPSILON_HI && pt_z >= ZERO) {
      a_intercepts->push_back(PtBase<ScalarType>(pt_0 + solution * pt_diff));
    }
  }
}

/**************** Flag: is vertex below aligned Cylinder?
 * *****************/
template <class VertexType>
bool vertexBelow(const VertexType& a_pt,
                 const AlignedCylinderBase<typename VertexType::value_type>&
                     a_cylinder) {
  using ScalarType = typename VertexType::value_type;
  const auto& pt = a_pt.getPt();
  return a_cylinder.r() - pt[2] * pt[2] - a_cylinder.b() * pt[1] * pt[1] > ScalarType(0);
}

// If centroid is outside of polygon, or any vertex on face
// is inside the ellipse, then the ellipse is not contained by the face.
template <class ScalarType, class HalfEdgeType>
bool ellipseContainedInFace(
    const AlignedCylinderBase<ScalarType>& a_aligned_cylinder,
    const PlaneBase<ScalarType>& a_face_plane, HalfEdgeType* const a_half_edge,
    bool* a_requires_nudge) {
  /* Defining constants and types */
  const ScalarType ZERO = ScalarType(0);
  const ScalarType TWO = ScalarType(2);
  const std::array<ScalarType, 2> ZERO_ZERO{{ZERO, ZERO}};

  /* Function */
  const auto& face_normal = a_face_plane.normal();

  // if nx is too small, the intersection is two parralelle lines
  // it can't be contain in the face
  if (fabs(face_normal[0]) < distance_epsilon<ScalarType>()) {
    return false;
  }

  const auto& face_distance = a_face_plane.distance();
  const std::array<ScalarType, 3> conic_center{
      {face_distance / (face_normal[0]),
       ZERO,
       ZERO}};
  // const ScalarType delta_face = face_distance / face_normal[2];
  // const ScalarType gamma_face =
  //     a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
  //     a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] - delta_face;
  // if (fabs(gamma_face) < distance_epsilon<ScalarType>()) {
  //   *a_requires_nudge = true;
  //   return false;
  // }
  // if (a_aligned_paraboloid.a() * gamma_face < ZERO) {
  //   return false;
  // }

  // Not true for a cylinder, an ellipse and ce in the yz plane

  // // First we will check if centroid is in the bounding box
  // // of the face polygon. Due to the fact we are assuming
  // // the paraboloid was in Z direction, we assume the ellipse
  // // lives on the x/y plane and project the polygon
  // // down to it as well (essentially neglecting the z component).
  auto current_half_edge = a_half_edge;
  std::array<ScalarType, 3> xyz_min{
    {ScalarType(DBL_MAX), ScalarType(DBL_MAX), ScalarType(DBL_MAX)}};
  std::array<ScalarType, 3> xyz_max{
      {-ScalarType(DBL_MAX), -ScalarType(DBL_MAX), -ScalarType(DBL_MAX)}};
  do {
    const PtBase<ScalarType>& location =
        current_half_edge->getVertex()->getLocation().getPt();
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      xyz_min[d] = minimum(xyz_min[d], location[d]);
      xyz_max[d] = maximum(xyz_max[d], location[d]);
    }
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_half_edge);
  if (conic_center[0] < xyz_min[0] || conic_center[0] > xyz_max[0] ||
      conic_center[1] < xyz_min[1] || conic_center[1] > xyz_max[1] ||
      conic_center[2] < xyz_min[2] || conic_center[2] > xyz_max[2]) {
    return false;
  }

  // If not outside bounding box, need to use
  // a more refined test based on the number
  // of edges a ray emitted from the point crosses.
  // Even edges means the point is outside the face.
  current_half_edge = a_half_edge;
  bool pt_internal_to_polygon = false;
  do {
    const PtBase<ScalarType>& location_0 =
        current_half_edge->getPreviousVertex()->getLocation().getPt();
    const PtBase<ScalarType>& location_1 =
        current_half_edge->getVertex()->getLocation().getPt();
    if (isPtBeforeIntersectionWithEdgeYZ<ScalarType>(ZERO_ZERO, location_0,
                                                   location_1)) {
      pt_internal_to_polygon = !pt_internal_to_polygon;
    }
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_half_edge);
  return pt_internal_to_polygon;
}

template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType>
ReturnType orientAndApplyType3Correction(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    HalfEdgeType* a_start, HalfEdgeType* a_end, bool* a_requires_nudge,
    SurfaceOutputType* a_surface) {
  // Defining constants and types
  using ReturnScalarType = typename ReturnType::value_type;
  using FloatType = float_type<ScalarType>;
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;
  const ScalarType MACHINE_EPSILON = machine_epsilon<ScalarType>();
  const ScalarType DISTANCE_EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ANGLE_EPSILON = angle_epsilon<ScalarType>();
  const ScalarType ZERO = ScalarType(0);
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType HALF = ONE / TWO;
  const ScalarType TEN = ScalarType(10);
  const ScalarType ONE_HUNDRED = ScalarType(100);
  
    #ifdef VALDEBUG
    std::cout << "Start computing M3" << std::endl;
    #endif

  // Store start and end points, end-start vector, and plane properties
  const auto& pt_0 = a_start->getVertex()->getLocation().getPt();
  const auto& pt_1 = a_end->getVertex()->getLocation().getPt();
  const auto edge_vector = Normal(pt_1 - pt_0);
  const auto& face_plane = a_end->getFace()->getPlane();
  const auto& face_normal = face_plane.normal();

    #ifdef VALDEBUG
    std::cout << "pt0 " << pt_0 << std::endl;
    std::cout << "pt1 " << pt_1 << std::endl;
    std::cout << "edge_vector " << edge_vector << std::endl;
    std::cout << "face_plane " << face_plane << std::endl;
    std::cout << "face_normal " << face_normal << std::endl;
    #endif

  // Compute tangents at start and end points. THEY ARE NOT YET NORMALIZED!
  Normal tgt_0 =
      computeTangentVectorAtPoint<ScalarType>(a_cylinder, face_normal, pt_0);
  Normal tgt_1 =
      computeTangentVectorAtPoint<ScalarType>(a_cylinder, face_normal, pt_1);
  
    #ifdef VALDEBUG
    std::cout << "tgt_0 " << tgt_0 << std::endl;
    std::cout << "tgt_1 " << tgt_1 << std::endl;
    #endif

  // Is the arc from an ellipse (hence could require splittin)?
  const bool elliptic_face = fabs(face_normal[0]) > MACHINE_EPSILON;

  // If the tangents could not be calculated, switch to QP and shake the
  // polytope
  if ((tgt_0[0] == ZERO && tgt_0[1] == ZERO && tgt_0[2] == ZERO) ||
      (tgt_1[0] == ZERO && tgt_1[1] == ZERO && tgt_1[2] == ZERO)) {
    *a_requires_nudge = true;
    return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
  }

  // If the start and end points almost coincide, switch to QP
  if constexpr (std::is_same_v<ScalarType, double>) {
    if (squaredMagnitude(edge_vector) < DISTANCE_EPSILON * DISTANCE_EPSILON) {
      *a_requires_nudge = true;
      return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    }
  }

  // CASE: The arc is a straight line
  if (!elliptic_face) {
    #ifdef VALDEBUG
    std::cout << "this is a stright line" << std::endl;
    #endif
    auto control_pt = Pt(HALF* pt_0 + HALF * pt_1);
    const auto arc = RationalBezierArcBase<ScalarType>(
        pt_1, control_pt, pt_0, HALF);
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      auto surface_arc =
          RationalBezierArc(pt_1.toDoublePt(), control_pt.toDoublePt(),
                            pt_0.toDoublePt(), double(1)/double(2));
      surface_arc.reset_start_point_id(
          reinterpret_cast<std::uintptr_t>(&pt_1));
      surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
      a_surface->addArc(surface_arc);
    }
    return computeType3Contribution<ReturnType, ScalarType>(a_cylinder, arc,
                                                            face_normal);

  }  // CASE: The arc is from an ellipse or hyperbola
  else {
    #ifdef VALDEBUG
    std::cout << "this is an ellipse" << std::endl;
    #endif
    // We need to normalize the tangents because we will compute their
    // angular orientation
    tgt_0.normalize();
    tgt_1.normalize();

    // Compute edge vectors and check if tangent is parallel to
    // edge
    Normal edge_0 = a_start->getVertex()->getLocation().getPt() -
                    a_start->getPreviousVertex()->getLocation().getPt();
    Normal edge_1 = a_end->getVertex()->getLocation().getPt() -
                    a_end->getPreviousVertex()->getLocation().getPt();
    edge_0.normalize();
    edge_1.normalize();
    bool tgt_0_parallel_edge_0 =
        fabs(ONE - fabs(tgt_0 * edge_0)) < ANGLE_EPSILON;
    bool tgt_1_parallel_edge_1 =
        fabs(ONE - fabs(tgt_1 * edge_1)) < ANGLE_EPSILON;

    // If both tangents are not parallel to the polytope edge from which
    // they originate, we can use the face normal + edge to orient the
    // tangent towards the inside of the face
    if (!tgt_0_parallel_edge_0 && !tgt_1_parallel_edge_1) {
      const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
      const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
      tgt_0 = (edge_normal_0 * tgt_0 < ZERO) ? -tgt_0 : tgt_0;
      tgt_1 = (edge_normal_1 * tgt_1 < ZERO) ? -tgt_1 : tgt_1;
    }
    // Otherwise, things get tricky...
    else {
      // If we are in DP, we directly switch to QP
      if constexpr (std::is_same_v<ScalarType, double>) {
        *a_requires_nudge = true;
        return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
      }
      // Compute the center of the ellipse
      const Pt conic_center = conicCenter<ScalarType>(face_plane, a_cylinder);
      // If the center coincides with one of the end points, the moment
      // contribution is 0
      if (squaredMagnitude(pt_0 - conic_center) <
              DISTANCE_EPSILON * DISTANCE_EPSILON ||
          squaredMagnitude(pt_1 - conic_center) <
              DISTANCE_EPSILON * DISTANCE_EPSILON) {
        return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
      }
      // If the start/end points are both the origin and the plane of the
      // face contains the origin, the contribution is 0
      if (squaredMagnitude(pt_0) < DISTANCE_EPSILON * DISTANCE_EPSILON &&
          squaredMagnitude(pt_1) < DISTANCE_EPSILON * DISTANCE_EPSILON &&
          face_plane.distance() < DISTANCE_EPSILON) {
        return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
      }
      // Use ellipse center to orient the tangents so as to intersect. They
      // may still both point in the wrong direction
      auto center_to_pt_0 = Normal(pt_0 - conic_center);
      auto center_to_pt_1 = Normal(pt_1 - conic_center);
      center_to_pt_0.normalize();
      center_to_pt_1.normalize();
      Normal dummy_tgt_0 = crossProduct(face_normal, center_to_pt_0);
      Normal dummy_tgt_1 = crossProduct(face_normal, center_to_pt_1);
      assert(fabs(tgt_0 * dummy_tgt_0) > ANGLE_EPSILON);
      assert(fabs(tgt_1 * dummy_tgt_1) > ANGLE_EPSILON);
      tgt_0 = (tgt_0 * dummy_tgt_0) < ZERO ? -tgt_0 : tgt_0;
      tgt_1 = (tgt_1 * dummy_tgt_1) > ZERO ? -tgt_1 : tgt_1;
      // At this point, the tangents form a valid arc (but they
      // may be oriented in the wrong direction)
      if (!tgt_0_parallel_edge_0) {
        const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
        if (edge_normal_0 * tgt_0 < ZERO) {
          tgt_0 = -tgt_0;
          tgt_1 = -tgt_1;
        }
      } else if (!tgt_1_parallel_edge_1) {
        const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
        if (edge_normal_1 * tgt_1 < ZERO) {
          tgt_0 = -tgt_0;
          tgt_1 = -tgt_1;
        }
      } else {
        *a_requires_nudge = true;
        return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
        // if (a_paraboloid.a() < ZERO) {
        //   tgt_0 = -tgt_0;
        //   tgt_1 = -tgt_1;
        // }
      }
    }
    UnsignedIndex_t split_counter = 0;
    return computeType3ContributionWithSplit<ReturnType, ScalarType>(
        a_cylinder, face_normal, pt_1, pt_1, pt_0, tgt_1, tgt_0,
        a_requires_nudge, &split_counter, a_surface);
  }

  return ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
}

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AligneCylinderType,
          class SurfaceOutputType>
ReturnType reformQuadraticIntersectionBases(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AligneCylinderType& a_aligned_cylinder,
    const UnsignedIndex_t a_nudge_iter, SurfaceOutputType* a_surface) {
  // Find out scalar type of polytope (DP or QP)
  using vertex_type = typename SegmentedHalfEdgePolyhedronType::vertex_type;
  using pt_type = typename vertex_type::pt_type;
  using ScalarType = typename pt_type::value_type;
  using FloatType = float_type<ScalarType>;

  if constexpr (std::is_same_v<FloatType, double>) {
    // This is needed to convert from DP to QP
    using QP_scalar_type = convert_to_quad<ScalarType>;
    using QP_pt_type = PtBase<QP_scalar_type>;
    using QP_vertex_type = VertexQuadratic<QP_pt_type>;
    using QP_halfedge_type = HalfEdgeQuadratic<QP_vertex_type>;
    using QP_face_type = FaceQuadratic<QP_halfedge_type>;
    const UnsignedIndex_t QP_kMaxHalfEdges = HalfEdgePolytopeType::maxHalfEdges;
    const UnsignedIndex_t QP_kMaxVertices = HalfEdgePolytopeType::maxVertices;
    const UnsignedIndex_t QP_kMaxFaces = HalfEdgePolytopeType::maxFaces;
    using QP_complete_polytope_type = HalfEdgePolyhedronQuadratic<
        QP_pt_type, QP_vertex_type, QP_halfedge_type, QP_face_type,
        QP_kMaxHalfEdges, QP_kMaxVertices, QP_kMaxFaces>;

    // Convert aligned cylinder to QP
    const auto QP_aligned_cylinder =
        AlignedCylinderBase<QP_scalar_type>(a_aligned_cylinder);

    // Convert polytope to QP
    QP_complete_polytope_type QP_polytope_cylinder;
    convertPolytopeFromDoubleToQuadPrecision(a_polytope, a_complete_polytope,
                                             &QP_polytope_cylinder);
    auto QP_segmented_cylinder =
        QP_polytope_cylinder.generateSegmentedPolyhedron();

    if (!QP_segmented_cylinder.checkValidHalfEdgeStructure()) {
      std::cout << "Polytope is not valid after conversion to QP!" << std::endl;
      std::cout << "PolytopeDP:" << std::endl;
      std::cout << *a_polytope << std::endl;
      std::cout << "PolytopeQP:" << std::endl;
      std::cout << QP_segmented_cylinder << std::endl;
      exit(-1);
    }

    assert(QP_segmented_cylinder.checkValidHalfEdgeStructure());

    // Nudge polytope and reset surface
    nudgePolyhedron(&QP_segmented_cylinder, &QP_polytope_cylinder,
                    a_nudge_iter, a_surface);
    // Try again!
    return formCylinderIntersectionBases<ReturnType>(
        &QP_segmented_cylinder, &QP_polytope_cylinder,
        QP_aligned_cylinder, a_nudge_iter + 1, a_surface);
  } else {
    // Nudge polytope (already QP) and reset surface
    nudgePolyhedron(a_polytope, a_complete_polytope, a_nudge_iter, a_surface);
    // Try again!
    return formCylinderIntersectionBases<ReturnType>(
        a_polytope, a_complete_polytope, a_aligned_cylinder, a_nudge_iter + 1,
        a_surface);
  }
}

// Assumes cylinder of function 0 = a*x^2 + b*y^2 + z.
// We will truncate the polyhedron to exist in the region
// q < a*x^2 + b*y^2 + z
template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AligneCylinderType,
          class SurfaceOutputType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
formCylinderIntersectionBases(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AligneCylinderType& a_aligned_cylinder,
    const UnsignedIndex_t a_nudge_iter, SurfaceOutputType* a_surface) {
  using vertex_type = typename SegmentedHalfEdgePolyhedronType::vertex_type;
  using pt_type = typename vertex_type::pt_type;
  using face_type = typename SegmentedHalfEdgePolyhedronType::face_type;
  using half_edge_type = typename HalfEdgePolytopeType::half_edge_type;

  // Defining type aliases (needed to ensure precision is consistent)
  using ScalarType = typename pt_type::value_type;
  using FloatType = float_type<ScalarType>;
  using ReturnScalarType = typename ReturnType::value_type;
  static_assert(std::is_same_v<typename pt_type::value_type, ScalarType>);
  static_assert(
      std::is_same_v<typename half_edge_type::value_type, ScalarType>);
  static_assert(std::is_same_v<typename face_type::value_type, ScalarType>);
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedCylinder = AlignedCylinderBase<ScalarType>;

  // Defining constants
  const ScalarType MACHINE_EPSILON = machine_epsilon<ScalarType>();
  const ScalarType DISTANCE_EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ANGLE_EPSILON = angle_epsilon<ScalarType>();
  const ScalarType ZERO = ScalarType(0);
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType THREE = ScalarType(3);
  const ScalarType FIVE = ScalarType(5);
  const ScalarType TEN = ScalarType(10);
  const ScalarType ONE_HUNDRED = ScalarType(100);
  
  #ifdef VALDEBUG
  std::cout << "========= try to compute moment ==========" << std::endl;
  std::cout << "========= nb iter : " << a_nudge_iter << " ==========" << std::endl;
  #endif

  assert(!(a_surface != nullptr &&
           std::is_same<SurfaceOutputType, NoSurfaceOutput>::value));

  // Initialize moments to 0
  ReturnType full_moments = ReturnType::fromScalarConstant(ReturnScalarType(0));

  // Initialising variables for handling degenerate cases
  bool requires_nudge = false;
  const ScalarType nudge_epsilon = DISTANCE_EPSILON;
  const ScalarType nudge_epsilon_sq = nudge_epsilon * nudge_epsilon;

  if constexpr (std::is_same_v<FloatType, Quad_t>) {
    // We only check this in QP for performance purposes
    // (i.e. when a_nudge_iter > 0)
    if (a_nudge_iter >= 100) {
      std::cout << "ERROR: Nudged more than 100 times. Moments returned "
                   "are wrong -> Context: b = "
                << a_aligned_cylinder.b()
                << ", r = " << a_aligned_cylinder.r() << std::endl;
      std::ofstream myfile("failed_nudge_comparison_cell.vtu");
      if (myfile.is_open()) {
        myfile << *a_polytope;
        myfile.close();
      }
      return ReturnType::fromScalarConstant(-ReturnScalarType(DBL_MAX));
    }
  }

  const auto& b = a_aligned_cylinder.b();
  const auto& r = a_aligned_cylinder.r();
  // Identify elliptic case
  const bool elliptic = b > ZERO;

  // First, triangulate faces (if necessary) and compute normals
  // The triangulation criterion is based on face planarity
  triangulatePolytopeAndComputeNormals(a_polytope, a_complete_polytope,
                                       nudge_epsilon, &requires_nudge);

  // Mark vertices clipped(= above paraboloid) or unclipped(= below
  // paraboloid)
  const auto starting_number_of_vertices = a_polytope->getNumberOfVertices();
  const auto starting_number_of_faces = a_polytope->getNumberOfFaces();
  UnsignedIndex_t number_of_vertices_above = 0;

  #ifdef VALDEBUG
  std::cout << "checking wich vertices is clipped" << std::endl;
  #endif
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    vertex.setAsUnnecessaryToSeek();  // Reset all
    vertex.markToBeClipped();

    #ifdef VALDEBUG
    std::cout << "point nb : " << v << ", vertex : " << vertex << std::endl;
    #endif

    const auto& pt = vertex.getLocation().getPt();
    const auto& hdist = sqrt(r/b) - fabs(pt[1]);
    const ScalarType dist_function = pt[2] * pt[2] + b * pt[1] * pt[1] - r;
    #ifdef VALDEBUG
    std::cout << "computed distance is " << dist_function << std::endl;
    #endif

    if (fabs(dist_function) < nudge_epsilon) {
      // If a polytope vertex lies within nudge_eps of the paraboloid
      // we directly require a nudge and switch to QP
  #ifdef VALDEBUG
  std::cout << "too close, gona nudge" << std::endl;
  #endif
      requires_nudge = true;
      break;
    } else if (dist_function < ZERO) {
  #ifdef VALDEBUG
  std::cout << "it's not clipped" << std::endl;
  #endif
      vertex.markToBeNotClipped();
    } else {
  #ifdef VALDEBUG
  std::cout << "it's clipped" << std::endl;
  #endif
      vertex.markToBeClipped();
      ++number_of_vertices_above;
    }
  }

  #ifdef VALDEBUG
  std::cout << "number_of_vertices above : " << number_of_vertices_above << std::endl;
  #endif
  if (!requires_nudge) {
    // Early termination cases, only possible with elliptic thanks to
    // convexity
    if (elliptic && number_of_vertices_above == 0) {
      // Whole volume below
      return ReturnType::calculateMoments(a_polytope);
    }

    if (!elliptic && number_of_vertices_above == starting_number_of_vertices) {
      // Zero volume - will be current value of full_moments
      return full_moments;
    }
  } else {
    // Nudge and try again!
    return reformQuadraticIntersectionBases<ReturnType>(
        a_polytope, a_complete_polytope, a_aligned_cylinder, a_nudge_iter,
        a_surface);
  }

  // Clear visitation knowledge from polytope.
  for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
    (*a_polytope)[f]->markAsNotVisited();
    (*a_polytope)[f]->clearIntersections();
  }

  // Will have 0, 1, or 2 intercepts.
  // If intercept exists, place into HalfEdgeStructure
  const bool check_from_clipped = true;
  const bool check_from_unclipped = true;
  #ifdef VALDEBUG
  std::cout << "\n========computing intersections==========" << std::endl;
  #endif

  //////// Compute intersections between edges and cylinder
  // Temporary stack vector to store single/double interesects
  StackVector<pt_type, 2> edge_intercepts;
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    vertex.setToSeek();
    if (check_from_clipped && vertex.isClipped()) {
      // CASE WHERE STARTING VERTEX IS CLIPPED
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        // If it has needsToSeek set, it means it is a new vertex
        // or already visited. Either way, do not need to check
        // for intersection
        const auto& vertex_start = current_edge->getPreviousVertex();
        if (vertex_start->needsToSeek()) {
          current_edge =
              current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
          continue;
        }
        const auto& vertex_end = current_edge->getVertex();
        if (vertex_start->isNotClipped() && vertex_end->isNotClipped()) {
          current_edge =
              current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
          continue;
        }

        // If previous vertex is not clipped, single-intercept
        // If previous vertex is clipped, need to check for
        // double-intercept Checking for double-intercept and
        // calculating single-intercept is the same routine, so
        // just always do it.
        const auto& edge_start = vertex_start->getLocation();
        const auto& edge_end = vertex_end->getLocation();
        checkAndFindIntercepts<pt_type, ScalarType>(
            a_aligned_cylinder, edge_start, edge_end, &edge_intercepts,
            nudge_epsilon, elliptic);

        // Size of returned intercepts indicates single or double
        // intercept (or none)
        if (edge_intercepts.size() == 1) {

  #ifdef VALDEBUG
  std::cout << "one intersection find at " << edge_intercepts[0] << std::endl;
  #endif
          // Check for intersection near end points
          if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[0] - edge_end) <
                  nudge_epsilon_sq) {
            requires_nudge = true;
  #ifdef VALDEBUG
  std::cout << "gonna nudge" << std::endl;
  #endif
            break;
          }

          // Add vertex to half edge structure, resetting
          // connectivity and creating a new half edge (and new
          // opposite half edge)
          placeSingleIntercept(edge_intercepts[0], current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->addIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(opposite_half_edge);
          opposite_face->addIntersection();
        } else if (edge_intercepts.size() == 2) {
  #ifdef VALDEBUG
  std::cout << "two intersection find at " << edge_intercepts[0] << " and " << edge_intercepts[1] << std::endl;
  #endif
          // Check for intersection near end points
          if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[1] - edge_end) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[0] - edge_intercepts[1]) <
                  nudge_epsilon_sq) {
            requires_nudge = true;
  #ifdef VALDEBUG
  std::cout << "gonna nudge" << std::endl;
  #endif
            break;
          }

          // Add the two vertices to half edge structure, resetting
          // connectivity and creating two new half edges (and new
          // opposite half edges)
          placeDoubleIntercept(edge_intercepts, current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge()->getPreviousHalfEdge());
          current_face->addDoubleIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(
              current_edge->getOppositeHalfEdge());
          opposite_face->addDoubleIntersection();
        }
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    } else if (check_from_unclipped && vertex.isNotClipped()) {
      // CASE WHERE STARTING VERTEX IS UNCLIPPED
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        // If it has needsToSeek set, it means it is a new vertex
        // or already visited. Either way, do not need to check
        // for intersection
        const auto& vertex_start = current_edge->getPreviousVertex();
        if (vertex_start->needsToSeek()) {
          current_edge =
              current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
          continue;
        }
        const auto& vertex_end = current_edge->getVertex();
        if (elliptic) {
          if (a_aligned_cylinder.b() > ZERO) {
            if (vertex_start->isNotClipped() && vertex_end->isNotClipped()) {
              current_edge =
                  current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
              continue;
            }
          } else if (vertex_start->isClipped() && vertex_end->isClipped()) {
            current_edge =
                current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
            continue;
          }
        }

        // If previous vertex is clipped, single-intercept
        // If previous vertex is not clipped, need to check for
        // double-intercept Checking for double-intercept and
        // calculating single-intercept is the same routine, so
        // just always do it.
        const auto& edge_start = vertex_start->getLocation();
        const auto& edge_end = vertex_end->getLocation();
        checkAndFindIntercepts<pt_type, ScalarType>(
            a_aligned_cylinder, edge_start, edge_end, &edge_intercepts,
            nudge_epsilon, elliptic);

        // Size of returned intercepts indicates single or double
        // intercept (or none)
        if (edge_intercepts.size() == 1) {
  #ifdef VALDEBUG
  std::cout << "one intersection find at " << edge_intercepts[0] << std::endl;
  #endif
          // Check for intersection near end points
          if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[0] - edge_end) <
                  nudge_epsilon_sq) {
            requires_nudge = true;
  #ifdef VALDEBUG
  std::cout << "gonna nudge" << std::endl;
  #endif
            break;
          }

          // Add vertex to half edge structure, resetting
          // connectivity and creating a new half edge (and new
          // opposite half edge)
          placeSingleIntercept(edge_intercepts[0], current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge());
          current_face->addIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->addIntersection();
        } else if (edge_intercepts.size() == 2) {
  #ifdef VALDEBUG
  std::cout << "two intersection find at " << edge_intercepts[0] << " and " << edge_intercepts[1] << std::endl;
  #endif
          // Check for intersection near end point
          if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[1] - edge_end) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[0] - edge_intercepts[1]) <
                  nudge_epsilon_sq) {
            requires_nudge = true;
  #ifdef VALDEBUG
  std::cout << "gonna nudge" << std::endl;
  #endif
            break;
          }

          // Add the two vertices to half edge structure, resetting
          // connectivity and creating two new half edges (and new
          // opposite half edges)
          placeDoubleIntercept(edge_intercepts, current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge());
          current_face->addDoubleIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(
              opposite_half_edge->getNextHalfEdge());
          opposite_face->addDoubleIntersection();
        }
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    }
    // If nudge is requested, we exit loop and mark unvisited vertices
    // as visited
    if (requires_nudge) {
      for (UnsignedIndex_t i = v + 1; i < starting_number_of_vertices; ++i) {
        a_polytope->getVertex(i)->setToSeek();
      }
      break;
    }
  }
  assert(a_polytope->checkValidHalfEdgeStructure());
  // After leaving above loop, all faces that were intersected
  // have a starting half edge that ends on a new (intersection)
  // vertex and the previous vertex is unclipped (or new). (i.e.,
  // the edge exists in the unclipped portion. All new vertices
  // also have their half edge being one that has its edge in the
  // unclipped portion.

  const auto vertices_after_intersection = a_polytope->getNumberOfVertices();
  const auto new_intersection_vertices =
      vertices_after_intersection - starting_number_of_vertices;

  #ifdef VALDEBUG
  std::cout << "\n========Done intersection=======" << std::endl;
  std::cout << "vertices after intersection : " << vertices_after_intersection << std::endl;
  std::cout << "nb new vertices : " << new_intersection_vertices << std::endl;
  #endif
  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    // Original vertices will be set as needsToSeek()
    // Reset newly created vertices will have doesNotNeedToSeek()
    a_polytope->getVertex(v)->setAsUnnecessaryToSeek();
  }

  // If intersection is too close to corners, switch to QP,
  // shift and rotate polyhedron randomly, and try again
  if (requires_nudge) {
    // Clean half-edge structure by removing intersections
    assert(a_polytope->checkValidHalfEdgeStructure());
    if (new_intersection_vertices > 0) {
      resetPolyhedron(a_polytope, a_complete_polytope);
      assert(a_polytope->getNumberOfVertices() == starting_number_of_vertices);
    }
    assert(a_polytope->checkValidHalfEdgeStructure());

    // Nudge and try again!
    return reformQuadraticIntersectionBases<ReturnType>(
        a_polytope, a_complete_polytope, a_aligned_cylinder, a_nudge_iter,
        a_surface);
  }

  // impossible for a cylinder
  // // Check for face-only intersections (that are: the paraboloid intersects
  // // the face but not the edges). Can only happen with elliptic paraboloids
  // if (elliptic) {
  //   if (a_aligned_paraboloid.a() > ZERO) {
  //     for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
  //       auto& face = *(*a_polytope)[f];
  //       if (face.getNumberOfIntersections() > 0) {
  //         continue;
  //       }
  //       // Face will be FaceOnly intersect if
  //       // 1 - Any vertex on the face is below the paraboloid
  //       // 2 - No z-component on the face is below the maximum for
  //       // a downward opening elliptic paraboloid, which is 0.
  //       bool face_valid = false;
  //       const auto starting_half_edge = face.getStartingHalfEdge();
  //       auto current_half_edge = starting_half_edge;
  //       do {
  //         const auto& vertex = *(current_half_edge->getVertex());
  //         if (vertex.isNotClipped()) {
  //           face_valid = false;
  //           break;
  //         }
  //         if (vertex.getLocation().getPt()[2] < ZERO) {
  //           face_valid = true;
  //         }
  //         current_half_edge = current_half_edge->getNextHalfEdge();
  //       } while (current_half_edge != starting_half_edge);
  //       if (!face_valid) {
  //         continue;
  //       }

  //       // Made it this far, check if there is a face-only
  //       // intersection
  //       const auto& face_plane = face.getPlane();
  //       if (fabs(face_plane.normal()[2]) > MACHINE_EPSILON) {
  //         // Get ellipse on this face
  //         if (ellipseContainedInFace<ScalarType>(a_aligned_paraboloid,
  //                                                face_plane, starting_half_edge,
  //                                                &requires_nudge)) {
  //           full_moments += computeFaceOnlyContribution<ReturnType, ScalarType>(
  //               a_aligned_paraboloid, face_plane,
  //               starting_half_edge->getVertex()->getLocation());
  //           // Return surface parametrization
  //           if constexpr (!std::is_same<SurfaceOutputType,
  //                                       NoSurfaceOutput>::value) {
  //             addEllipseToSurfaceOutput<ScalarType>(a_aligned_paraboloid,
  //                                                   face_plane, a_surface);
  //           }
  //         }
  //       }
  //     }
  //   } else {
  //     // Case where a_aligned_paraboloid.a() < 0.0)
  //     for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
  //       auto& face = *(*a_polytope)[f];
  //       if (face.getNumberOfIntersections() > 0) {
  //         continue;
  //       }
  //       // Face will be FaceOnly intersect if
  //       // 1 - Any vertex on the face is above the paraboloid
  //       // 2 - No z-component on the face is above the minimum for
  //       // an upward opening elliptic paraboloid, which is 0.
  //       bool face_valid = false;
  //       const auto starting_half_edge = face.getStartingHalfEdge();
  //       auto current_half_edge = starting_half_edge;
  //       do {
  //         const auto& vertex = *(current_half_edge->getVertex());
  //         if (vertex.isClipped()) {
  //           face_valid = false;
  //           break;
  //         }
  //         if (vertex.getLocation().getPt()[2] > ZERO) {
  //           face_valid = true;
  //         }
  //         current_half_edge = current_half_edge->getNextHalfEdge();
  //       } while (current_half_edge != starting_half_edge);
  //       if (!face_valid) {
  //         continue;
  //       }

  //       // Made it this far, check if there is a face-only
  //       // intersection
  //       const auto& face_plane = face.getPlane();
  //       if (fabs(face_plane.normal()[2]) > MACHINE_EPSILON) {
  //         // Get ellipse on this face
  //         if (ellipseContainedInFace<ScalarType>(a_aligned_paraboloid,
  //                                                face_plane, starting_half_edge,
  //                                                &requires_nudge)) {
  //           full_moments += computeFaceOnlyContribution<ReturnType, ScalarType>(
  //               a_aligned_paraboloid, face_plane,
  //               starting_half_edge->getVertex()->getLocation());
  //           if constexpr (!std::is_same<SurfaceOutputType,
  //                                       NoSurfaceOutput>::value) {
  //             addEllipseToSurfaceOutput<ScalarType>(a_aligned_paraboloid,
  //                                                   face_plane, a_surface);
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  // If no edges are intersected by the cylinder, our job is done
  if (new_intersection_vertices == 0) {

    if (number_of_vertices_above == starting_number_of_vertices) {
  #ifdef VALDEBUG
  std::cout << "All vertices above, fast end" << std::endl;
  #endif
      // All points above
      return full_moments;
    } else if (number_of_vertices_above == 0) {
  #ifdef VALDEBUG
  std::cout << "All vertices below, fast end" << std::endl;
  #endif
      // All points below
      return ReturnType::calculateMoments(a_polytope) + full_moments;
    }
    // ELSE: The polytope has 0 intersections, but a number of clipped
    // vertices > 0 and a number of unclipped vertices > 0. This is very
    // rare, and happens when the polytope actually consists of several
    // manifolds that are either entirely above or entirely below the
    // paraboloid. In which case, the remainder of this function will
    // compute the volume of the unclipped elements by looping over their
    // faces. else {
  }

  #ifdef VALDEBUG
  std::cout << "\n=======Starting computing the faces=========" << std::endl;
  #endif
  // There are edge-intersections. We then need to calculate the
  // contributions of all unclipped faces to the moments
  for (UnsignedIndex_t f = 0; f < a_polytope->getNumberOfFaces(); ++f) {
  #ifdef VALDEBUG
  std::cout << "\ndoing face number : " << f << std::endl;
  #endif
    auto& face = *(*a_polytope)[f];
  #ifdef VALDEBUG
  std::cout << "doing face : " << face << std::endl;
  #endif
    const auto& face_normal = face.getPlane().normal();
    // If magnitude(normal) ~ 0, the face area is ~ 0 so we can skip
    if (squaredMagnitude(face_normal) < MACHINE_EPSILON) {
  #ifdef VALDEBUG
  std::cout << "face is small, skip to the next one" << std::endl;
  #endif
      continue;
    }

    // The starting half-edge is, by construction, an intersection which is
    // an entry (i.e. goes from above to below paraboloid, following the
    // half-edge structure)
    auto starting_half_edge = face.getStartingHalfEdge();

    // Count intersections, this cannot be an odd number otherwise it
    // disobeys the Jordan theorem
    const auto intersection_size = face.getNumberOfIntersections();
    if (intersection_size % 2 == 1) {
      // Discrete topology is ambiguous, let's shake things up
      requires_nudge = true;
      break;
    }

    // Find main normal direction
    UnsignedIndex_t max_component_index = 0;
    ScalarType max_component = fabs(face_normal[0]);
    for (UnsignedIndex_t d = 1; d < 3; ++d) {
      if (fabs(face_normal[d]) > max_component) {
        max_component_index = d;
        max_component = fabs(face_normal[d]);
      }
    }
  #ifdef VALDEBUG
  std::cout << "the main component of normal is " << max_component_index << ", with value of " << max_component << std::endl;
  #endif
    // #ifdef NUMERICAL_INTEGRATION
    //     // Find main normal direction
    //     UnsignedIndex_t max_component_index = 0;
    //     ScalarType max_component = fabs(face_normal[0]);
    //     for (UnsignedIndex_t d = 1; d < 3; ++d) {
    //       if (fabs(face_normal[d]) > max_component) {
    //         max_component_index = d;
    //         max_component = fabs(face_normal[d]);
    //       }
    //     }

    //     if (intersection_size == 0) {
    //       if (starting_half_edge->getVertex()->isNotClipped()) {
    //         auto type1_m = full_moments;
    //         auto current_half_edge = starting_half_edge;
    //         auto prev_pt =
    //         current_half_edge->getPreviousVertex()->getLocation(); do {
    //           const auto& curr_pt =
    //           current_half_edge->getVertex()->getLocation(); full_moments +=
    //               computeType1ContributionQuadrature<ReturnType, ScalarType>(
    //                   prev_pt, curr_pt, face_normal, max_component_index);
    //           prev_pt = curr_pt;
    //           current_half_edge = current_half_edge->getNextHalfEdge();
    //         } while (current_half_edge != starting_half_edge);
    //         type1_m = full_moments - type1_m;

    //         auto exact_m =
    //         ReturnType::fromScalarConstant(ReturnScalarType(0)); const auto&
    //         ref_pt = starting_half_edge->getVertex()->getLocation();
    //         current_half_edge =
    //             starting_half_edge->getNextHalfEdge()->getNextHalfEdge();
    //         prev_pt = current_half_edge->getPreviousVertex()->getLocation();
    //         do {
    //           const auto& curr_pt =
    //           current_half_edge->getVertex()->getLocation(); exact_m +=
    //           computeType1Contribution<ReturnType, ScalarType>(
    //               ref_pt, prev_pt, curr_pt);
    //           prev_pt = curr_pt;
    //           current_half_edge = current_half_edge->getNextHalfEdge();
    //         } while (current_half_edge != starting_half_edge);

    //         std::cout << "        Max normal = " << max_component_index
    //                   << std::endl;
    //         std::cout << "Type 1 total exact = " << exact_m << std::endl;
    //         std::cout << "Type 1 total       = " << type1_m << std::endl;
    //       }
    //     }
    // #else
    // This face has not intersections so the moment contribution is only of
    // type 1
    if (intersection_size == 0) {
  #ifdef VALDEBUG
  std::cout << "the face has no intersection" << std::endl;
  #endif
      // The face is entirely below (= unclipped)
      if (starting_half_edge->getVertex()->isNotClipped()) {
        // We need a reference point for the type 1 moment contribution
        const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
        auto current_half_edge = starting_half_edge;
        auto prev_pt = current_half_edge->getPreviousVertex()->getLocation();
        bool skip_first = true;  // This avoid calculating a type 1 equal to 0
  #ifdef VALDEBUG
  std::cout << "it is below the cylinder" << std::endl;
  std::cout << "ref point : " << ref_pt << std::endl;
  #endif
        do {
          const auto& curr_pt = current_half_edge->getVertex()->getLocation();

  #ifdef VALDEBUG
  std::cout << "\ncomputing contribution of edge" << std::endl;
  std::cout << "current point : " << curr_pt << std::endl;
  std::cout << "previous point : " << prev_pt << std::endl;
  std::cout << "face normal : " << face_normal << std::endl;
  #endif
          const auto& type1contribution = computeType1Contribution<ReturnType, ScalarType>(
              ref_pt, prev_pt, curr_pt, &skip_first, false, face_normal,
              max_component_index);
  #ifdef VALDEBUG
  std::cout << "\nType 1 contribution of the face : " << type1contribution << std::endl;
  #endif
          full_moments += type1contribution;
          prev_pt = curr_pt;
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
      }
      else {
  #ifdef VALDEBUG
  std::cout << "it is above the cylinder, skip" << std::endl;
  #endif
      }
    }
    // #endif
    // The face has 2 intersections (i.e. 1 arc): we start from the entry
    else if (intersection_size == 2) {
  #ifdef VALDEBUG
  std::cout << "the face has 2 intersections" << std::endl;
  #endif
      // We need a reference point for the type 1 moment contribution
      const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
  #ifdef VALDEBUG
  std::cout << "the ref point will be : " << ref_pt << std::endl;
  #endif
      // We first traverse straight boundary arcs and return the second
      // intersection (i.e. the exit)
      half_edge_type* exit_half_edge;
      bool skip_first = true;  // This avoid calculating a type 1 equal to 0
          const auto& type1contribution = computeUnclippedSegmentType1Contribution<ReturnType, ScalarType>(
              ref_pt, starting_half_edge, exit_half_edge,
              &skip_first, face_normal, max_component_index);
  #ifdef VALDEBUG
  std::cout << "\nType 1 contribution of the face : " << type1contribution << std::endl;
  #endif
      full_moments += type1contribution;
      // We can now compute the type 2 and 3 moment contribution
          const auto& type3contribution = computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
          a_aligned_cylinder, ref_pt, starting_half_edge, exit_half_edge,
          &skip_first, face_normal, max_component_index, false, &requires_nudge,
          a_surface);
  #ifdef VALDEBUG
  std::cout << "\nType 3 contribution of the face : " << type3contribution << std::endl;
  #endif
      full_moments += type3contribution;
    }
    // The face has more than 2 intersections (i.e. more than 1 arc). We
    // need to discriminate elliptic/hyperbolic/parabolic cases
    else {
  #ifdef VALDEBUG
  std::cout << "the face has more than 2 intersections (i.e. more than 1 arc)" << std::endl;
  #endif
      // These flags identify the type of the conic section arcs in the face
      const bool hyperbolic_face = b < MACHINE_EPSILON;
      const bool rectangle_face = fabs(face_normal[0]) < MACHINE_EPSILON;

      // If the face is convex and we do not want to output the parametrized
      // surface, we don't need to sort the intersections to avoid
      // overlapping conic section arcs
      if (face.isTriangle() &&
          std::is_same_v<SurfaceOutputType, NoSurfaceOutput>) {
        // CASE: the conic section arcs are arcs of an ellipse
        if (!rectangle_face && !hyperbolic_face) {
          // The sign of coeff_a gives us the direction of traversal of the
          // intersection list
          // const bool reverse = a_aligned_paraboloid.a() < ZERO;
          const bool reverse = false;
          // We need a reference point for the type 1 moment contribution
          const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
          half_edge_type* exit_half_edge;
          auto current_edge = starting_half_edge;
          bool skip_first = true;  // This avoid calculating a type 1
                                   // contrib that we know is equal to 0
          std::size_t found_intersections = 0;
          // We traverse the half-edge structure of the face
          do {
            full_moments +=
                computeUnclippedSegmentType1Contribution<ReturnType,
                                                         ScalarType>(
                    ref_pt, current_edge, exit_half_edge,
                    &skip_first, face_normal, max_component_index);

            // From the exit intersection, we move to the next entry and
            // compute type 2 and 3 moment contributions
            skip_first = false;
            if (reverse) {
              full_moments +=
                  computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
                      a_aligned_cylinder, ref_pt, current_edge,
                      exit_half_edge, &skip_first, face_normal,
                      max_component_index, false, &requires_nudge, a_surface);
              current_edge = exit_half_edge->getNextHalfEdge();
              while (current_edge->getVertex()->needsToSeek()) {
                current_edge = current_edge->getNextHalfEdge();
              }
            } else {
              current_edge = exit_half_edge->getNextHalfEdge();
              while (current_edge->getVertex()->needsToSeek()) {
                current_edge = current_edge->getNextHalfEdge();
              }
              full_moments +=
                  computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
                      a_aligned_cylinder, ref_pt, current_edge,
                      exit_half_edge, &skip_first, face_normal,
                      max_component_index, false, &requires_nudge, a_surface);
            }
            found_intersections += 2;
          } while (found_intersections != intersection_size);
        }
        // CASE: the intersection are parrallllele lines
        else {
          // i don't know how to do that for now :)
          // let's just nudge
          requires_nudge = true;
          break;
        }
        // else {
        //   // The half-edges ending at an intersection will be stores in this
        //   // vector
        //   SmallVector<half_edge_type*, 6> intersections;
        //   intersections.resize(0);
        //   // Find intersections and add to list
        //   auto current_edge = starting_half_edge;
        //   bool reverse_order = false;
        //   half_edge_type* exit_half_edge;
        //   std::size_t found_intersections = 0;
        //   const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
        //   bool skip_first = true;
        //   do {
        //     full_moments +=
        //         computeUnclippedSegmentType1Contribution<ReturnType,
        //                                                  ScalarType>(
        //             ref_pt, current_edge, exit_half_edge,
        //             &skip_first, face_normal, max_component_index);
        //     intersections.push_back(current_edge);
        //     intersections.push_back(exit_half_edge);
        //     current_edge = exit_half_edge->getNextHalfEdge();
        //     while (current_edge->getVertex()->needsToSeek()) {
        //       current_edge = current_edge->getNextHalfEdge();
        //     }
        //     found_intersections += 2;
        //     skip_first = false;
        //   } while (found_intersections != intersection_size);

        //   // Now, we discriminate hyperbolic and parabolic cases
        //   if (hyperbolic_face) {
        //     // We need the conic center to determine nappe ownership. Note
        //     // that, by this point, normal[2] has to be different than 0
        //     const std::array<ScalarType, 2> conic_center{
        //         {face_normal[0] /
        //              (TWO * a_aligned_paraboloid.a() * face_normal[2]),
        //          face_normal[1] /
        //              (TWO * a_aligned_paraboloid.b() * face_normal[2])}};
        //     // We need some point that lives in the plane of the face (e.g.,
        //     // the first intersection)
        //     const auto& pt_in_plane =
        //         intersections[0]->getVertex()->getLocation().getPt();
        //     // The next lines calculate whether the conic center on the face
        //     // lies above or below the paraboloid
        //     const ScalarType delta_face = (face_normal[0] * pt_in_plane[0] +
        //                                    face_normal[1] * pt_in_plane[1] +
        //                                    face_normal[2] * pt_in_plane[2]) /
        //                                   face_normal[2];
        //     const ScalarType gamma_face =
        //         a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
        //         a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] -
        //         delta_face;
        //     const ScalarType z_diff =
        //         face_normal[0] * conic_center[0] / face_normal[2] +
        //         face_normal[1] * conic_center[1] / face_normal[2] - gamma_face -
        //         TWO * delta_face;
        //     const std::size_t split_ind =
        //         a_aligned_paraboloid.a() * gamma_face > ZERO ? 0 : 1;

        //     // Check if all intersections belong to the same branch of the
        //     // hyperbola
        //     bool vertices_on_same_branch = true;
        //     for (std::size_t i = 1; i < intersection_size; ++i) {
        //       const auto& curr_pt =
        //           intersections[i]->getVertex()->getLocation().getPt();
        //       if ((pt_in_plane[split_ind] - conic_center[split_ind]) *
        //               (curr_pt[split_ind] - conic_center[split_ind]) <
        //           ZERO) {
        //         vertices_on_same_branch = false;
        //         break;
        //       }
        //     }

        //     // Update traversal order
        //     reverse_order = z_diff < ZERO ? !vertices_on_same_branch
        //                                   : vertices_on_same_branch;
        //   }
        //   // The conic section arcs are arcs of a parabola
        //   else {
        //     // The polygon formed by the intersection must be convex, so the
        //     // average of the intersections must lie inside that polygon
        //     Pt intersection_avg = Pt(ZERO, ZERO, ZERO);
        //     for (std::size_t i = 0; i < intersection_size; ++i) {
        //       const auto& curr_pt =
        //           intersections[i]->getVertex()->getLocation().getPt();
        //       intersection_avg += curr_pt;
        //     }
        //     intersection_avg /= ScalarType(intersection_size);

        //     // Is this average point below or above the paraboloid?
        //     const ScalarType z_diff =
        //         intersection_avg[2] +
        //         a_aligned_paraboloid.a() * intersection_avg[0] *
        //             intersection_avg[0] +
        //         a_aligned_paraboloid.b() * intersection_avg[1] *
        //             intersection_avg[1];

        //     // Update traversal order
        //     reverse_order = (z_diff > ZERO);
        //   }

        //   // Traverse face from entry->exit until all intersections have
        //   // been traversed
        //   bool first_entry = true;
        //   // Loop over vertices, only use entry vertices
        //   for (std::size_t i = 0; i < intersection_size; i += 2) {
        //     // Identify entry vertex
        //     auto entry_half_edge = intersections[i];
        //     // Loop over internal portion from entry->exit
        //     int exit_index;
        //     if (reverse_order) {
        //       exit_index = (i != intersection_size - 1) ? i + 1 : 0;
        //     } else {
        //       exit_index = (i != 0) ? i - 1 : intersection_size - 1;
        //     }
        //     auto exit_half_edge = intersections[exit_index];
        //     full_moments +=
        //         computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
        //             a_aligned_paraboloid, ref_pt, entry_half_edge,
        //             exit_half_edge, &first_entry, face_normal,
        //             max_component_index, false, &requires_nudge, a_surface);
        //     first_entry = false;
        //   }
        //   // Clear list of intersections
        //   intersections.clear();
        // }
      }
      // The face is not convex, or we want to output the clipped surface.
      // We then need to sort the intersections
      else {
  #ifdef VALDEBUG
  std::cout << "we are going to sort the intersection" << std::endl;
  #endif
        // This flag is there to prevent degenerate configuration with
        // overlapping arcs
        bool ignore_type3_contributions = false;
        // The intersections and some scalar needed for sorting will be
        // stored in this vector
        using stype = std::pair<half_edge_type*, ScalarType>;
        SmallVector<stype, 6> intersections;
        intersections.resize(0);
        // Find intersections and determine status
        auto current_edge = starting_half_edge;
        bool reverse_order = false;
        half_edge_type* exit_half_edge;
        std::size_t found_intersections = 0;
        const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
  #ifdef VALDEBUG
  std::cout << "we starting with the edge : " << *current_edge << std::endl;
  std::cout << "the ref point is : " << ref_pt << std::endl;
  #endif
        bool skip_first = true;
        do {
          auto type1contribution = computeUnclippedSegmentType1Contribution<ReturnType, ScalarType>(
                  ref_pt, current_edge, exit_half_edge,
                  &skip_first, face_normal, max_component_index);
  #ifdef VALDEBUG
  std::cout << "type 1 contribution : " << type1contribution << std::endl;
  #endif
          full_moments += type1contribution;
          intersections.push_back(std::pair<half_edge_type*, ScalarType>(
              {current_edge, ScalarType(0)}));
          intersections.push_back(std::pair<half_edge_type*, ScalarType>(
              {exit_half_edge, ScalarType(0)}));
  #ifdef VALDEBUG
  std::cout << "adding intersections : " << *current_edge << std::endl;
  std::cout << "and exit : " << *exit_half_edge << std::endl;
  #endif
          current_edge = exit_half_edge->getNextHalfEdge();
  #ifdef VALDEBUG
  std::cout << "next edge is : " << *current_edge << std::endl;
  #endif
          while (current_edge->getVertex()->needsToSeek()) {
  #ifdef VALDEBUG
  std::cout << "next edge is is skiped " << *current_edge << std::endl;
  #endif
            current_edge = current_edge->getNextHalfEdge();
  #ifdef VALDEBUG
  std::cout << "next edge is  : " << *current_edge << std::endl;
  #endif
          }
          found_intersections += 2;
          skip_first = false;
        } while (found_intersections != intersection_size);
  #ifdef VALDEBUG
  std::cout << "ended with " << found_intersections << " intersections" << std::endl;
  #endif

        // Here also, we need to discriminate elliptic/hyperbolic/parabolic
        // conic sections First: we start with elliptic conic section arcs
        if (true) {
          // Compute the center of this ellipse
          // const std::array<ScalarType, 2> conic_center{
          //     { face.getPlane().distance() /
          //          face_normal[0],
          //       ScalarType(0)}};
  #ifdef VALDEBUG
  std::cout << "we don't have rectangle intersection" << std::endl;
  std::cout << "these are our intersections (befor sorting) : " << std::endl;
  #endif
          // const ScalarType invert =
          //     a_aligned_paraboloid.a() < ZERO ? -normal_invert : normal_invert;
          // Store angular position of intersection on the ellipse
          const ScalarType invert = rectangle_face ? copysign(ONE, face_normal[2]) : copysign(ONE, face_normal[0]);
          if (rectangle_face) {
            for (auto& element : intersections) {
              const auto& pt = element.first->getVertex()->getLocation().getPt();
              element.second = invert * atan2(pt[1],
                                              pt[0]);
    #ifdef VALDEBUG
    std::cout << "point : " << pt << ", sorting value : " << element.second << std::endl;
    #endif
            }
          } else {
            for (auto& element : intersections) {
              const auto& pt = element.first->getVertex()->getLocation().getPt();
              element.second = invert * atan2(pt[2],
                                              pt[1]);
    #ifdef VALDEBUG
    std::cout << "point : " << pt << ", sorting value : " << element.second << std::endl;
    #endif
            }
          }
          // Sort intersections
          std::sort(intersections.begin(), intersections.end(),
                    [](const stype& a, const stype& b) {
                      return a.second < b.second;
                    });

          // Backdoor for case where two angles are equal (because
          // of floating-point errors)
          bool restart_sort = false;
          for (std::size_t i = 0; i < intersection_size - 1; ++i) {
            if (fabs(intersections[i].second - intersections[i + 1].second) <
                ONE_HUNDRED * MACHINE_EPSILON) {
              restart_sort = true;
              break;
            }
          }
          // restart_sort = true;
          // If this sorting failed, we now sort based on Y positions
          if (restart_sort) {
            SmallVector<stype, 6> intersection_copy;
            intersection_copy.resize(intersections.size());
            // We split in X direction
            UnsignedIndex_t split_ind = 0;
            UnsignedIndex_t store_ind = 1;
            UnsignedIndex_t pos_end = 0;
            UnsignedIndex_t neg_end = 0;
            // Each half of the ellipse (X positive and X negative) will be
            // sorted separately
            for (auto& element : intersections) {
              const auto& pt =
                  element.first->getVertex()->getLocation().getPt();
              if (pt[split_ind] > ZERO) {
                const auto loc = pos_end++;
                intersection_copy[loc].first = element.first;
                intersection_copy[loc].second = invert * pt[store_ind];
              } else {
                const auto loc = intersection_size - 1 - (neg_end++);
                intersection_copy[loc].first = element.first;
                intersection_copy[loc].second = -invert * pt[store_ind];
              }
            }
            // Sort each ellipse half
            std::sort(intersection_copy.begin(),
                      intersection_copy.begin() + pos_end,
                      [](const stype& a, const stype& b) {
                        return a.second < b.second;
                      });
            std::sort(intersection_copy.begin() + pos_end,
                      intersection_copy.end(),
                      [](const stype& a, const stype& b) {
                        return a.second < b.second;
                      });
            intersections = intersection_copy;
            // Let's check if that sort was successful
            bool re_restart_sort = false;
            if (pos_end > 0) {
              for (UnsignedIndex_t i = 0; i < pos_end - 1; ++i) {
                if (fabs(intersections[i].second -
                         intersections[i + 1].second) <
                    ONE_HUNDRED * MACHINE_EPSILON) {
                  re_restart_sort = true;
                  break;
                }
              }
            }
            for (UnsignedIndex_t i = pos_end; i < intersection_size - 1; ++i) {
              if (fabs(intersections[i].second - intersections[i + 1].second) <
                  ONE_HUNDRED * MACHINE_EPSILON) {
                re_restart_sort = true;
                break;
              }
            }

            // If this sorting failed, we now sort based on X positions
            if (re_restart_sort) {
              // We split in Y direction
              split_ind = 1;
              store_ind = 0;
              pos_end = 0;
              neg_end = 0;
              // Each half of the ellipse (Y positive and Y negative) will
              // be sorted separately
              for (auto& element : intersections) {
                const auto& pt =
                    element.first->getVertex()->getLocation().getPt();
                if (pt[split_ind] > ZERO) {
                  const auto loc = pos_end++;
                  intersection_copy[loc].first = element.first;
                  intersection_copy[loc].second = -invert * pt[store_ind];
                } else {
                  const auto loc = intersection_size - 1 - (neg_end++);
                  intersection_copy[loc].first = element.first;
                  intersection_copy[loc].second = invert * pt[store_ind];
                }
              }
              // Sort each ellipse half
              std::sort(intersection_copy.begin(),
                        intersection_copy.begin() + pos_end,
                        [](const stype& a, const stype& b) {
                          return a.second < b.second;
                        });
              std::sort(intersection_copy.begin() + pos_end,
                        intersection_copy.end(),
                        [](const stype& a, const stype& b) {
                          return a.second < b.second;
                        });
              intersections = intersection_copy;
              // Let's check if that sort was successful
              bool re_re_restart_sort = false;
              if (pos_end > 0) {
                for (UnsignedIndex_t i = 0; i < pos_end - 1; ++i) {
                  if (fabs(intersections[i].second -
                           intersections[i + 1].second) <
                      ONE_HUNDRED * MACHINE_EPSILON) {
                    re_re_restart_sort = true;
                    break;
                  }
                }
              }
              for (UnsignedIndex_t i = pos_end; i < intersection_size - 1;
                   ++i) {
                if (fabs(intersections[i].second -
                         intersections[i + 1].second) <
                    ONE_HUNDRED * MACHINE_EPSILON) {
                  re_re_restart_sort = true;
                  break;
                }
              }
              // If all sorts have failed, switch to QP and try again
              if (re_re_restart_sort) {
                requires_nudge = true;
                break;
              }
            }
          }
        }
        // The section is a hyperbolique
        else 
        {
          // pour l'instant on skip
          requires_nudge = true;
          break;
        }
        // if (hyperbolic_face) {
        //   // We will sort each hyperbola branch separately
        //   std::size_t pos_end = 0;
        //   std::size_t neg_end = 0;
        //   SmallVector<stype, 6> intersection_copy;
        //   intersection_copy.resize(intersections.size());
        //   // Compute center of the hyperbola
        //   const std::array<ScalarType, 2> conic_center{
        //       {face_normal[0] /
        //            (TWO * a_aligned_paraboloid.a() * face_normal[2]),
        //        face_normal[1] /
        //            (TWO * a_aligned_paraboloid.b() * face_normal[2])}};
        //   // We need some point that lives in the plane of the face (e.g.,
        //   // the first intersection)
        //   const auto& pt_in_plane =
        //       intersections[0].first->getVertex()->getLocation().getPt();
        //   const ScalarType delta_face = (face_normal[0] * pt_in_plane[0] +
        //                                  face_normal[1] * pt_in_plane[1] +
        //                                  face_normal[2] * pt_in_plane[2]) /
        //                                 face_normal[2];
        //   const ScalarType gamma_face =
        //       a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
        //       a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] -
        //       delta_face;
        //   const std::size_t split_ind =
        //       a_aligned_paraboloid.a() * gamma_face > ZERO ? 0 : 1;
        //   const std::size_t store_ind = split_ind == 0 ? 1 : 0;
        //   const ScalarType z_center_plane =
        //       -face_normal[0] * conic_center[0] / face_normal[2] -
        //       face_normal[1] * conic_center[1] / face_normal[2] + delta_face;
        //   const ScalarType z_center_paraboloid =
        //       -a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] -
        //       a_aligned_paraboloid.b() * conic_center[1] * conic_center[1];
        //   ScalarType total_invert = copysign(ONE, face_normal[2]);
        //   total_invert = z_center_plane < z_center_paraboloid ? total_invert
        //                                                       : -total_invert;
        //   total_invert = split_ind == 0 ? total_invert : -total_invert;
        //   // Prepare sorting scalar for each hyperbola branch
        //   for (auto& element : intersections) {
        //     const auto& pt = element.first->getVertex()->getLocation().getPt();
        //     if (pt[split_ind] > conic_center[split_ind]) {
        //       const auto loc = pos_end++;
        //       intersection_copy[loc].first = element.first;
        //       intersection_copy[loc].second = total_invert * pt[store_ind];
        //     } else {
        //       const auto loc = intersection_size - 1 - (neg_end++);
        //       intersection_copy[loc].first = element.first;
        //       intersection_copy[loc].second = -total_invert * pt[store_ind];
        //     }
        //   }
        //   // Sort each hyperbola branch separately
        //   assert(pos_end + neg_end == intersection_size);
        //   std::sort(intersection_copy.begin(),
        //             intersection_copy.begin() + pos_end,
        //             [](const stype& a, const stype& b) {
        //               return a.second < b.second;
        //             });
        //   std::sort(intersection_copy.begin() + pos_end,
        //             intersection_copy.end(),
        //             [](const stype& a, const stype& b) {
        //               return a.second < b.second;
        //             });
        //   // If the determination of the traversal order was too close for
        //   // confort, we swith to QP
        //   if (fabs(gamma_face) < ONE_HUNDRED * MACHINE_EPSILON ||
        //       fabs(z_center_plane - z_center_paraboloid) <
        //           ONE_HUNDRED * MACHINE_EPSILON) {
        //     requires_nudge = true;
        //     break;
        //   }
        //   // If the sorted intersections are too close to each other, we
        //   // switch to QP
        //   else {
        //     for (std::size_t i = 0; i < intersection_size - 1; ++i) {
        //       if (fabs(intersection_copy[i].second -
        //                intersection_copy[i + 1].second) <
        //           ONE_HUNDRED * MACHINE_EPSILON) {
        //         requires_nudge = true;
        //         break;
        //       }
        //     }
        //   }
        //   intersections = intersection_copy;
        // }
        // // The conic section arc is a parabola
        // else {
        //   // We exploit the convexity of the parabola: the average of the
        //   // intersections must lie inside the polygon of intersectionsd
        //   Pt intersection_avg = Pt(ZERO, ZERO, ZERO);
        //   for (std::size_t i = 0; i < intersection_size; ++i) {
        //     const auto& curr_pt =
        //         intersections[i].first->getVertex()->getLocation().getPt();
        //     intersection_avg += curr_pt;
        //   }
        //   intersection_avg /= ScalarType(intersection_size);
        //   // Is this average point above of below the paraboloid?
        //   const ScalarType z_diff =
        //       intersection_avg[2] +
        //       a_aligned_paraboloid.a() * intersection_avg[0] *
        //           intersection_avg[0] +
        //       a_aligned_paraboloid.b() * intersection_avg[1] *
        //           intersection_avg[1];
        //   // If we can't determine that accurately enough, we switch to QP
        //   if (fabs(z_diff) < ONE_HUNDRED * MACHINE_EPSILON) {
        //     requires_nudge = true;
        //     break;
        //   }
        //   // Find dominant face normal direction
        //   std::size_t dir = 0;
        //   ScalarType max_normal = fabs(face_normal[dir]);
        //   for (std::size_t d = 1; d < 3; ++d) {
        //     if (fabs(face_normal[d]) > max_normal) {
        //       dir = d;
        //       max_normal = fabs(face_normal[d]);
        //     }
        //   }
        //   const ScalarType normal_invert = copysign(ONE, face_normal[dir]);
        //   const ScalarType invert =
        //       z_diff < ZERO ? normal_invert : -normal_invert;
        //   // Compute convex hull on projected plane
        //   const std::size_t x_id = (dir + 1) % 3;
        //   const std::size_t y_id = (dir + 2) % 3;
        //   // Find the leftmost point
        //   std::size_t left_id = 0;
        //   auto left_pt =
        //       intersections[left_id].first->getVertex()->getLocation().getPt();
        //   for (UnsignedIndex_t i = 1; i < intersection_size; ++i) {
        //     const auto pt =
        //         intersections[i].first->getVertex()->getLocation().getPt();
        //     if (pt[x_id] < left_pt[x_id]) {
        //       left_id = i;
        //       left_pt = pt;
        //     }
        //   }
        //   // Start from leftmost point, keep moving
        //   // counterclockwise until reach the start point again.
        //   std::size_t p = left_id;
        //   auto p_pt = left_pt;
        //   std::size_t hull_size = 0;
        //   SmallVector<stype, 6> intersection_copy;
        //   intersection_copy.resize(intersections.size());
        //   bool is_flat = false;
        //   bool has_flat = false;
        //   do {
        //     p_pt = intersections[p].first->getVertex()->getLocation().getPt();
        //     intersection_copy[hull_size++] = intersections[p];
        //     if (hull_size == intersection_size) {
        //       break;
        //     }
        //     std::size_t q = (p + 1) % intersection_size;
        //     std::size_t flatness_counter = 0;
        //     for (std::size_t i = 0; i < intersection_size; ++i) {
        //       if (q != i && p != i) {
        //         const auto i_pt =
        //             intersections[i].first->getVertex()->getLocation().getPt();
        //         const auto q_pt =
        //             intersections[q].first->getVertex()->getLocation().getPt();
        //         ScalarType dot_product =
        //             (i_pt[y_id] - p_pt[y_id]) * (q_pt[x_id] - i_pt[x_id]) -
        //             (i_pt[x_id] - p_pt[x_id]) * (q_pt[y_id] - i_pt[y_id]);
        //         // If the points are aligned, keep closest and
        //         // increase flatness counter
        //         if (fabs(dot_product) < ANGLE_EPSILON) {
        //           has_flat = true;
        //           flatness_counter++;
        //         } else if (dot_product < ZERO) {
        //           q = i;
        //         }
        //       }
        //     }
        //     // Check if all intersections are aligned
        //     if (flatness_counter == intersection_size - 2) {
        //       is_flat = true;
        //       break;
        //     }
        //     p = q;
        //   } while (p != left_id);
        //   // If the convex-hull of the intersections does not have a 0 area
        //   // (i.e. is flat)
        //   if (!is_flat) {
        //     // Are some intersection aligned though?
        //     if (!is_flat && has_flat) {
        //       // Let's sort based on angular position
        //       for (auto& element : intersections) {
        //         const auto& pt =
        //             element.first->getVertex()->getLocation().getPt();
        //         element.second =
        //             invert *
        //             copysign(
        //                 ONE - (pt[x_id] - intersection_avg[x_id]) /
        //                           (fabs(pt[x_id] - intersection_avg[x_id]) +
        //                            fabs(pt[y_id] - intersection_avg[y_id])),
        //                 (pt[y_id] - intersection_avg[y_id]));
        //       }
        //       std::sort(intersections.begin(), intersections.end(),
        //                 [](const stype& a, const stype& b) {
        //                   return a.second < b.second;
        //                 });
        //       // If there remains aligned intersections, switch to QP
        //       for (std::size_t i = 0; i < intersection_size - 1; ++i) {
        //         if (fabs(intersections[i].second -
        //                  intersections[i + 1].second) <if (hyperbolic_face) {
        //   // We will sort each hyperbola branch separately
        //   std::size_t pos_end = 0;
        //   std::size_t neg_end = 0;
        //   SmallVector<stype, 6> intersection_copy;
        //   intersection_copy.resize(intersections.size());
        //   // Compute center of the hyperbola
        //   const std::array<ScalarType, 2> conic_center{
        //       {face_normal[0] /
        //            (TWO * a_aligned_paraboloid.a() * face_normal[2]),
        //        face_normal[1] /
        //            (TWO * a_aligned_paraboloid.b() * face_normal[2])}};
        //   // We need some point that lives in the plane of the face (e.g.,
        //   // the first intersection)
        //   const auto& pt_in_plane =
        //       intersections[0].first->getVertex()->getLocation().getPt();
        //   const ScalarType delta_face = (face_normal[0] * pt_in_plane[0] +
        //                                  face_normal[1] * pt_in_plane[1] +
        //                                  face_normal[2] * pt_in_plane[2]) /
        //                                 face_normal[2];
        //   const ScalarType gamma_face =
        //       a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
        //       a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] -
        //       delta_face;
        //   const std::size_t split_ind =
        //       a_aligned_paraboloid.a() * gamma_face > ZERO ? 0 : 1;
        //   const std::size_t store_ind = split_ind == 0 ? 1 : 0;
        //   const ScalarType z_center_plane =
        //       -face_normal[0] * conic_center[0] / face_normal[2] -
        //       face_normal[1] * conic_center[1] / face_normal[2] + delta_face;
        //   const ScalarType z_center_paraboloid =
        //       -a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] -
        //       a_aligned_paraboloid.b() * conic_center[1] * conic_center[1];
        //   ScalarType total_invert = copysign(ONE, face_normal[2]);
        //   total_invert = z_center_plane < z_center_paraboloid ? total_invert
        //                                                       : -total_invert;
        //   total_invert = split_ind == 0 ? total_invert : -total_invert;
        //   // Prepare sorting scalar for each hyperbola branch
        //   for (auto& element : intersections) {
        //     const auto& pt = element.first->getVertex()->getLocation().getPt();
        //     if (pt[split_ind] > conic_center[split_ind]) {
        //       const auto loc = pos_end++;
        //       intersection_copy[loc].first = element.first;
        //       intersection_copy[loc].second = total_invert * pt[store_ind];
        //     } else {
        //       const auto loc = intersection_size - 1 - (neg_end++);
        //       intersection_copy[loc].first = element.first;
        //       intersection_copy[loc].second = -total_invert * pt[store_ind];
        //     }
        //   }
        //   // Sort each hyperbola branch separately
        //   assert(pos_end + neg_end == intersection_size);
        //   std::sort(intersection_copy.begin(),
        //             intersection_copy.begin() + pos_end,
        //             [](const stype& a, const stype& b) {
        //               return a.second < b.second;
        //             });
        //   std::sort(intersection_copy.begin() + pos_end,
        //             intersection_copy.end(),
        //             [](const stype& a, const stype& b) {
        //               return a.second < b.second;
        //             });
        //   // If the determination of the traversal order was too close for
        //   // confort, we swith to QP
        //   if (fabs(gamma_face) < ONE_HUNDRED * MACHINE_EPSILON ||
        //       fabs(z_center_plane - z_center_paraboloid) <
        //           ONE_HUNDRED * MACHINE_EPSILON) {
        //     requires_nudge = true;
        //     break;
        //   }
        //   // If the sorted intersections are too close to each other, we
        //   // switch to QP
        //   else {
        //     for (std::size_t i = 0; i < intersection_size - 1; ++i) {
        //       if (fabs(intersection_copy[i].second -
        //                intersection_copy[i + 1].second) <
        //           ONE_HUNDRED * MACHINE_EPSILON) {
        //         requires_nudge = true;
        //         break;
        //       }
        //     }
        //   }
        //   intersections = intersection_copy;
        // }
        // // The conic section arc is a parabola
        // else {
        //   // We exploit the convexity of the parabola: the average of the
        //   // intersections must lie inside the polygon of intersectionsd
        //   Pt intersection_avg = Pt(ZERO, ZERO, ZERO);
        //   for (std::size_t i = 0; i < intersection_size; ++i) {
        //     const auto& curr_pt =
        //         intersections[i].first->getVertex()->getLocation().getPt();
        //     intersection_avg += curr_pt;
        //   }
        //   intersection_avg /= ScalarType(intersection_size);
        //   // Is this average point above of below the paraboloid?
        //   const ScalarType z_diff =
        //       intersection_avg[2] +
        //       a_aligned_paraboloid.a() * intersection_avg[0] *
        //           intersection_avg[0] +
        //       a_aligned_paraboloid.b() * intersection_avg[1] *
        //           intersection_avg[1];
        //   // If we can't determine that accurately enough, we switch to QP
        //   if (fabs(z_diff) < ONE_HUNDRED * MACHINE_EPSILON) {
        //     requires_nudge = true;
        //     break;
        //   }
        //   // Find dominant face normal direction
        //   std::size_t dir = 0;
        //   ScalarType max_normal = fabs(face_normal[dir]);
        //   for (std::size_t d = 1; d < 3; ++d) {
        //     if (fabs(face_normal[d]) > max_normal) {
        //       dir = d;
        //       max_normal = fabs(face_normal[d]);
        //     }
        //   }
        //   const ScalarType normal_invert = copysign(ONE, face_normal[dir]);
        //   const ScalarType invert =
        //       z_diff < ZERO ? normal_invert : -normal_invert;
        //   // Compute convex hull on projected plane
        //   const std::size_t x_id = (dir + 1) % 3;
        //   const std::size_t y_id = (dir + 2) % 3;
        //   // Find the leftmost point
        //   std::size_t left_id = 0;
        //   auto left_pt =
        //       intersections[left_id].first->getVertex()->getLocation().getPt();
        //   for (UnsignedIndex_t i = 1; i < intersection_size; ++i) {
        //     const auto pt =
        //         intersections[i].first->getVertex()->getLocation().getPt();
        //     if (pt[x_id] < left_pt[x_id]) {
        //       left_id = i;
        //       left_pt = pt;
        //     }
        //   }
        //   // Start from leftmost point, keep moving
        //   // counterclockwise until reach the start point again.
        //   std::size_t p = left_id;
        //   auto p_pt = left_pt;
        //   std::size_t hull_size = 0;
        //   SmallVector<stype, 6> intersection_copy;
        //   intersection_copy.resize(intersections.size());
        //   bool is_flat = false;
        //   bool has_flat = false;
        //   do {
        //     p_pt = intersections[p].first->getVertex()->getLocation().getPt();
        //     intersection_copy[hull_size++] = intersections[p];
        //     if (hull_size == intersection_size) {
        //       break;
        //     }
        //     std::size_t q = (p + 1) % intersection_size;
        //     std::size_t flatness_counter = 0;
        //     for (std::size_t i = 0; i < intersection_size; ++i) {
        //       if (q != i && p != i) {
        //         const auto i_pt =
        //             intersections[i].first->getVertex()->getLocation().getPt();
        //         const auto q_pt =
        //             intersections[q].first->getVertex()->getLocation().getPt();
        //         ScalarType dot_product =
        //             (i_pt[y_id] - p_pt[y_id]) * (q_pt[x_id] - i_pt[x_id]) -
        //             (i_pt[x_id] - p_pt[x_id]) * (q_pt[y_id] - i_pt[y_id]);
        //         // If the points are aligned, keep closest and
        //         // increase flatness counter
        //         if (fabs(dot_product) < ANGLE_EPSILON) {
        //           has_flat = true;
        //           flatness_counter++;
        //         } else if (dot_product < ZERO) {
        //           q = i;
        //         }
        //       }
        //     }
        //     // Check if all intersections are aligned
        //     if (flatness_counter == intersection_size - 2) {
        //       is_flat = true;
        //       break;
        //     }
        //     p = q;
        //   } while (p != left_id);
        //   // If the convex-hull of the intersections does not have a 0 area
        //   // (i.e. is flat)
        //   if (!is_flat) {
        //     // Are some intersection aligned though?
        //     if (!is_flat && has_flat) {
        //       // Let's sort based on angular position
        //       for (auto& element : intersections) {
        //         const auto& pt =
        //             element.first->getVertex()->getLocation().getPt();
        //         element.second =
        //             invert *
        //             copysign(
        //                 ONE - (pt[x_id] - intersection_avg[x_id]) /
        //                           (fabs(pt[x_id] - intersection_avg[x_id]) +
        //                            fabs(pt[y_id] - intersection_avg[y_id])),
        //                 (pt[y_id] - intersection_avg[y_id]));
        //       }
        //       std::sort(intersections.begin(), intersections.end(),
        //                 [](const stype& a, const stype& b) {
        //                   return a.second < b.second;
        //                 });
        //       // If there remains aligned intersections, switch to QP
        //       for (std::size_t i = 0; i < intersection_size - 1; ++i) {
        //         if (fabs(intersections[i].second -
        //                  intersections[i + 1].second) <
        //             ONE_HUNDRED * MACHINE_EPSILON) {
        //           requires_nudge = true;
        //         }
        //       }
        //     } else {
        //       // Is convex hull incomplete? Then error out
        //       if (hull_size != intersection_size) {
        //         requires_nudge = true;
        //       }
        //       // Else: update intersection with ordered list
        //       if (invert > ZERO) {
        //         intersections = intersection_copy;
        //       } else {
        //         for (std::size_t i = 0; i < intersection_size; ++i) {
        //           intersections[i] =
        //               intersection_copy[intersection_size - 1 - i];
        //         }
        //       }
        //     }
        //   } else {
        //     requires_nudge = true;
        //     break;
        //   }
        // }
        //         }
        //       }
        //     } else {
        //       // Is convex hull incomplete? Then error out
        //       if (hull_size != intersection_size) {
        //         requires_nudge = true;
        //       }
        //       // Else: update intersection with ordered list
        //       if (invert > ZERO) {
        //         intersections = intersection_copy;
        //       } else {
        //         for (std::size_t i = 0; i < intersection_size; ++i) {
        //           intersections[i] =
        //               intersection_copy[intersection_size - 1 - i];
        //         }
        //       }
        //     }
        //   } else {
        //     requires_nudge = true;
        //     break;
        //   }
        // }
        // Traverse face from entry->exit until all intersections have been
        // traversed
        auto prev_vertex = intersections[0].first->getPreviousVertex();
        auto entry_half_edge = intersections[0].first;
        bool entry_first =
            (prev_vertex->isClipped() ||
             (prev_vertex->doesNotNeedToSeek() &&
              entry_half_edge->getNextHalfEdge()->getVertex()->isNotClipped()));
        std::size_t start_id = entry_first ? 0 : 1;
        entry_first = false;
        for (std::size_t i = start_id; i < intersection_size; i += 2) {
          // Identify entry vertex
          const auto entry_half_edge = intersections[i].first;
          // Loop over internal portion from entry->exit
          const UnsignedIndex_t exit_index =
              i != 0 ? i - 1 : intersection_size - 1;
          const auto exit_half_edge = intersections[exit_index].first;
          full_moments +=
              computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
                  a_aligned_cylinder, ref_pt, entry_half_edge, exit_half_edge,
                  &entry_first, face_normal, max_component_index,
                  ignore_type3_contributions, &requires_nudge, a_surface);
        }
        // Clear list of intersections
        intersections.clear();
      }
    }
    // #endif
    // If some ambiguous configuration has been detected, switch to QP and
    // shake things up
    if (requires_nudge) {
      break;
    }
  }  // End loop over faces

  // Last call for the nudge
  if (requires_nudge) {
    resetPolyhedron(a_polytope, a_complete_polytope);
    assert(a_polytope->getNumberOfVertices() == starting_number_of_vertices);

    // Nudge and try again!
    return reformQuadraticIntersectionBases<ReturnType>(
        a_polytope, a_complete_polytope, a_aligned_cylinder, a_nudge_iter,
        a_surface);
  }

  return full_moments;
}
}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_INTERSECTION_TPP_
