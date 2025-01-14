// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_INTERSECTION_H_
#define IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_INTERSECTION_H_

#include <float.h>
#include <cassert>
#include <cmath>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/generic_cutting/quadratic_intersection/surface_output.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/mymath.h"
#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/quadratic_reconstruction/rational_bezier_arc.h"

namespace IRL {

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class CylinderType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithCylinder(SegmentedHalfEdgePolyhedronType* a_polytope,
                                  HalfEdgePolytopeType* a_complete_polytope,
                                  const CylinderType& a_cylinder);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AlignedCylinderType,
          class ScalarType, class SurfaceOutputType = NoSurfaceOutput>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithAlignedCylinder(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AlignedCylinderType& a_cylinder,
    const ScalarType a_inv_volume_scale,
    SurfaceOutputType* a_surface = nullptr);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AligneCylinderType,
          class ScalarType, class SurfaceOutputType = NoSurfaceOutput>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
formCylinderIntersectionBases(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AligneCylinderType& a_aligned_cylinder,
    const UnsignedIndex_t a_nudge_iter, const ScalarType a_inv_volume_scale,
    SurfaceOutputType* a_surface = nullptr);

template <class ScalarType>
inline NormalBase<ScalarType> computeTangentVectorAtPoint(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_plane_normal,
    const PtBase<ScalarType>& a_pt);

template <class ScalarType>
inline NormalBase<ScalarType> computeAndCorrectTangentVectorAtPt(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_plane_normal,
    const PtBase<ScalarType>& a_origin_pt, const PtBase<ScalarType>& a_end_pt,
    const NormalBase<ScalarType>& a_end_tangent,
    const PtBase<ScalarType>& a_intersection_pt);

// template <class ScalarType, class PtTypeWithGradient>
// inline PtTypeWithGradient computeAndCorrectTangentVectorAndGradientAtPt(
//     const AlignedParaboloidBase<ScalarType>& a_paraboloid,
//     const PtTypeWithGradient& a_plane_normal,
//     const PtTypeWithGradient& a_origin_pt, const PtTypeWithGradient& a_end_pt,
//     const PtTypeWithGradient& a_end_tangent,
//     const PtTypeWithGradient& a_intersection_pt);

template <class ReturnType, class ScalarType,
          class SurfaceOutputType = NoSurfaceOutput, class PtType>
ReturnType computeType3ContributionWithSplit(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_plane_normal, const PtType& a_pt_ref,
    const PtType& a_pt_0, const PtType& a_pt_1,
    const NormalBase<ScalarType>& a_tangent_0,
    const NormalBase<ScalarType>& a_tangent_1, bool* a_requires_nudge,
    UnsignedIndex_t* a_split_counter, SurfaceOutputType* a_surface = nullptr);

template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType = NoSurfaceOutput, class PtType>
ReturnType computeNewEdgeSegmentContribution(
    const AlignedCylinderBase<ScalarType>& a_aligned_cylinder,
    const PtType& a_ref_pt, const HalfEdgeType a_entry_half_edge,
    const HalfEdgeType a_exit_half_edge, const bool skip_first,
    const bool a_ignore_type3, bool* a_requires_nudge,
    SurfaceOutputType* a_surface = nullptr);

template <class PtType, class ScalarType>
void checkAndFindIntercepts(
    const AlignedCylinderBase<ScalarType>& a_cylinder, const PtType& a_pt_0,
    const PtType& a_pt_1, StackVector<PtType, 2>* a_intercepts,
    const ScalarType a_nudge_epsilon, const bool a_elliptic);

template <class VertexType>
bool vertexBelow(
    const VertexType& a_pt,
    const AlignedCylinderBase<typename VertexType::value_type>& a_cylinder);

template <class ScalarType, class HalfEdgeType>
bool ellipseContainedInFace(
    const AlignedCylinderBase<ScalarType>& a_aligned_cylinder,
    const PlaneBase<ScalarType>& a_face_plane, HalfEdgeType* const a_half_edge);

template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType = NoSurfaceOutput>
ReturnType orientAndApplyType3Correction(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    HalfEdgeType* a_start, HalfEdgeType* a_end, bool* a_requires_nudge,
    SurfaceOutputType* a_surface = nullptr);

// template <class ReturnType, class ScalarType, class HalfEdgeType,
//           class SurfaceOutputType>
// ReturnType orientAndApplyType3CorrectionWithGradients(
//     const AlignedParaboloid& a_paraboloid, HalfEdgeType* a_start,
//     HalfEdgeType* a_end, SurfaceOutputType* a_surface, bool* a_requires_nudge);

template <class ScalarType, class HalfEdgeType, class PtType>
Normal determineNudgeDirection(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    HalfEdgeType* a_current_edge, const PtType& a_inter_pt);

template <class ScalarType>
AlignedCylinderBase<Quad_t> nudgeCylinder(
    AlignedCylinderBase<ScalarType>& a_cylinder,
    const UnsignedIndex_t a_nudge_iter);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AligneCylinderType,
          class SurfaceOutputType>
ReturnType reformQuadraticIntersectionBases(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AligneCylinderType& a_aligned_cylinder,
    const UnsignedIndex_t a_nudge_iter, SurfaceOutputType* a_surface);

}  // namespace IRL

#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection.tpp"

#endif  // IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_INTERSECTION_H_
