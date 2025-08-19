// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_H_
#define IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_H_

#include <float.h>
#include <cassert>
#include <cmath>
#include "irl/data_structures/stack_vector.h"

namespace IRL {

template <class ReturnType, class ScalarType, class HalfEdgeType, class PtType>
ReturnType computeUnclippedSegmentType1Contribution(
    const PtType& a_ref_pt, const HalfEdgeType a_entry_half_edge,
    HalfEdgeType& a_exit_half_edge, const bool skip_first);

template <class PtType, class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType>
void placeSingleIntercept(const PtType& a_intersection_location,
                          HalfEdgeType* a_half_edge_with_intersection,
                          SegmentedHalfEdgePolytopeType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope);

template <class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType, class VertexType>
void placeDoubleIntercept(
    const StackVector<VertexType, 2>& a_intersection_location,
    HalfEdgeType* a_half_edge_with_intersection,
    SegmentedHalfEdgePolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope);

template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdgeXY(
    const std::array<ScalarType, 2>& a_test_pt,
    const PtBase<ScalarType>& a_vertex_0, const PtBase<ScalarType>& a_vertex_1);

template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdgeYZ(
    const std::array<ScalarType, 2>& a_test_pt,
    const PtBase<ScalarType>& a_vertex_0, const PtBase<ScalarType>& a_vertex_1);

template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdgeWithComponent(
    const PtBase<ScalarType>& a_test_pt, const PtBase<ScalarType>& a_vertex_0,
    const PtBase<ScalarType>& a_vertex_1, const UnsignedIndex_t a_index);

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, void>
resetPolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                HalfEdgePolytopeType* a_complete_polytope);

template <class ScalarType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class SurfaceOutputType>
void nudgePolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                     HalfEdgePolytopeType* a_complete_polytope,
                     const UnsignedIndex_t a_nudge_iter,
                     SurfaceOutputType* a_surface);

}  // namespace IRL

#include "irl/generic_cutting/quadratic_intersection/quadratic_intersection.tpp"

#endif  // IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_H_