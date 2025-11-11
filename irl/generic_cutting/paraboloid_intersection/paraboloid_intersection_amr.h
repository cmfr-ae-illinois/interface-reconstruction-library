// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_AMR_H_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_AMR_H_

#include <float.h>
#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/quadratic_intersection/moment_contributions.h"
#include "irl/generic_cutting/quadratic_intersection/quadratic_intersection_amr.h"
#include "irl/geometry/general/geometry_type_traits.h"
#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/quadratic_reconstruction/ellipse.h"

namespace IRL {

template <class ReturnType, UnsignedIndex_t strategyID>
inline ReturnType computeMomentContributionClippedTriangle(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    std::vector<double>& a_amr_triangles_clipped,
    std::vector<double>& a_amr_triangles_unclipped, const bool a_print);

template <class ReturnType, UnsignedIndex_t strategyID>
inline ReturnType computeMomentContributionUnclippedTriangle(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    std::vector<double>& a_amr_triangles_clipped,
    std::vector<double>& a_amr_triangles_unclipped, const bool a_print);

template <class ReturnType>
inline void computeMomentContributionMixedTriangle(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const Normal& a_normal,
    const double a_signed_area,
    std::array<ReturnType, N_AMR_STRATEGIES>& a_moments_to_add,
    std::vector<double>& a_amr_triangles_clipped,
    std::vector<double>& a_amr_triangles_unclipped, const bool a_print);

inline std::pair<bool, bool> computeZBounds(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2);

template <class ReturnType>
inline void computeMomentContributionAMR(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const Normal& a_normal,
    const double a_signed_area, const UnsignedIndex_t a_amr_level,
    const UnsignedIndex_t a_max_amr_level,
    std::array<std::pair<ReturnType, ReturnType>, N_AMR_STRATEGIES>&
        a_full_moments_ref,
    std::array<std::pair<ReturnType, ReturnType>, N_AMR_STRATEGIES>&
        a_full_moments,
    std::vector<double>& a_amr_triangles_clipped,
    std::vector<double>& a_amr_triangles_unclipped, const bool a_print);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithParaboloidAMR(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AlignedParaboloid& a_paraboloid,
    const UnsignedIndex_t a_max_amr_level,
    const std::string& a_filename = no_amr_output);
}  // namespace IRL

#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection_amr.tpp"

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_AMR_H_
