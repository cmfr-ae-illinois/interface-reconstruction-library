// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_INTERSECTION_AMR_H_
#define IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_INTERSECTION_AMR_H_

#include <float.h>
#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/quadratic_intersection/moment_contributions.h"
#include "irl/generic_cutting/quadratic_intersection/quadratic_intersection_amr.h"
#include "irl/geometry/general/geometry_type_traits.h"
#include "irl/cylinder_reconstruction/aligned_cylinder.h"
#include "irl/quadratic_reconstruction/ellipse.h"

namespace IRL {

template <class ReturnType>
ReturnType computeMomentContributionClippedTriangle(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print);

template <class ReturnType>
ReturnType computeMomentContributionUnclippedTriangle(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print);

template <class ReturnType>
void computeMomentContributionMixedTriangle(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const Normal& a_normal,
    const double a_signed_area, std::array<ReturnType, 1>& a_moments_to_add,
    const bool a_print);

std::pair<bool, bool> computeZBounds(const AlignedCylinder& a_aligned_cylinder,
                                     const Pt& a_pt_0, const Pt& a_pt_1,
                                     const Pt& a_pt_2);

template <class ReturnType>
void computeMomentContributionAMR(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const Normal& a_normal,
    const double a_signed_area,
    const UnsignedIndex_t a_max_amr_level, const double max_length,
    std::array<std::pair<ReturnType, ReturnType>, 1>& a_full_moments_ref,
    std::array<std::pair<ReturnType, ReturnType>, 1>& a_full_moments,
    const bool a_print);

template <class SegmentedHalfEdgePolyhedronType>
void printClippedPolytope(SegmentedHalfEdgePolyhedronType* a_polytope);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithCylinderAMR(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AlignedCylinder& a_cylinder, const UnsignedIndex_t a_max_amr_level,
    const std::string& a_filename = no_amr_output);
}  // namespace IRL

#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection_amr.tpp"

#endif  // IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_INTERSECTION_AMR_H_
