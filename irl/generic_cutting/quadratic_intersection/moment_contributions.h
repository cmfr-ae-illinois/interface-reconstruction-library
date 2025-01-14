// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_MOMENT_CONTRIBUTIONS_H_
#define IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_MOMENT_CONTRIBUTIONS_H_

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
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/mymath.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/quadratic_reconstruction/rational_bezier_arc.h"

namespace IRL {

// The definition seems to be commented, soooo useless???
// template <class ReturnType, class ScalarType>
// ReturnType computeType1ContributionQuadrature(
//     const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1,
//     const NormalBase<ScalarType>& a_face_normal,
//     const UnsignedIndex_t a_max_normal_index);

template <class ReturnType, class ScalarType>
ReturnType computeType1Contribution(const PtBase<ScalarType>& a_ref_pt,
                                    const PtBase<ScalarType>& a_pt_0,
                                    const PtBase<ScalarType>& a_pt_1);
}  // namespace IRL

#include "irl/generic_cutting/quadratic_intersection/moment_contributions.tpp"

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_H_
