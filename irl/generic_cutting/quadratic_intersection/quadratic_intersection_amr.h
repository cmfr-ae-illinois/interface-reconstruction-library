// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_AMR_H_
#define IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_AMR_H_

#include <float.h>
#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/quadratic_intersection/moment_contributions.h"
#include "irl/geometry/general/geometry_type_traits.h"
#include "irl/quadratic_reconstruction/ellipse.h"

namespace IRL {

#ifndef AMR_DBL_EPSILON
#define AMR_DBL_EPSILON 10.0 * DBL_EPSILON
#endif

#ifndef no_amr_output
#define no_amr_output ""
#endif

#define ONE_AMR_STRATEGY

#ifdef ONE_AMR_STRATEGY
static constexpr UnsignedIndex_t N_AMR_STRATEGIES = 1;
#elif defined(TWO_AMR_STRATEGY)
static constexpr UnsignedIndex_t N_AMR_STRATEGIES = 2;
#elif defined(THREE_AMR_STRATEGY)
static constexpr UnsignedIndex_t N_AMR_STRATEGIES = 3;
#else
static constexpr UnsignedIndex_t N_AMR_STRATEGIES = 1;
#endif

template <class ReturnType>
inline void kahanSummationMoments(
    std::array<std::pair<ReturnType, ReturnType>, N_AMR_STRATEGIES>&
        a_full_moments,
    std::array<std::pair<ReturnType, ReturnType>, N_AMR_STRATEGIES>&
        a_full_moments_ref,
    std::array<ReturnType, N_AMR_STRATEGIES>& a_moments_to_add);

inline const double triangleSignedArea(const Pt& a_pt_0, const Pt& a_pt_1,
                                       const Pt& a_pt_2);

}  // namespace IRL

#include "irl/generic_cutting/quadratic_intersection/quadratic_intersection_amr.tpp"

#endif  // IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_AMR_H_
