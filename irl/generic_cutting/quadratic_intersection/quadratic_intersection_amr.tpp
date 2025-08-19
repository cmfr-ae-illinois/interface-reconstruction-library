// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_AMR_TPP_
#define IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_AMR_TPP_

#include <float.h>
#include <cassert>
#include <cmath>
#include <fstream>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"

namespace IRL {

const double triangleSignedArea(const Pt& a_pt_0, const Pt& a_pt_1,
                                const Pt& a_pt_2) {
  return 0.5 * ((a_pt_1[0] - a_pt_0[0]) * (a_pt_2[1] - a_pt_0[1]) -
                (a_pt_2[0] - a_pt_0[0]) * (a_pt_1[1] - a_pt_0[1]));
}

template <>
void kahanSummationMoments(
    std::array<std::pair<Volume, Volume>, N_AMR_STRATEGIES>& a_full_moments,
    std::array<std::pair<Volume, Volume>, N_AMR_STRATEGIES>& a_full_moments_ref,
    std::array<Volume, N_AMR_STRATEGIES>& a_moments_to_add) {
  for (UnsignedIndex_t m = 0; m < N_AMR_STRATEGIES; ++m) {
    const double y_first = a_moments_to_add[m] - a_full_moments_ref[m].first;
    const double y_second =
        std::fabs(a_moments_to_add[m]) - a_full_moments_ref[m].second;
    const double t_first = a_full_moments[m].first + y_first;
    const double t_second = a_full_moments[m].first + y_second;
    a_full_moments_ref[m].first = (t_first - a_full_moments[m].first) - y_first;
    a_full_moments_ref[m].second =
        (t_second - a_full_moments[m].second) - y_second;
    a_full_moments[m].first = t_first;
    a_full_moments[m].second = t_second;
  }
}

template <>
void kahanSummationMoments(
    std::array<std::pair<VolumeMoments, VolumeMoments>, N_AMR_STRATEGIES>&
        a_full_moments,
    std::array<std::pair<VolumeMoments, VolumeMoments>, N_AMR_STRATEGIES>&
        a_full_moments_ref,
    std::array<VolumeMoments, N_AMR_STRATEGIES>& a_moments_to_add) {
  for (UnsignedIndex_t m = 0; m < N_AMR_STRATEGIES; ++m) {
    const double y_first_m0 =
        a_moments_to_add[m].volume() - a_full_moments_ref[m].first.volume();
    const double y_first_m1x = a_moments_to_add[m].centroid()[0] -
                               a_full_moments_ref[m].first.centroid()[0];
    const double y_first_m1y = a_moments_to_add[m].centroid()[1] -
                               a_full_moments_ref[m].first.centroid()[1];
    const double y_first_m1z = a_moments_to_add[m].centroid()[2] -
                               a_full_moments_ref[m].first.centroid()[2];
    const double y_second_m0 = std::fabs(a_moments_to_add[m].volume()) -
                               a_full_moments_ref[m].second.volume();
    const double y_second_m1x = std::fabs(a_moments_to_add[m].centroid()[0]) -
                                a_full_moments_ref[m].second.centroid()[0];
    const double y_second_m1y = std::fabs(a_moments_to_add[m].centroid()[1]) -
                                a_full_moments_ref[m].second.centroid()[1];
    const double y_second_m1z = std::fabs(a_moments_to_add[m].centroid()[2]) -
                                a_full_moments_ref[m].second.centroid()[2];
    const double t_first_m0 = a_full_moments[m].first.volume() + y_first_m0;
    const double t_first_m1x =
        a_full_moments[m].first.centroid()[0] + y_first_m1x;
    const double t_first_m1y =
        a_full_moments[m].first.centroid()[1] + y_first_m1y;
    const double t_first_m1z =
        a_full_moments[m].first.centroid()[2] + y_first_m1z;
    const double t_second_m0 = a_full_moments[m].first.volume() + y_second_m0;
    const double t_second_m1x =
        a_full_moments[m].first.centroid()[0] + y_second_m1x;
    const double t_second_m1y =
        a_full_moments[m].first.centroid()[1] + y_second_m1y;
    const double t_second_m1z =
        a_full_moments[m].first.centroid()[2] + y_second_m1z;
    a_full_moments_ref[m].first.volume() =
        (t_first_m0 - a_full_moments[m].first.volume()) - y_first_m0;
    a_full_moments_ref[m].first.centroid()[0] =
        (t_first_m1x - a_full_moments[m].first.centroid()[0]) - y_first_m1x;
    a_full_moments_ref[m].first.centroid()[1] =
        (t_first_m1y - a_full_moments[m].first.centroid()[1]) - y_first_m1y;
    a_full_moments_ref[m].first.centroid()[2] =
        (t_first_m1z - a_full_moments[m].first.centroid()[2]) - y_first_m1z;
    a_full_moments_ref[m].second.volume() =
        (t_second_m0 - a_full_moments[m].second.volume()) - y_second_m0;
    a_full_moments_ref[m].second.centroid()[0] =
        (t_second_m1x - a_full_moments[m].second.centroid()[0]) - y_second_m1x;
    a_full_moments_ref[m].second.centroid()[1] =
        (t_second_m1y - a_full_moments[m].second.centroid()[1]) - y_second_m1y;
    a_full_moments_ref[m].second.centroid()[2] =
        (t_second_m1z - a_full_moments[m].second.centroid()[2]) - y_second_m1z;
    a_full_moments[m].first.volume() = t_first_m0;
    a_full_moments[m].first.centroid()[0] = t_first_m1x;
    a_full_moments[m].first.centroid()[1] = t_first_m1y;
    a_full_moments[m].first.centroid()[2] = t_first_m1z;
    a_full_moments[m].second.volume() = t_second_m0;
    a_full_moments[m].second.centroid()[0] = t_second_m1x;
    a_full_moments[m].second.centroid()[1] = t_second_m1y;
    a_full_moments[m].second.centroid()[2] = t_second_m1z;
  }
}
}  // namespace IRL
#endif  // IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_AMR_TPP_
