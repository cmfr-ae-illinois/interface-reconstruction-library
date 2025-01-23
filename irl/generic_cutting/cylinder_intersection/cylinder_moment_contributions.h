// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_MOMENT_CONTRIBUTIONS_H_
#define IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_MOMENT_CONTRIBUTIONS_H_

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

template <class ReturnType, class ScalarType>
ReturnType computeType2Contribution(
    const AlignedCylinderBase<ScalarType>& a_aligned_cylinder,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1);

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 2> coeffsVC3SeriesOne(
    const ContainerType& a_weight);

// template <class ContainerType, class ScalarType>
// inline std::array<ContainerType, 12> coeffsV3andC3SeriesOne(
//     const ContainerType& a_weight);

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 2> coeffsVC3Exact(
    const ContainerType& a_weight);

// template <class ContainerType, class ScalarType>
// inline std::array<ContainerType, 12> coeffsV3andC3Exact(
//     const ContainerType& a_weight);

template <class ReturnType, class ScalarType>
ReturnType computeType3Contribution(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const RationalBezierArcBase<ScalarType>& a_arc);

// Don't know what these two function do sooooooo.

// template <class ReturnType, class ScalarType>
// ReturnType computeFaceOnlyContribution(
//     const AlignedParaboloidBase<ScalarType>& a_paraboloid,
//     const PlaneBase<ScalarType>& a_face_plane,
//     const PtBase<ScalarType>& a_pt_ref);

// template <class ReturnType, class ScalarType>
// ReturnType computeTriangleCorrection(
//     const AlignedParaboloidBase<ScalarType>& a_paraboloid,
//     const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1,
//     const PtBase<ScalarType>& a_pt_2);

// renaming this vc3 for volume cylinder 3
// to distinguish with v3 for paraboloid
static double vc3Series[41][2] = {
    {2.6666666666666666667e-1,1.3333333333333333333e-1},
    {7.6190476190476190476e-2,1.1428571428571428571e-1},
    {-6.3492063492063492063e-2,-4.4444444444444444444e-2},
    {4.6176046176046176046e-2,9.2352092352092352092e-3},
    {-3.1080031080031080031e-2,4.4400044400044400044e-3},
    {1.9891219891219891220e-2,-7.8144078144078144078e-3},
    {-1.2285753462224050459e-2,7.1666895196306961013e-3},
    {7.3899268945708574192e-3,-5.4192797226852954407e-3},
    {-4.3547783485863981220e-3,3.7213560433374674861e-3},
    {2.5245091875863177519e-3,-2.4097587699687578541e-3},
    {-1.4440192552993737541e-3,1.4995584574262727446e-3},
    {8.1681897269459525483e-4,-9.0657929936433099712e-4},
    {-4.5770028642369561693e-4,5.3616319266775772269e-4},
    {2.5440412942657026847e-4,-3.1164505854754857887e-4},
    {-1.4041786364453553779e-4,1.7861978243018123557e-4},
    {7.7029228056430923588e-5,-1.0119526038786023295e-4},
    {-4.2027771388897277566e-5,5.6774357841141936361e-5},
    {2.2820509351437435782e-5,-3.1588389260147608478e-5},
    {-1.2337917657484467313e-5,1.7449340687013746629e-5},
    {6.6446557763686237672e-6,-9.5786596256742498463e-6},
    {-3.5659652666511614217e-6,5.2291427427967624011e-6},
    {1.9076288356553325539e-6,-2.8407081574432669552e-6},
    {-1.0175200190183823919e-6,1.5364552287177574118e-6},
    {5.4128942188190163048e-7,-8.2775643900093880107e-7},
    {-2.8724084887600911995e-7,4.4436404826117650180e-7},
    {1.5208097307762519223e-7,-2.3777739441501716564e-7},
    {-8.0350473630081326261e-8,1.2685875270168012348e-7},
    {4.2369364626031584469e-8,-6.7498780886988248361e-8},
    {-2.2300971193678216884e-8,3.5825431143392748414e-8},
    {1.1718080101768521187e-8,-1.8971004680887343777e-8},
    {-6.1474850995431780381e-9,1.0024592179368705096e-8},
    {3.2202521850279141577e-9,-5.2867241754201585369e-9},
    {-1.6845069174398735743e-9,2.7829752938544465942e-9},
    {8.7999849208126558044e-10,-1.4624736844588651789e-9},
    {-4.5914506335424614290e-10,7.6730999326317711268e-10},
    {2.3928017015946998990e-10,-4.0197707106448656909e-10},
    {-1.2456034831894487569e-10,2.1029014270850207759e-10},
    {6.4773085680737670563e-11,-1.0986511840463581815e-10},
    {-3.3649176089506021842e-11,5.7326706093951112821e-11},
    {1.7463891204933616526e-11,-2.9877249413318486809e-11},
    {-9.0555412336170488045e-12,1.5553869826511675189e-11}};

// static double cx3Series[41][4] = {
//     {..., ..., ..., ...}, {...}, ...};

// static double cz3Series[41][5] = {
//     {..., ..., ..., ..., ...}, {...}, ...};
}  // namespace IRL

#include "irl/generic_cutting/cylinder_intersection/cylinder_moment_contributions.tpp"

#endif  // IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_MOMENT_CONTRIBUTIONS_H_
