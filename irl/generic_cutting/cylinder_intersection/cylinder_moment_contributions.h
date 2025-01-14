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

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsVC3SeriesOne(
    const ContainerType& a_weight);

// template <class ContainerType, class ScalarType>
// inline std::array<ContainerType, 12> coeffsV3andC3SeriesOne(
//     const ContainerType& a_weight);

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsVC3Exact(
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
    {1.33333333333333333333e-1,-6.66666666666666666667e-2},
    {3.80952380952380952381e-2,1.90476190476190476190e-2},
    {-3.17460317460317460317e-2,9.52380952380952380952e-3},
    {2.30880230880230880231e-2,-1.84704184704184704185e-2},
    {-1.55400155400155400155e-2,1.77600177600177600178e-2},
    {9.94560994560994560995e-3,-1.38528138528138528139e-2},
    {-6.14287673111202522967e-3,9.72622149092737328031e-3},
    {3.69496344728542870958e-3,-6.40460330862807642993e-3},
    {-2.17738917429319906100e-3,4.03806719596193280404e-3},
    {1.26225459379315887594e-3,-2.46713397877753780298e-3},
    {-7.22009627649686877039e-4,1.47178885636282324935e-3},
    {4.08409486347297627416e-4,-8.61699136029463125977e-4},
    {-2.28850143211847808466e-4,4.96931739545726669812e-4},
    {1.27202064713285134234e-4,-2.83024593987059423671e-4},
    {-7.02089318222677688955e-5,1.59518823037358386682e-4},
    {3.85146140282154617941e-5,-8.91122442221455782687e-5},
    {-2.10138856944486387829e-5,4.94010646150196069634e-5},
    {1.14102546757187178912e-5,-2.72044493057925221300e-5},
    {-6.16895882874223365661e-6,1.48936291722491069710e-5},
    {3.32232788818431188361e-6,-8.11165770102143680674e-6},
    {-1.78298263332558071087e-6,4.39755400472396191140e-6},
    {9.53814417827666276940e-7,-2.37416849654929975456e-6},
    {-5.08760009509191195956e-7,1.27698762386806990185e-6},
    {2.70644710940950815240e-7,-6.84522930441420215776e-7},
    {-1.43620424438004559974e-7,3.65802448568592810874e-7},
    {7.60404865388125961172e-8,-1.94929183746321178935e-7},
    {-4.01752368150406631307e-8,1.03604613165880724872e-7},
    {2.11846823130157922347e-8,-5.49340727565099164154e-8},
    {-1.11504855968391084420e-8,2.90632011685354826489e-8},
    {5.85904005088426059351e-9,-1.53445423913279324818e-8},
    {-3.07374254977158901905e-9,8.08603863945594156717e-9},
    {1.61012609251395707883e-9,-4.25348818022403634729e-9},
    {-8.42253458719936787160e-10,2.23374110564716008427e-9},
    {4.39999246040632790219e-10,-1.17123608827006537968e-9},
    {-2.29572531677123071448e-10,6.13227528308711627787e-10},
    {1.19640085079734994949e-10,-3.20628620611978279493e-10},
    {-6.22801741594724378468e-11,1.67425245513723476641e-10},
    {3.23865428403688352813e-11,-8.73191020426867443545e-11},
    {-1.68245880447530109210e-11,4.54879410917285673314e-11},
    {8.73194560246680826297e-12,-2.36705703091260516676e-11},
    {-4.52777061680852440224e-12,1.23047055300643619968e-11}};

// static double cx3Series[41][4] = {
//     {..., ..., ..., ...}, {...}, ...};

// static double cz3Series[41][5] = {
//     {..., ..., ..., ..., ...}, {...}, ...};
}  // namespace IRL

#include "irl/generic_cutting/cylinder_intersection/cylinder_moment_contributions.tpp"

#endif  // IRL_GENERIC_CUTTING_CYLINDER_INTERSECTION_CYLINDER_MOMENT_CONTRIBUTIONS_H_
