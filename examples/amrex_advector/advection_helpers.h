// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_H_
#define EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_H_

#include "irl/amrex/sepunion_multifab.h"

using namespace amrex;

inline IRL::Vec3<double> GetVelocity(const IRL::Pt& pt,
                                     Array4<Real const> const& vx,
                                     Array4<Real const> const& vy,
                                     Array4<Real const> const& vz,
                                     const Box& bx, const Geometry& a_geom);

inline IRL::Pt ProjectVertex(const IRL::Pt& pt, const double dt,
                             Array4<Real const> const& vx,
                             Array4<Real const> const& vy,
                             Array4<Real const> const& vz, const Box& bx,
                             const Geometry& a_geom);

// inline IRL::Pt ProjectVertex(const IRL::Pt& pt, const double dt,
//                              const double time);

// inline IRL::Vec3<double> GetDeformationTestCaseVelocity(const IRL::Pt& pt,
//                                                         const double time);

#include "examples/amrex_advector/advection_helpers.tpp"

#endif  // EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_H_
