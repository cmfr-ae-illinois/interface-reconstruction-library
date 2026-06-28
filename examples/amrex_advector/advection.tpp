// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_ADVECTION_H_
#define EXAMPLES_AMREX_ADVECTOR_ADVECTION_H_

#include "irl/amrex/sepunion_multifab.h"

#include "examples/amrex_advector/advection_remap.h"

using namespace amrex;

void AmrCoreAdv::TransportMoments(
    const SepUnionMultiFab& a_interface_with_ghost,
    const Array<MultiFab, AMREX_SPACEDIM>& a_facevel, const MultiFab& a_band_id,
    MultiFab& a_moments, const Geometry& a_geom, const double a_dt) {
  if (advection_name == "remap" || advection_name == "default") {
    LagrangianRemap::TransportMoments(a_interface_with_ghost, a_facevel,
                                      a_band_id, a_moments, a_geom, a_dt);
  } else {
    std::ostringstream oss;
    oss << "Unknown advection method: " << advection_name << '\n';
    throw std::runtime_error(oss.str());
  }
}

// void AmrCoreAdv::TransportMoments(
//     const SepUnionMultiFab& a_interface_with_ghost,
//     const Array<MultiFab, AMREX_SPACEDIM>& a_facevel, const MultiFab&
//     a_band_id, MultiFab& a_moments, const Geometry& a_geom, const double
//     a_dt, const double a_time) {
//   if (advection_name == "remap" || advection_name == "default") {
//     LagrangianRemap::TransportMoments(a_interface_with_ghost, a_facevel,
//                                       a_band_id, a_moments, a_geom, a_dt,
//                                       a_time);
//   } else {
//     std::ostringstream oss;
//     oss << "Unknown advection method: " << advection_name << '\n';
//     throw std::runtime_error(oss.str());
//   }
// }

#endif  // EXAMPLES_AMREX_ADVECTOR_ADVECTION_H_
