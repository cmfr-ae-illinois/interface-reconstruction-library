// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTIONS_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTIONS_H_

#include "irl/amrex/sepunion_multifab.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"

#include "examples/amrex_advector/reconstruction_cf.h"
#include "examples/amrex_advector/reconstruction_elvira.h"
#include "examples/amrex_advector/reconstruction_hybrid.h"
#include "examples/amrex_advector/reconstruction_ivf.h"
#include "examples/amrex_advector/reconstruction_lvira.h"
#include "examples/amrex_advector/reconstruction_mof1.h"
#include "examples/amrex_advector/reconstruction_mof2.h"
#include "examples/amrex_advector/reconstruction_plicnet.h"
#include "examples/amrex_advector/reconstruction_pu.h"
#include "examples/amrex_advector/reconstruction_vf.h"
#include "examples/amrex_advector/reconstruction_vf2.h"

using namespace amrex;

void AmrCoreAdv::GetReconstruction(const int lev) {
  SepUnionMultiFab interface_with_ghost(grids[lev], dmap[lev],
                                        interface[lev].nComp(), num_grow);
  InitializeSepUnionMultiFab(interface_with_ghost);
  MultiFab moments_with_ghost(grids[lev], dmap[lev], ncomp_moments, num_grow);
  moments_with_ghost.setVal(0.0);
  MultiFab::Copy(moments_with_ghost, moments_new[lev], 0, 0,
                 moments_new[lev].nComp(), 0);
  moments_with_ghost.FillBoundary(geom[lev].periodicity());
  GetReconstruction(interface[lev], interface_with_ghost, moments_with_ghost,
                    geom[lev], &interface_scalar_fields[lev]);
}

void AmrCoreAdv::GetReconstruction(
    SepUnionMultiFab& a_interface, SepUnionMultiFab& a_interface_with_ghost,
    const MultiFab& a_moments, const Geometry& a_geom,
    std::vector<InterfaceScalarField>* scalar_fields) {
  RecordNumberMixedCells();
  ReconstructionLoopTimer loop_timer(&reconstruction_time);
  loop_timer.start();
  if (reconstruction_name == "elvira" || reconstruction_name == "default") {
    ELVIRA::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                              a_geom, scalar_fields, &reconstruction_loop_time);
  } else if (reconstruction_name == "lvira") {
    LVIRA::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                             a_geom, scalar_fields, &reconstruction_loop_time);
  } else if (reconstruction_name == "plicnet") {
    PLICNet::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                               a_geom, scalar_fields,
                               &reconstruction_loop_time);
  } else if (reconstruction_name == "mof" || reconstruction_name == "mof1") {
    MOF1::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                            a_geom, scalar_fields, &reconstruction_loop_time);
  } else if (reconstruction_name == "vf") {
    VF::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                          a_geom, scalar_fields, &reconstruction_loop_time);
  } else if (reconstruction_name == "vf2") {
    VF2::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                           a_geom, scalar_fields, &reconstruction_loop_time);
  } else if (reconstruction_name == "ivf") {
    iVF::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                           a_geom, scalar_fields, &reconstruction_loop_time);
  } else if (reconstruction_name == "pu") {
    PU::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                          a_geom, scalar_fields, &reconstruction_loop_time);
  } else if (reconstruction_name == "cf") {
    CF::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                          a_geom, scalar_fields, &reconstruction_loop_time);
  } else if (reconstruction_name == "mof2") {
    number_mof2_iterations += MOF2::GetReconstruction(
        a_interface, a_interface_with_ghost, a_moments, a_geom, scalar_fields,
        &reconstruction_loop_time);
  } else if (reconstruction_name == "supermof2") {
    number_mof2_iterations += SuperMOF2::GetReconstruction(
        a_interface, a_interface_with_ghost, a_moments, a_geom, scalar_fields,
        &reconstruction_loop_time);
  } else if (reconstruction_name == "hybrid") {
    HYBRID::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                              a_geom, scalar_fields, &reconstruction_loop_time);
  } else {
    std::ostringstream oss;
    oss << "Unknown reconstruction method: " << reconstruction_name << '\n';
    throw std::runtime_error(oss.str());
  }
  loop_timer.stop();
}

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTIONS_H_
