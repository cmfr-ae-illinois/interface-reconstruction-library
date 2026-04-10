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

#include "examples/amrex_advector/reconstruction_elvira.h"
#include "examples/amrex_advector/reconstruction_jibben.h"
#include "examples/amrex_advector/reconstruction_lvira.h"
#include "examples/amrex_advector/reconstruction_mof1.h"
#include "examples/amrex_advector/reconstruction_plicnet.h"

using namespace amrex;

void AmrCoreAdv::GetReconstruction(const int lev) {
  SepUnionMultiFab interface_with_ghost(grids[lev], dmap[lev],
                                        interface[lev].nComp(), num_grow);
  MultiFab moments_with_ghost(grids[lev], dmap[lev], ncomp_moments, num_grow);
  MultiFab::Copy(moments_with_ghost, moments_new[lev], 0, 0,
                 moments_new[lev].nComp(), moments_new[lev].nGrow());
  moments_with_ghost.FillBoundary(geom[lev].periodicity());
  GetReconstruction(interface[lev], interface_with_ghost, moments_with_ghost,
                    geom[lev]);
}

void AmrCoreAdv::GetReconstruction(SepUnionMultiFab& a_interface,
                                   SepUnionMultiFab& a_interface_with_ghost,
                                   const MultiFab& a_moments,
                                   const Geometry& a_geom) {
  if (reconstruction_name == "elvira" || reconstruction_name == "default") {
    ELVIRA::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                              a_geom);
  } else if (reconstruction_name == "lvira") {
    LVIRA::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                             a_geom);
  } else if (reconstruction_name == "plicnet") {
    PLICNet::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                               a_geom);
  } else if (reconstruction_name == "mof" || reconstruction_name == "mof1") {
    MOF1::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                            a_geom);
  } else if (reconstruction_name == "jibben") {
    Jibben::GetReconstruction(a_interface, a_interface_with_ghost, a_moments,
                              a_geom);
  } else {
    std::ostringstream oss;
    oss << "Unknown reconstruction method: " << reconstruction_name << '\n';
    throw std::runtime_error(oss.str());
  }
}

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTIONS_H_
