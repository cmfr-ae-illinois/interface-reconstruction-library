// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This class is heavily inspired from the iMultiFab class in AmReX. Notable
// additions include the geometric shifting of the interface across periodic
// boundaries during parallel copy

#ifndef IRL_SEPUMULTIFAB_H_
#define IRL_SEPUMULTIFAB_H_

#include <AMReX_FabArray.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include "irl/amrex/sepunion_arraybox.h"

namespace amrex {

class SepUnionMultiFab : public FabArray<SepUnionArrayBox> {
 public:
  using FabArray<SepUnionArrayBox>::FabArray;

  /**
   * \brief Copy from src to dst including nghost ghost cells.
   * The two SepUnionMultiFabs MUST have the same underlying BoxArray.
   * The copy is local.
   */
  static void Copy(SepUnionMultiFab& dst, const SepUnionMultiFab& src,
                   int srccomp, int dstcomp, int numcomp, int nghost) {
    amrex::Copy(dst, src, srccomp, dstcomp, numcomp, IntVect(nghost));
  }

  static void Copy(SepUnionMultiFab& dst, const SepUnionMultiFab& src,
                   int srccomp, int dstcomp, int numcomp,
                   const IntVect& nghost) {
    // don't have to BL_ASSERT(dst.boxArray() == src.boxArray());
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost));

    BL_PROFILE("SepUnionMultiFab::Copy()");

    amrex::Copy(dst, src, srccomp, dstcomp, numcomp, nghost);
  }
  void ParallelCopyWithPeriodicShift(
      const SepUnionMultiFab& src, int scomp, int dcomp, int ncomp,
      const IntVect& snghost, const IntVect& dnghost, const Geometry& geom,
      CpOp op = FabArrayBase::COPY, const FabArrayBase::CPC* a_cpc = nullptr,
      bool deterministic = false);

  void ParallelCopyWithPeriodicShift(const SepUnionMultiFab& src, int src_comp,
                                     int dest_comp, int num_comp,
                                     int src_nghost, int dst_nghost,
                                     const Geometry& geom);

  void FillBoundaryWithPeriodicShift(int scomp, int ncomp, const Geometry& geom,
                                     bool cross);

  void FillBoundaryWithPeriodicShift(const Geometry& geom);

  void PeriodicShift(const Geometry& geom);
};

static_assert(sizeof(SepUnionMultiFab::value_type) ==
                  sizeof(IRL::SeparatorUnion),
              "SeparatorUnion size must match AMReX FabArray container");

static_assert(
    sizeof(IRL::SeparatorUnion) % sizeof(double) == 0,
    "SeparatorUnion must be a multiple of double size for AMReX FabArray");

// Write SepUnionMultiFab checkpoint
static void SepUnionMultiFab_Write(const SepUnionMultiFab& mf,
                                   const std::string& name,
                                   bool set_ghost = false);

// Read SepUnionMultiFab checkpoint
static void SepUnionMultiFab_Read(
    SepUnionMultiFab& mf, const std::string& name,
    const char* faHeader = nullptr,
    int coordinatorProc = ParallelDescriptor::IOProcessorNumber(),
    int allow_empty_mf = 0);

}  // namespace amrex

#include "irl/amrex/sepunion_multifab.tpp"

#endif  // IRL_SEPUMULTIFAB_H_
