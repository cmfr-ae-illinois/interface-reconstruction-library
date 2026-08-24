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

void AmrCoreAdv::ResetMoments(const SepUnionMultiFab& a_interface,
                              MultiFab& a_moments, const Geometry& a_geom) {
  if (!transport_m1 && !transport_m2) return;

  const auto dx = a_geom.CellSizeArray();
  const auto problo = a_geom.ProbLoArray();
  const Real vol = dx[0] * dx[1] * dx[2];
  const Real vol_inv = 1.0 / vol;

  for (MFIter mfi(a_interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    Array4<IRL::SeparatorUnion const> interface_array =
        a_interface.const_array(mfi);
    Array4<Real> moments_array = a_moments.array(mfi);
    const Box& bx = mfi.tilebox();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      // Compute cell volume fraction
      const double vfrac = moments_array(i, j, k) * vol_inv;
      // Skip cell if vfrac ~0 or vfrac~1
      if (vfrac < IRL::global_constants::VF_LOW ||
          vfrac > IRL::global_constants::VF_HIGH) {
        return;
      }
      // Construct local cell as IRL object
      const double x = problo[0] + i * dx[0];
      const double y = problo[1] + j * dx[1];
      const double z = problo[2] + k * dx[2];
      const IRL::Pt lower_cell_pt(x, y, z);
      const IRL::Pt upper_cell_pt(x + dx[0], y + dx[1], z + dx[2]);
      const auto cell =
          IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt, upper_cell_pt);
      auto moments_exact =
          IRL::getVolumeMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
              cell, interface_array(i, j, k));
      moments_array(i, j, k, 1) = moments_exact[0].centroid()[0];
      moments_array(i, j, k, 2) = moments_exact[0].centroid()[1];
      moments_array(i, j, k, 3) = moments_exact[0].centroid()[2];
      moments_array(i, j, k, 4) = moments_exact[1].centroid()[0];
      moments_array(i, j, k, 5) = moments_exact[1].centroid()[1];
      moments_array(i, j, k, 6) = moments_exact[1].centroid()[2];
      if (transport_m2) {
        auto moments_full = IRL::getVolumeMoments<
            IRL::SeparatedMoments<IRL::GeneralMoments3D<2>>>(
            cell, interface_array(i, j, k));
        moments_array(i, j, k, 7) = moments_full[0][4];
        moments_array(i, j, k, 8) = moments_full[0][5];
        moments_array(i, j, k, 9) = moments_full[0][6];
        moments_array(i, j, k, 10) = moments_full[0][7];
        moments_array(i, j, k, 11) = moments_full[0][8];
        moments_array(i, j, k, 12) = moments_full[0][9];
        moments_array(i, j, k, 13) = moments_full[1][4];
        moments_array(i, j, k, 14) = moments_full[1][5];
        moments_array(i, j, k, 15) = moments_full[1][6];
        moments_array(i, j, k, 16) = moments_full[1][7];
        moments_array(i, j, k, 17) = moments_full[1][8];
        moments_array(i, j, k, 18) = moments_full[1][9];
      }
    });
  }
}

void AmrCoreAdv::TransportMoments(
    const SepUnionMultiFab& a_interface_with_ghost,
    const Array<MultiFab, AMREX_SPACEDIM>& a_facevel, const MultiFab& a_band_id,
    MultiFab& a_moments, const Geometry& a_geom, const double a_dt,
    const double a_time) {
  if (reset_moments) {
    ResetMoments(a_interface_with_ghost, a_moments, a_geom);
  }
  if (advection_name == "remap" || advection_name == "default") {
    LagrangianRemap::TransportMoments(
        a_interface_with_ghost, a_facevel, a_band_id, a_moments, a_geom, a_dt,
        a_time, velocity_field_type, transport_m1, transport_m2);
  } else {
    std::ostringstream oss;
    oss << "Unknown advection method: " << advection_name << '\n';
    throw std::runtime_error(oss.str());
  }
}

#endif  // EXAMPLES_AMREX_ADVECTOR_ADVECTION_H_
