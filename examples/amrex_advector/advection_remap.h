// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_ADVECTION_REMAP_H_
#define EXAMPLES_AMREX_ADVECTOR_ADVECTION_REMAP_H_

#include "irl/amrex/sepunion_multifab.h"

#include "examples/amrex_advector/advection_helpers.h"

using namespace amrex;

// Lookup tables for construction of flux-corrected Poly24
static constexpr std::array<std::array<int, 4>, 6> flux_id_table = {
    {{4, 5, 6, 7},
     {0, 1, 2, 3},
     {1, 5, 4, 0},
     {2, 6, 7, 3},
     {6, 5, 1, 2},
     {7, 4, 0, 3}}};
static constexpr std::array<int, 6> face_center_id_table = {
    {13, 8, 9, 11, 10, 12}};

struct LagrangianRemap {
  static void TransportMoments(const SepUnionMultiFab& a_interface_with_ghost,
                               const Array<MultiFab, AMREX_SPACEDIM>& a_facevel,
                               const MultiFab& a_band_id, MultiFab& a_moments,
                               const Geometry& a_geom, const double a_dt,
                               const double a_time,
                               const VelocityFieldType velocity_field_type) {
    const auto dx = a_geom.CellSizeArray();
    const auto problo = a_geom.ProbLoArray();
    const double new_time = a_time;                // t^{n+1}
    const double old_time = a_time - a_dt;         // t^n
    const double half_time = a_time - 0.5 * a_dt;  // t^{n+1/2}

    for (MFIter mfi(a_interface_with_ghost, TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      // Valid tile and tile including ghost cells.
      const Box& bx = mfi.tilebox();
      const Box& grown_bx = mfi.growntilebox();

      // only used when velocity field type is interpolated
      const Array4<Real const> velx = a_facevel[0].const_array(mfi);
      const Array4<Real const> vely = a_facevel[1].const_array(mfi);
      const Array4<Real const> velz = a_facevel[2].const_array(mfi);

      const Array4<IRL::SeparatorUnion const> interface_array =
          a_interface_with_ghost.const_array(mfi);
      const Array4<Real const> band_id_array = a_band_id.const_array(mfi);
      const Array4<Real> moments_array = a_moments.array(mfi);

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (band_id_array(i, j, k) == 0.0) {
          return;
        }
        std::array<double, 6> flux_volumes;
        std::array<IRL::Pt, 8> cell;
        std::array<IRL::Pt, 14> preimage;
        IRL::CappedDodecahedron flux;
        const double x = problo[0] + i * dx[0];
        const double y = problo[1] + j * dx[1];
        const double z = problo[2] + k * dx[2];
        // cell corners
        cell[0] = IRL::Pt(x + dx[0], y, z + dx[2]);
        cell[1] = IRL::Pt(x + dx[0], y, z);
        cell[2] = IRL::Pt(x + dx[0], y + dx[1], z);
        cell[3] = IRL::Pt(x + dx[0], y + dx[1], z + dx[2]);
        cell[4] = IRL::Pt(x, y, z + dx[2]);
        cell[5] = IRL::Pt(x, y, z);
        cell[6] = IRL::Pt(x, y + dx[1], z);
        cell[7] = IRL::Pt(x, y + dx[1], z + dx[2]);
        // projecting cell corners
        for (int n = 0; n < 8; ++n) {
          preimage[n] =
              ProjectVertex(cell[n], -a_dt, new_time, velocity_field_type, velx,
                            vely, velz, grown_bx, a_geom);
        }
        // face center locations
        const IRL::Pt xlo_face(x, y + 0.5 * dx[1], z + 0.5 * dx[2]);
        const IRL::Pt xhi_face(x + dx[0], y + 0.5 * dx[1], z + 0.5 * dx[2]);
        const IRL::Pt ylo_face(x + 0.5 * dx[0], y, z + 0.5 * dx[2]);
        const IRL::Pt yhi_face(x + 0.5 * dx[0], y + dx[1], z + 0.5 * dx[2]);
        const IRL::Pt zlo_face(x + 0.5 * dx[0], y + 0.5 * dx[1], z);
        const IRL::Pt zhi_face(x + 0.5 * dx[0], y + 0.5 * dx[1], z + dx[2]);
        // Evaluate face velocities at t^{n+1/2}
        const auto u_xlo = GetVelocity(xlo_face, half_time, velocity_field_type,
                                       velx, vely, velz, grown_bx, a_geom);
        const auto u_xhi = GetVelocity(xhi_face, half_time, velocity_field_type,
                                       velx, vely, velz, grown_bx, a_geom);
        const auto u_ylo = GetVelocity(ylo_face, half_time, velocity_field_type,
                                       velx, vely, velz, grown_bx, a_geom);
        const auto u_yhi = GetVelocity(yhi_face, half_time, velocity_field_type,
                                       velx, vely, velz, grown_bx, a_geom);
        const auto u_zlo = GetVelocity(zlo_face, half_time, velocity_field_type,
                                       velx, vely, velz, grown_bx, a_geom);
        const auto u_zhi = GetVelocity(zhi_face, half_time, velocity_field_type,
                                       velx, vely, velz, grown_bx, a_geom);
        // Compute soleinoidal flux volumes
        flux_volumes[0] = a_dt * u_xlo[0] * dx[1] * dx[2];
        flux_volumes[1] = a_dt * u_xhi[0] * dx[1] * dx[2];
        flux_volumes[2] = a_dt * u_ylo[1] * dx[0] * dx[2];
        flux_volumes[3] = a_dt * u_yhi[1] * dx[0] * dx[2];
        flux_volumes[4] = a_dt * u_zlo[2] * dx[0] * dx[1];
        flux_volumes[5] = a_dt * u_zhi[2] * dx[0] * dx[1];
        // Create face flux hexahedra to compute correction
        for (int f = 0; f < 6; ++f) {
          for (int n = 0; n < 4; ++n) {
            flux[n] = cell[flux_id_table[f][n]];
            flux[n + 4] = preimage[flux_id_table[f][n]];
          }
          // Back-project the destination face center from t^{n+1} to t^n.
          flux[8] = ProjectVertex(
              0.25 * (flux[0] + flux[1] + flux[2] + flux[3]), -a_dt, new_time,
              velocity_field_type, velx, vely, velz, grown_bx, a_geom);
          flux.adjustCapToMatchVolume(flux_volumes[f]);
          preimage[face_center_id_table[f]] = flux[8];
        }
        const auto preimage_cell =
            IRL::Polyhedron24::fromRawPtPointer(14, preimage.data());

        // Compute bounding box of the preimage polyhedron
        Real xlo = preimage[0][0];
        Real ylo = preimage[0][1];
        Real zlo = preimage[0][2];
        Real xhi = preimage[0][0];
        Real yhi = preimage[0][1];
        Real zhi = preimage[0][2];
        for (int n = 1; n < 14; ++n) {
          xlo = preimage[n][0] < xlo ? preimage[n][0] : xlo;
          ylo = preimage[n][1] < ylo ? preimage[n][1] : ylo;
          zlo = preimage[n][2] < zlo ? preimage[n][2] : zlo;
          xhi = preimage[n][0] > xhi ? preimage[n][0] : xhi;
          yhi = preimage[n][1] > yhi ? preimage[n][1] : yhi;
          zhi = preimage[n][2] > zhi ? preimage[n][2] : zhi;
        }
        const int ilo =
            static_cast<int>(amrex::Math::floor((xlo - problo[0]) / dx[0]));
        const int jlo =
            static_cast<int>(amrex::Math::floor((ylo - problo[1]) / dx[1]));
        const int klo =
            static_cast<int>(amrex::Math::floor((zlo - problo[2]) / dx[2]));
        const int ihi =
            static_cast<int>(amrex::Math::floor((xhi - problo[0]) / dx[0]));
        const int jhi =
            static_cast<int>(amrex::Math::floor((yhi - problo[1]) / dx[1]));
        const int khi =
            static_cast<int>(amrex::Math::floor((zhi - problo[2]) / dx[2]));

        if (!grown_bx.contains(ilo, jlo, klo)) {
          std::ostringstream oss;
          oss << "Cell " << i << " " << j << " " << k << " has lo " << ilo
              << " " << jlo << " " << klo << '\n';
          throw std::runtime_error(oss.str());
        }
        if (!grown_bx.contains(ihi, jhi, khi)) {
          std::ostringstream oss;
          oss << "Cell " << i << " " << j << " " << k << " has hi " << ihi
              << " " << jhi << " " << khi << '\n';
          throw std::runtime_error(oss.str());
        }

        // Intersect preimage
        for (int n = 0; n < 4; ++n) {
          moments_array(i, j, k, n) = 0.0;
        }
        for (int ii = ilo; ii <= ihi; ++ii) {
          for (int jj = jlo; jj <= jhi; ++jj) {
            for (int kk = klo; kk <= khi; ++kk) {
              const double xloc = problo[0] + ii * dx[0];
              const double yloc = problo[1] + jj * dx[1];
              const double zloc = problo[2] + kk * dx[2];
              if (xlo > xloc + dx[0] || xhi < xloc || ylo > yloc + dx[1] ||
                  yhi < yloc || zlo > zloc + dx[2] || zhi < zloc) {
                continue;
              }
              const auto cell_loc = IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(xloc, yloc, zloc),
                  IRL::Pt(xloc + dx[0], yloc + dx[1], zloc + dx[2]));
              IRL::PlanarLocalizer localizer = cell_loc.getLocalizer();
              IRL::LocalizedSeparatorUnion local_sep(
                  &localizer, &interface_array(ii, jj, kk));
              const auto cut_moments =
                  IRL::getVolumeMoments<IRL::VolumeMoments>(preimage_cell,
                                                            local_sep);
              const double m0 = cut_moments.volume();
              IRL::Pt centroid =
                  (1.0 / IRL::safelyEpsilon(m0)) * cut_moments.centroid();
              centroid[0] = std::max(centroid[0], xloc);
              centroid[0] = std::min(centroid[0], xloc + dx[0]);
              centroid[1] = std::max(centroid[1], yloc);
              centroid[1] = std::min(centroid[1], yloc + dx[1]);
              centroid[2] = std::max(centroid[2], zloc);
              centroid[2] = std::min(centroid[2], zloc + dx[2]);
              // forward projecting the centroid
              centroid =
                  ProjectVertex(centroid, a_dt, old_time, velocity_field_type,
                                velx, vely, velz, grown_bx, a_geom);
              moments_array(i, j, k, 0) += m0;
              moments_array(i, j, k, 1) += m0 * centroid[0];
              moments_array(i, j, k, 2) += m0 * centroid[1];
              moments_array(i, j, k, 3) += m0 * centroid[2];
            }
          }
        }
      });
    }
    a_moments.FillBoundary(a_geom.periodicity());
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_ADVECTION_REMAP_H_
