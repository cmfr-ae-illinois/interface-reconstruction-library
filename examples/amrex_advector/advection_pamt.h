// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_ADVECTION_PAMT_H_
#define EXAMPLES_AMREX_ADVECTOR_ADVECTION_PAMT_H_

#include <array>
#include <cmath>
#include <sstream>
#include <utility>

#include "examples/amrex_advector/advection_helpers.h"
#include "irl/amrex/sepunion_multifab.h"

using namespace amrex;

static constexpr std::array<std::array<int, 4>, 6> pamt_flux_id_table = {{
    {{4, 5, 6, 7}},
    {{0, 1, 2, 3}},
    {{1, 5, 4, 0}},
    {{2, 6, 7, 3}},
    {{6, 5, 1, 2}},
    {{7, 4, 0, 3}},
}};

static constexpr std::array<int, 6> pamt_face_center_id_table = {
    {13, 8, 9, 11, 10, 12}};

static constexpr int pamt_cell_center_id = 14;

struct PAMT {
  static Eigen::Matrix3d MomentMatrix(
      const IRL::GeneralMoments3D<2>& a_moments) {
    return Eigen::Matrix3d({{a_moments[4], a_moments[5], a_moments[6]},
                            {a_moments[5], a_moments[7], a_moments[8]},
                            {a_moments[6], a_moments[8], a_moments[9]}});
  }

  static Eigen::Matrix3d AffineMatrix(
      const std::array<IRL::Pt, 4>& a_source,
      const std::array<IRL::Pt, 4>& a_destination, double& a_jacobian) {
    const Eigen::Vector3d x0(a_source[0][0], a_source[0][1], a_source[0][2]);
    const Eigen::Vector3d x1(a_source[1][0], a_source[1][1], a_source[1][2]);
    const Eigen::Vector3d x2(a_source[2][0], a_source[2][1], a_source[2][2]);
    const Eigen::Vector3d x3(a_source[3][0], a_source[3][1], a_source[3][2]);
    const Eigen::Vector3d y0(a_destination[0][0], a_destination[0][1],
                             a_destination[0][2]);
    const Eigen::Vector3d y1(a_destination[1][0], a_destination[1][1],
                             a_destination[1][2]);
    const Eigen::Vector3d y2(a_destination[2][0], a_destination[2][1],
                             a_destination[2][2]);
    const Eigen::Vector3d y3(a_destination[3][0], a_destination[3][1],
                             a_destination[3][2]);

    const Eigen::Vector3d e1 = x1 - x0;
    const Eigen::Vector3d e2 = x2 - x0;
    const Eigen::Vector3d e3 = x3 - x0;
    const Eigen::Vector3d f1 = y1 - y0;
    const Eigen::Vector3d f2 = y2 - y0;
    const Eigen::Vector3d f3 = y3 - y0;

    const double source_det = e1.dot(e2.cross(e3));
    const double destination_det = f1.dot(f2.cross(f3));
    a_jacobian = std::abs(destination_det / source_det);

    return (f1 * e2.cross(e3).transpose() + f2 * e3.cross(e1).transpose() +
            f3 * e1.cross(e2).transpose()) /
           source_det;
  }

  static double SignedTetDeterminant(const std::array<IRL::Pt, 4>& a_tet) {
    const Eigen::Vector3d x0(a_tet[0][0], a_tet[0][1], a_tet[0][2]);
    const Eigen::Vector3d x1(a_tet[1][0], a_tet[1][1], a_tet[1][2]);
    const Eigen::Vector3d x2(a_tet[2][0], a_tet[2][1], a_tet[2][2]);
    const Eigen::Vector3d x3(a_tet[3][0], a_tet[3][1], a_tet[3][2]);
    return (x0 - x3).dot((x1 - x3).cross(x2 - x3));
  }

  static void MakeSourceTetPositive(std::array<IRL::Pt, 4>& a_source,
                                    std::array<IRL::Pt, 4>& a_destination) {
    if (SignedTetDeterminant(a_source) < 0.0) {
      std::swap(a_source[1], a_source[2]);
      std::swap(a_destination[1], a_destination[2]);
    }
  }

  static void AddMappedVolumeMoments(const IRL::VolumeMoments& a_source_moments,
                                     const Eigen::Matrix3d& a_affine_matrix,
                                     const Eigen::Vector3d& a_translation,
                                     const double /*a_jacobian*/, double& a_m0,
                                     Eigen::Vector3d& a_m1) {
    const double source_m0 = a_source_moments.volume();
    const IRL::Pt& source_centroid = a_source_moments.centroid();
    const Eigen::Vector3d source_m1(source_centroid[0], source_centroid[1],
                                    source_centroid[2]);

    a_m0 += source_m0;
    a_m1 += a_affine_matrix * source_m1 + a_translation * source_m0;
  }

  static void AddMappedSecondMoment(
      const IRL::VolumeMoments& a_source_volume_moments,
      const IRL::GeneralMoments3D<2>& a_source_general_moments,
      const Eigen::Matrix3d& a_affine_matrix,
      const Eigen::Vector3d& a_translation, const double /*a_jacobian*/,
      Eigen::Matrix3d& a_m2) {
    const double source_m0 = a_source_volume_moments.volume();
    const IRL::Pt& source_centroid = a_source_volume_moments.centroid();
    const Eigen::Vector3d source_m1(source_centroid[0], source_centroid[1],
                                    source_centroid[2]);
    const Eigen::Matrix3d source_m2 = MomentMatrix(a_source_general_moments);

    a_m2 +=
        a_affine_matrix * source_m2 * a_affine_matrix.transpose() +
        (a_affine_matrix * source_m1) * a_translation.transpose() +
        a_translation * source_m1.transpose() * a_affine_matrix.transpose() +
        source_m0 * a_translation * a_translation.transpose();
  }

  static void TransportMoments(const SepUnionMultiFab& a_interface_with_ghost,
                               const Array<MultiFab, AMREX_SPACEDIM>& a_facevel,
                               const MultiFab& a_band_id, MultiFab& a_moments,
                               const Geometry& a_geom, const double a_dt,
                               const double a_time,
                               const VelocityFieldType velocity_field_type,
                               const bool transport_m1,
                               const bool transport_m2) {
    const auto dx = a_geom.CellSizeArray();
    const auto problo = a_geom.ProbLoArray();
    const double new_time = a_time;
    const double cell_volume = dx[0] * dx[1] * dx[2];
    const int ncomp = a_moments.nComp();
    const int ngrow = a_interface_with_ghost.nGrow();

    for (MFIter mfi(a_interface_with_ghost, TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const Box& bx = mfi.tilebox();
      const Box& grown_bx = grow(bx, ngrow);
      const Array4<Real const> velx = a_facevel[0].const_array(mfi);
      const Array4<Real const> vely = a_facevel[1].const_array(mfi);
      const Array4<Real const> velz = a_facevel[2].const_array(mfi);
      const Array4<IRL::SeparatorUnion const> interface_array =
          a_interface_with_ghost.const_array(mfi);
      const Array4<Real const> band_id_array = a_band_id.const_array(mfi);
      const Array4<Real> moments_array = a_moments.array(mfi);

      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (band_id_array(i, j, k) == 0.0) {
          return;
        }

        std::array<IRL::Pt, 15> cell;
        std::array<IRL::Pt, 15> preimage;
        std::array<double, 6> flux_volumes;
        IRL::CappedDodecahedron flux;
        const double x = problo[0] + i * dx[0];
        const double y = problo[1] + j * dx[1];
        const double z = problo[2] + k * dx[2];

        cell[0] = IRL::Pt(x + dx[0], y, z + dx[2]);
        cell[1] = IRL::Pt(x + dx[0], y, z);
        cell[2] = IRL::Pt(x + dx[0], y + dx[1], z);
        cell[3] = IRL::Pt(x + dx[0], y + dx[1], z + dx[2]);
        cell[4] = IRL::Pt(x, y, z + dx[2]);
        cell[5] = IRL::Pt(x, y, z);
        cell[6] = IRL::Pt(x, y + dx[1], z);
        cell[7] = IRL::Pt(x, y + dx[1], z + dx[2]);

        for (int n = 0; n < 8; ++n) {
          preimage[n] =
              ProjectVertex(cell[n], -a_dt, new_time, velocity_field_type, velx,
                            vely, velz, grown_bx, a_geom);
        }

        flux_volumes[0] = a_dt * velx(i, j, k) * dx[1] * dx[2];
        flux_volumes[1] = a_dt * velx(i + 1, j, k) * dx[1] * dx[2];
        flux_volumes[2] = a_dt * vely(i, j, k) * dx[0] * dx[2];
        flux_volumes[3] = a_dt * vely(i, j + 1, k) * dx[0] * dx[2];
        flux_volumes[4] = a_dt * velz(i, j, k) * dx[0] * dx[1];
        flux_volumes[5] = a_dt * velz(i, j, k + 1) * dx[0] * dx[1];

        // Correct each preimage face center to match the Eulerian face flux
        for (int f = 0; f < 6; ++f) {
          const int face_center_id = pamt_face_center_id_table[f];
          for (int n = 0; n < 4; ++n) {
            const int vertex_id = pamt_flux_id_table[f][n];
            flux[n] = cell[vertex_id];
            flux[n + 4] = preimage[vertex_id];
          }
          cell[face_center_id] = 0.25 * (flux[0] + flux[1] + flux[2] + flux[3]);
          flux[8] = ProjectVertex(cell[face_center_id], -a_dt, new_time,
                                  velocity_field_type, velx, vely, velz,
                                  grown_bx, a_geom);
          flux.adjustCapToMatchVolume(flux_volumes[f]);
          preimage[face_center_id] = flux[8];
        }

        // project the cell center from t^{n+1} to t^n
        cell[pamt_cell_center_id] =
            IRL::Pt(x + 0.5 * dx[0], y + 0.5 * dx[1], z + 0.5 * dx[2]);
        preimage[pamt_cell_center_id] = ProjectVertex(
            cell[pamt_cell_center_id], -a_dt, new_time, velocity_field_type,
            velx, vely, velz, grown_bx, a_geom);

        for (int n = 0; n < ncomp; ++n) {
          moments_array(i, j, k, n) = 0.0;
        }

        double M0_l = 0.0;
        double M0_g = 0.0;
        Eigen::Vector3d M1_l = Eigen::Vector3d::Zero();
        Eigen::Vector3d M1_g = Eigen::Vector3d::Zero();
        Eigen::Matrix3d M2_l = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d M2_g = Eigen::Matrix3d::Zero();

        // Decompose the 15-point geometry into 24 face-center tetrahedra
        for (int f = 0; f < 6; ++f) {
          const int face_center_id = pamt_face_center_id_table[f];
          for (int tri = 0; tri < 4; ++tri) {
            const int vertex0 = pamt_flux_id_table[f][tri];
            const int vertex1 = pamt_flux_id_table[f][(tri + 1) % 4];

            std::array<IRL::Pt, 4> source_tet_pts = {
                {preimage[pamt_cell_center_id], preimage[vertex0],
                 preimage[vertex1], preimage[face_center_id]}};
            std::array<IRL::Pt, 4> destination_tet_pts = {
                {cell[pamt_cell_center_id], cell[vertex0], cell[vertex1],
                 cell[face_center_id]}};

            MakeSourceTetPositive(source_tet_pts, destination_tet_pts);

            double jacobian = 1.0;
            const Eigen::Matrix3d affine_matrix =
                AffineMatrix(source_tet_pts, destination_tet_pts, jacobian);
            const Eigen::Vector3d y0(destination_tet_pts[0][0],
                                     destination_tet_pts[0][1],
                                     destination_tet_pts[0][2]);
            const Eigen::Vector3d x0(source_tet_pts[0][0], source_tet_pts[0][1],
                                     source_tet_pts[0][2]);
            const Eigen::Vector3d translation = y0 - affine_matrix * x0;
            const auto source_tet =
                IRL::Tet::fromRawPtPointer(4, source_tet_pts.data());

            Real xlo = source_tet_pts[0][0];
            Real ylo = source_tet_pts[0][1];
            Real zlo = source_tet_pts[0][2];
            Real xhi = source_tet_pts[0][0];
            Real yhi = source_tet_pts[0][1];
            Real zhi = source_tet_pts[0][2];
            for (int n = 1; n < 4; ++n) {
              xlo = source_tet_pts[n][0] < xlo ? source_tet_pts[n][0] : xlo;
              ylo = source_tet_pts[n][1] < ylo ? source_tet_pts[n][1] : ylo;
              zlo = source_tet_pts[n][2] < zlo ? source_tet_pts[n][2] : zlo;
              xhi = source_tet_pts[n][0] > xhi ? source_tet_pts[n][0] : xhi;
              yhi = source_tet_pts[n][1] > yhi ? source_tet_pts[n][1] : yhi;
              zhi = source_tet_pts[n][2] > zhi ? source_tet_pts[n][2] : zhi;
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

            // Intersect the preimage tet with old cells and phase geometry.
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

                  const auto cut_volume_moments = IRL::getVolumeMoments<
                      IRL::SeparatedMoments<IRL::VolumeMoments>>(source_tet,
                                                                 local_sep);
                  AddMappedVolumeMoments(cut_volume_moments[0], affine_matrix,
                                         translation, jacobian, M0_l, M1_l);
                  AddMappedVolumeMoments(cut_volume_moments[1], affine_matrix,
                                         translation, jacobian, M0_g, M1_g);

                  if (transport_m2) {
                    const auto cut_general_moments = IRL::getVolumeMoments<
                        IRL::SeparatedMoments<IRL::GeneralMoments3D<2>>>(
                        source_tet, local_sep);
                    AddMappedSecondMoment(cut_volume_moments[0],
                                          cut_general_moments[0], affine_matrix,
                                          translation, jacobian, M2_l);
                    AddMappedSecondMoment(cut_volume_moments[1],
                                          cut_general_moments[1], affine_matrix,
                                          translation, jacobian, M2_g);
                  }
                }
              }
            }
          }
        }

        if (M0_l < 0.0) {
          M0_l = 0.0;
        } else if (M0_l > cell_volume) {
          M0_l = cell_volume;
        }
        moments_array(i, j, k, comp_m0) = M0_l;
        moments_array(i, j, k, comp_vf) = M0_l / cell_volume;

        if (transport_m1 || transport_m2) {
          moments_array(i, j, k, comp_m1_l) = M1_l[0];
          moments_array(i, j, k, comp_m1_l + 1) = M1_l[1];
          moments_array(i, j, k, comp_m1_l + 2) = M1_l[2];
          moments_array(i, j, k, comp_m1_g) = M1_g[0];
          moments_array(i, j, k, comp_m1_g + 1) = M1_g[1];
          moments_array(i, j, k, comp_m1_g + 2) = M1_g[2];
        }

        if (transport_m2) {
          moments_array(i, j, k, comp_m2_l) = M2_l(0, 0);
          moments_array(i, j, k, comp_m2_l + 1) = M2_l(0, 1);
          moments_array(i, j, k, comp_m2_l + 2) = M2_l(0, 2);
          moments_array(i, j, k, comp_m2_l + 3) = M2_l(1, 1);
          moments_array(i, j, k, comp_m2_l + 4) = M2_l(1, 2);
          moments_array(i, j, k, comp_m2_l + 5) = M2_l(2, 2);
          moments_array(i, j, k, comp_m2_g) = M2_g(0, 0);
          moments_array(i, j, k, comp_m2_g + 1) = M2_g(0, 1);
          moments_array(i, j, k, comp_m2_g + 2) = M2_g(0, 2);
          moments_array(i, j, k, comp_m2_g + 3) = M2_g(1, 1);
          moments_array(i, j, k, comp_m2_g + 4) = M2_g(1, 2);
          moments_array(i, j, k, comp_m2_g + 5) = M2_g(2, 2);
        }
      });
    }

    a_moments.FillBoundary(a_geom.periodicity());
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_ADVECTION_PAMT_H_
