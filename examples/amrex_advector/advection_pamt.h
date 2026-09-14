// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_ADVECTION_MOMENTS_H_
#define EXAMPLES_AMREX_ADVECTOR_ADVECTION_MOMENTS_H_

#include <array>
#include <cmath>
#include <sstream>

#include "irl/amrex/sepunion_multifab.h"

#include "examples/amrex_advector/advection_helpers.h"

using namespace amrex;

static constexpr std::array<std::array<int, 4>, 6> tet_id_table = {
    {{{0, 1, 2, 6}},
     {{0, 2, 3, 6}},
     {{0, 3, 7, 6}},
     {{0, 7, 4, 6}},
     {{0, 4, 5, 6}},
     {{0, 5, 1, 6}}}};

// piecewise affine moment transport
struct PAMT {
  static Eigen::Vector3d ToEigen(const IRL::Pt& a_pt) {
    return Eigen::Vector3d({a_pt[0], a_pt[1], a_pt[2]});
  }

  static IRL::Pt ToPt(const Eigen::Vector3d& a_pt) {
    return IRL::Pt(a_pt[0], a_pt[1], a_pt[2]);
  }

  static Eigen::Matrix3d MomentMatrix(
      const IRL::GeneralMoments3D<2>& a_moments) {
    return Eigen::Matrix3d({{a_moments[4], a_moments[5], a_moments[6]},
                            {a_moments[5], a_moments[7], a_moments[8]},
                            {a_moments[6], a_moments[8], a_moments[9]}});
  }

  static Eigen::Matrix3d AffineMatrix(
      const std::array<IRL::Pt, 4>& a_source,
      const std::array<IRL::Pt, 4>& a_destination, double& a_jacobian) {
    const Eigen::Vector3d x0 = ToEigen(a_source[0]);
    const Eigen::Vector3d e1 = ToEigen(a_source[1]) - x0;
    const Eigen::Vector3d e2 = ToEigen(a_source[2]) - x0;
    const Eigen::Vector3d e3 = ToEigen(a_source[3]) - x0;

    const Eigen::Vector3d y0 = ToEigen(a_destination[0]);
    const Eigen::Vector3d f1 = ToEigen(a_destination[1]) - y0;
    const Eigen::Vector3d f2 = ToEigen(a_destination[2]) - y0;
    const Eigen::Vector3d f3 = ToEigen(a_destination[3]) - y0;

    const double source_det = e1.dot(e2.cross(e3));
    const double destination_det = f1.dot(f2.cross(f3));
    a_jacobian = std::abs(destination_det / source_det);

    return (f1 * e2.cross(e3).transpose() + f2 * e3.cross(e1).transpose() +
            f3 * e1.cross(e2).transpose()) /
           source_det;
  }

  static void AddMappedVolumeMoments(
      const IRL::VolumeMoments& a_source_moments,
      const Eigen::Matrix3d& a_affine_matrix,
      const Eigen::Vector3d& a_translation, const double a_jacobian,
      double& a_m0, Eigen::Vector3d& a_m1) {
    const double source_m0 = a_source_moments.volume();
    const Eigen::Vector3d source_m1 = ToEigen(a_source_moments.centroid());
    a_m0 += a_jacobian * source_m0;
    a_m1 +=
        a_jacobian * (a_affine_matrix * source_m1 + a_translation * source_m0);
  }

  static void AddMappedSecondMoment(
      const IRL::VolumeMoments& a_source_volume_moments,
      const IRL::GeneralMoments3D<2>& a_source_general_moments,
      const Eigen::Matrix3d& a_affine_matrix,
      const Eigen::Vector3d& a_translation, const double a_jacobian,
      Eigen::Matrix3d& a_m2) {
    const double source_m0 = a_source_volume_moments.volume();
    const Eigen::Vector3d source_m1 =
        ToEigen(a_source_volume_moments.centroid());
    const Eigen::Matrix3d source_m2 = MomentMatrix(a_source_general_moments);

    a_m2 +=
        a_jacobian *
        (a_affine_matrix * source_m2 * a_affine_matrix.transpose() +
         (a_affine_matrix * source_m1) * a_translation.transpose() +
         a_translation * source_m1.transpose() * a_affine_matrix.transpose() +
         source_m0 * a_translation * a_translation.transpose());
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

      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (band_id_array(i, j, k) == 0.0) {
          return;
        }

        std::array<IRL::Pt, 8> cell;
        std::array<IRL::Pt, 8> preimage;
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

        for (int n = 0; n < ncomp; ++n) {
          moments_array(i, j, k, n) = 0.0;
        }

        double M0_l = 0.0;
        double M0_g = 0.0;
        Eigen::Vector3d M1_l = Eigen::Vector3d::Zero();
        Eigen::Vector3d M1_g = Eigen::Vector3d::Zero();
        Eigen::Matrix3d M2_l = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d M2_g = Eigen::Matrix3d::Zero();

        for (int ktet = 0; ktet < 6; ++ktet) {
          std::array<IRL::Pt, 4> source_tet_pts;
          std::array<IRL::Pt, 4> destination_tet_pts;
          for (int n = 0; n < 4; ++n) {
            const int vertex_id = tet_id_table[ktet][n];
            source_tet_pts[n] = preimage[vertex_id];
            destination_tet_pts[n] = cell[vertex_id];
          }

          double jacobian = 1.0;
          const Eigen::Matrix3d affine_matrix =
              AffineMatrix(source_tet_pts, destination_tet_pts, jacobian);
          const Eigen::Vector3d translation =
              ToEigen(destination_tet_pts[0]) -
              affine_matrix * ToEigen(source_tet_pts[0]);

          const auto source_tet =
              IRL::Tet::fromRawPtPointer(4, source_tet_pts.data());

          Real xlo = source_tet_pts[0][0], ylo = source_tet_pts[0][1],
               zlo = source_tet_pts[0][2];
          Real xhi = source_tet_pts[0][0], yhi = source_tet_pts[0][1],
               zhi = source_tet_pts[0][2];
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

#ifndef NDEBUG
          if (!grown_bx.contains(ilo, jlo, klo) ||
              !grown_bx.contains(ihi, jhi, khi)) {
            std::ostringstream oss;
            oss << "Preimage tetrahedron for cell " << i << " " << j << " " << k
                << " is outside of the grown box\n";
            throw std::runtime_error(oss.str());
          }
#endif

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

#endif  // EXAMPLES_AMREX_ADVECTOR_ADVECTION_MOMENTS_H_
