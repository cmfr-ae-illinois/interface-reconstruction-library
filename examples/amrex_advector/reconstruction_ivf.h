// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_IVF_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_IVF_H_

#include "irl/amrex/sepunion_multifab.h"

#include "examples/amrex_advector/reconstruction_lvira.h"
#include "examples/amrex_advector/reconstruction_plicnet.h"
#include "irl/interface_reconstruction_methods/iterative_jibben.h"

using namespace amrex;

struct iVF {
  // some helper functions
  static IRL::Normal normalFromArray(const Array4<Real const>& normal_array,
                                     const int i, const int j, const int k) {
    IRL::Normal n(normal_array(i, j, k, 0), normal_array(i, j, k, 1),
                  normal_array(i, j, k, 2));
    n.normalize();
    return n;
  }

  static void setNormalArray(const Array4<Real>& normal_array, const int i,
                             const int j, const int k, const IRL::Normal& n) {
    normal_array(i, j, k, 0) = n[0];
    normal_array(i, j, k, 1) = n[1];
    normal_array(i, j, k, 2) = n[2];
  }

  static IRL::Normal localToGlobalNormal(
      const IRL::Normal& a_local_normal,
      const IRL::ReferenceFrame& a_local_frame) {
    const IRL::Normal& e_xi = a_local_frame[0];
    const IRL::Normal& e_eta = a_local_frame[1];
    const IRL::Normal& e_zeta = a_local_frame[2];

    IRL::Normal global_normal(
        a_local_normal[0] * e_xi[0] + a_local_normal[1] * e_eta[0] +
            a_local_normal[2] * e_zeta[0],
        a_local_normal[0] * e_xi[1] + a_local_normal[1] * e_eta[1] +
            a_local_normal[2] * e_zeta[1],
        a_local_normal[0] * e_xi[2] + a_local_normal[1] * e_eta[2] +
            a_local_normal[2] * e_zeta[2]);

    global_normal.normalize();
    return global_normal;
  }

  static IRL::RectangularCuboid makeCell(
      const int i, const int j, const int k,
      const GpuArray<Real, AMREX_SPACEDIM>& problo,
      const GpuArray<Real, AMREX_SPACEDIM>& dx) {
    const double x = problo[0] + i * dx[0];
    const double y = problo[1] + j * dx[1];
    const double z = problo[2] + k * dx[2];

    return IRL::RectangularCuboid::fromBoundingPts(
        IRL::Pt(x, y, z), IRL::Pt(x + dx[0], y + dx[1], z + dx[2]));
  }

  static IRL::PlanarSeparator makeVolumeMatchedPLIC(
      const int i, const int j, const int k, const IRL::Normal& normal,
      const double vfrac, const GpuArray<Real, AMREX_SPACEDIM>& problo,
      const GpuArray<Real, AMREX_SPACEDIM>& dx) {
    auto cell = makeCell(i, j, k, problo, dx);

    IRL::PlanarSeparator plic =
        IRL::PlanarSeparator::fromOnePlane(IRL::Plane(normal, 0.0));

    IRL::setDistanceToMatchVolumeFraction(cell, vfrac, &plic, 1.0e-14);

    return plic;
  }

  static IRL::Polygon polygonFromPLIC(
      const int i, const int j, const int k, const IRL::PlanarSeparator& plic,
      const GpuArray<Real, AMREX_SPACEDIM>& problo,
      const GpuArray<Real, AMREX_SPACEDIM>& dx) {
    auto cell = makeCell(i, j, k, problo, dx);

    return IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(cell, plic,
                                                                plic[0]);
  }

  static void GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr) {
    // initial plic reconstruction
    LVIRA::GetReconstruction(interface, interface_with_ghost, moments, geom,
                             nullptr);

    // scalar field for interface type
    if (scalar_fields) {
      scalar_fields->clear();

      scalar_fields->emplace_back("interface_type", moments.boxArray(),
                                  moments.DistributionMap(), 0);
    }

    // some parameters
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const Real vol = dx[0] * dx[1] * dx[2];
    const Real vol_inv = 1.0 / vol;

    constexpr int nlayers = 1;
    constexpr int max_iter = 5;

    constexpr double alpha_deg = 50.0;
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double alpha_rad = alpha_deg * pi / 180.0;
    const double cos_alpha = std::cos(alpha_rad);

    constexpr double eps_kappa = 1.0e-3;

    constexpr int nstencil =
        (1 + 2 * nlayers) * (1 + 2 * nlayers) * (1 + 2 * nlayers);

    // some fields for storing data
    const BoxArray& ba = interface.boxArray();
    const DistributionMapping& dm = interface.DistributionMap();

    MultiFab init_normal(ba, dm, 3, interface_with_ghost.nGrowVect());
    MultiFab current_normal(ba, dm, 3, interface_with_ghost.nGrowVect());

    iMultiFab marked(ba, dm, 1, 0);
    iMultiFab fallback_to_plic(ba, dm, 1, 0);

    MultiFab kappa_r(ba, dm, 1, 0);
    MultiFab kappa_rm1(ba, dm, 1, 0);
    MultiFab kappa_rm2(ba, dm, 1, 0);

    init_normal.setVal(0.0);
    current_normal.setVal(0.0);

    marked.setVal(0);
    fallback_to_plic.setVal(0);

    kappa_r.setVal(0.0);
    kappa_rm1.setVal(0.0);
    kappa_rm2.setVal(0.0);

    // initializing fields using plic reconstruction
    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);
      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      Array4<Real> init_normal_array = init_normal.array(mfi);
      Array4<Real> current_normal_array = current_normal.array(mfi);

      Array4<int> marked_array = marked.array(mfi);
      Array4<int> fallback_array = fallback_to_plic.array(mfi);

      Array4<Real> kappa_r_array = kappa_r.array(mfi);
      Array4<Real> kappa_rm1_array = kappa_rm1.array(mfi);
      Array4<Real> kappa_rm2_array = kappa_rm2.array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        marked_array(i, j, k) = 0;
        fallback_array(i, j, k) = 0;

        kappa_r_array(i, j, k) = 0.0;
        kappa_rm1_array(i, j, k) = 0.0;
        kappa_rm2_array(i, j, k) = 0.0;

        const double vfrac = moments_array(i, j, k) * vol_inv;

        if (vfrac < IRL::global_constants::VF_LOW) {
          interface_array(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
          return;
        } else if (vfrac > IRL::global_constants::VF_HIGH) {
          interface_array(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          return;
        }

        const IRL::Normal n =
            interface_with_ghost_array(i, j, k).getPlane().normal();

        setNormalArray(init_normal_array, i, j, k, n);
        setNormalArray(current_normal_array, i, j, k, n);

        marked_array(i, j, k) = 1;
      });
    }

    init_normal.FillBoundary(geom.periodicity());
    current_normal.FillBoundary(geom.periodicity());

    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);

    // main ippic iteration loop
    for (int iter = 1; iter <= max_iter; ++iter) {
      // if iter > 1 and cell is marked then use plic with updated normal
      if (iter > 1) {
        for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);
          Array4<Real const> moments_array = moments.const_array(mfi);
          Array4<Real const> current_normal_array =
              current_normal.const_array(mfi);
          Array4<int const> marked_array = marked.const_array(mfi);

          const Box& bx = mfi.tilebox();

          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (marked_array(i, j, k) == 0) {
              return;
            }

            const double vfrac = moments_array(i, j, k) * vol_inv;

            const IRL::Normal n =
                normalFromArray(current_normal_array, i, j, k);

            IRL::PlanarSeparator plic =
                makeVolumeMatchedPLIC(i, j, k, n, vfrac, problo, dx);

            interface_array(i, j, k) = plic;
          });
        }

        interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                       interface.nGrowVect());
        interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
      }

      // build filtered stencil and update normal for each marked cell
      for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);
        Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
            interface_with_ghost.const_array(mfi);

        Array4<Real const> moments_array = moments.const_array(mfi);

        Array4<Real const> init_normal_array = init_normal.const_array(mfi);
        Array4<Real> current_normal_array = current_normal.array(mfi);

        Array4<int> marked_array = marked.array(mfi);
        Array4<int> fallback_array = fallback_to_plic.array(mfi);

        Array4<Real> kappa_r_array = kappa_r.array(mfi);
        Array4<Real> kappa_rm1_array = kappa_rm1.array(mfi);
        Array4<Real> kappa_rm2_array = kappa_rm2.array(mfi);

        const Box& bx = mfi.tilebox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          if (marked_array(i, j, k) == 0) {
            return;
          }

          const double vfrac = moments_array(i, j, k) * vol_inv;

          const IRL::Normal init_n =
              normalFromArray(init_normal_array, i, j, k);
          const IRL::Normal current_n =
              normalFromArray(current_normal_array, i, j, k);

          // under resolved condition based on normal deviation
          if ((init_n * current_n) < cos_alpha) {
            fallback_array(i, j, k) = 1;
            marked_array(i, j, k) = 0;

            IRL::PlanarSeparator initial_plic =
                makeVolumeMatchedPLIC(i, j, k, init_n, vfrac, problo, dx);

            interface_array(i, j, k) = initial_plic;
            return;
          }

          // building neighborhood
          IRL::JibbenNeighborhood neighborhood;
          neighborhood.reserve(nstencil);

          int count = 0;
          int center_count = -1;

          IRL::Pt local_datum;
          IRL::ReferenceFrame local_frame;

          for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
            for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
              for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
                const double vfrac_local = moments_array(ii, jj, kk) * vol_inv;

                if (vfrac_local < IRL::global_constants::VF_LOW ||
                    vfrac_local > IRL::global_constants::VF_HIGH) {
                  continue;
                }

                const IRL::Normal neighbor_init_n =
                    normalFromArray(init_normal_array, ii, jj, kk);

                const double alignment = neighbor_init_n * current_n;

                if (alignment < cos_alpha) {
                  continue;
                }

                const auto neighbor_plic = IRL::PlanarSeparator::fromOnePlane(
                    interface_with_ghost_array(ii, jj, kk).getPlane());

                const IRL::Polygon neighbor_polygon =
                    polygonFromPLIC(ii, jj, kk, neighbor_plic, problo, dx);

                if (neighbor_polygon.getNumberOfVertices() == 0) {
                  continue;
                }

                neighborhood.addMember(neighbor_polygon);

                if (ii == i && jj == j && kk == k) {
                  center_count = count;
                  neighborhood.setCenterOfStencil(count);

                  local_datum = neighbor_polygon.calculateCentroid();
                  local_frame = IRL::ReferenceFrame::fromNormal(current_n);
                }

                ++count;
              }
            }
          }

          // underresolved condition based on neighborhood size
          if (count < 6 || center_count < 0) {
            fallback_array(i, j, k) = 1;
            marked_array(i, j, k) = 0;

            IRL::PlanarSeparator initial_plic =
                makeVolumeMatchedPLIC(i, j, k, init_n, vfrac, problo, dx);

            interface_array(i, j, k) = initial_plic;
            return;
          }

          // paraboloid fit
          neighborhood.localize();

          IRL::iJibben_3D ijibben(&neighborhood);
          ijibben.getParaboloidCoefficients2();

          const std::pair<double, IRL::Normal> fit_result =
              ijibben.computeAveragedCurvatureAndNormal();

          const double kappa_new = fit_result.first;
          double kappa_accepted = kappa_new;

          if (iter > 2) {
            const double new_change =
                std::abs(kappa_new - kappa_rm1_array(i, j, k));
            const double old_change =
                std::abs(kappa_rm1_array(i, j, k) - kappa_rm2_array(i, j, k));

            if (new_change < old_change) {
              kappa_accepted = kappa_new;
            } else {
              kappa_accepted = kappa_rm1_array(i, j, k);
            }
          }
          kappa_r_array(i, j, k) = kappa_accepted;

          // Update global normal using averaged local normal
          const IRL::Normal updated_global_normal =
              localToGlobalNormal(fit_result.second, local_frame);

          setNormalArray(current_normal_array, i, j, k, updated_global_normal);

          // Convergence condition
          if (iter > 2) {
            const double denom =
                std::max(std::abs(kappa_r_array(i, j, k)), 1.0e-14);

            const double rel_change =
                std::abs(kappa_r_array(i, j, k) - kappa_rm1_array(i, j, k)) /
                denom;

            if (rel_change < eps_kappa) {
              marked_array(i, j, k) = 0;
            }
          }

          // Shift curvature history
          kappa_rm2_array(i, j, k) = kappa_rm1_array(i, j, k);
          kappa_rm1_array(i, j, k) = kappa_r_array(i, j, k);
        });
      }
      current_normal.FillBoundary(geom.periodicity());

      interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                     interface.nGrowVect());
      interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
    }

    // final plic reconstruction with updated normals (if any)
    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);
      Array4<Real const> current_normal_array = current_normal.const_array(mfi);

      Array4<int const> fallback_array = fallback_to_plic.const_array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double vfrac = moments_array(i, j, k) * vol_inv;

        if (vfrac < IRL::global_constants::VF_LOW ||
            vfrac > IRL::global_constants::VF_HIGH) {
          return;
        }

        if (fallback_array(i, j, k) != 0) {
          return;
        }

        const IRL::Normal current_n =
            normalFromArray(current_normal_array, i, j, k);

        IRL::PlanarSeparator final_plic =
            makeVolumeMatchedPLIC(i, j, k, current_n, vfrac, problo, dx);

        interface_array(i, j, k) = final_plic;
      });
    }

    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);

    // final ppic reconstruction
    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);
      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      Array4<Real const> init_normal_array = init_normal.const_array(mfi);
      Array4<Real const> current_normal_array = current_normal.const_array(mfi);

      Array4<int const> fallback_array = fallback_to_plic.const_array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double vfrac = moments_array(i, j, k) * vol_inv;

        if (vfrac < IRL::global_constants::VF_LOW ||
            vfrac > IRL::global_constants::VF_HIGH) {
          return;
        }

        if (fallback_array(i, j, k) != 0) {
          return;
        }

        const IRL::Normal current_n =
            normalFromArray(current_normal_array, i, j, k);

        IRL::JibbenNeighborhood neighborhood;
        neighborhood.reserve(nstencil);

        int count = 0;
        int center_count = -1;

        IRL::Pt local_datum;
        IRL::ReferenceFrame local_frame;

        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              const double vfrac_local = moments_array(ii, jj, kk) * vol_inv;

              if (vfrac_local < IRL::global_constants::VF_LOW ||
                  vfrac_local > IRL::global_constants::VF_HIGH) {
                continue;
              }

              const IRL::Normal neighbor_init_n =
                  normalFromArray(init_normal_array, ii, jj, kk);

              const double alignment = neighbor_init_n * current_n;

              if (alignment < cos_alpha) {
                continue;
              }

              const auto neighbor_plic = IRL::PlanarSeparator::fromOnePlane(
                  interface_with_ghost_array(ii, jj, kk).getPlane());

              const IRL::Polygon neighbor_polygon =
                  polygonFromPLIC(ii, jj, kk, neighbor_plic, problo, dx);

              if (neighbor_polygon.getNumberOfVertices() == 0) {
                continue;
              }

              neighborhood.addMember(neighbor_polygon);

              if (ii == i && jj == j && kk == k) {
                center_count = count;
                neighborhood.setCenterOfStencil(count);

                local_datum = neighbor_polygon.calculateCentroid();
                local_frame = IRL::ReferenceFrame::fromNormal(current_n);
              }

              ++count;
            }
          }
        }

        if (count < 6 || center_count < 0) {
          const IRL::Normal init_n =
              normalFromArray(init_normal_array, i, j, k);

          IRL::PlanarSeparator initial_plic =
              makeVolumeMatchedPLIC(i, j, k, init_n, vfrac, problo, dx);

          interface_array(i, j, k) = initial_plic;
          return;
        }

        neighborhood.localize();

        IRL::iJibben_3D ijibben(&neighborhood);
        ijibben.getParaboloidCoefficients();

        interface_array(i, j, k) = ijibben.makeParaboloid2(&neighborhood);

        const auto cell = makeCell(i, j, k, problo, dx);

        IRL::setDistanceToMatchVolumeFraction(
            cell, vfrac, &interface_array(i, j, k), 1.0e-14);
      });
    }
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_IVF_H_
