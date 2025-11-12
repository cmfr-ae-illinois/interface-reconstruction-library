// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <iostream>

#include "examples/2d_advector/deformation_2d.h"
#include "examples/2d_advector/irl2d.h"
#include "examples/2d_advector/oscillation_2d.h"
#include "examples/2d_advector/reconstruction_types.h"
#include "examples/2d_advector/rotation_2d.h"
#include "examples/2d_advector/vof_advection.h"
#include "examples/2d_advector/vtk.h"

void advectVOF(const std::string& a_simulation_type,
               const std::string& a_advection_method,
               const std::string& a_reconstruction_method, const double a_dt,
               const double a_time, const Data<double>& a_U,
               const Data<double>& a_V, Data<IRL2D::Moments>* a_liquid_moments,
               Data<IRL2D::Moments>* a_gas_moments,
               Data<IRL2D::Parabola>* a_interface) {
  if (a_advection_method == "SemiLagL" || a_advection_method == "SemiLagQ") {
    SemiLag::advectVOF(a_simulation_type, a_advection_method,
                       a_reconstruction_method, a_dt, a_time, a_U, a_V,
                       a_liquid_moments, a_gas_moments, a_interface);
  } else if (a_advection_method == "FullLagL" ||
             a_advection_method == "FullLagQ") {
    FullLag::advectVOF(a_simulation_type, a_advection_method,
                       a_reconstruction_method, a_dt, a_time, a_U, a_V,
                       a_liquid_moments, a_gas_moments, a_interface);
  } else if (a_advection_method == "ReSyFullLagL"){
    ReSyFullLagL::advectVOF(a_simulation_type, a_advection_method,
                       a_reconstruction_method, a_dt, a_time, a_U, a_V,
                       a_liquid_moments, a_gas_moments, a_interface);
  } else {
    std::cout << "Unknown advection method of : " << a_reconstruction_method
              << '\n';
    std::cout
        << "Valid entries are: SemiLagL, SemiLagQ, FullLagL, FullLagQ, ReSyFullLagL. \n";
    std::exit(-1);
  }
}

void SemiLag::advectVOF(const std::string& a_simulation_type,
                        const std::string& a_advection_method,
                        const std::string& a_reconstruction_method,
                        const double a_dt, const double a_time,
                        const Data<double>& a_U, const Data<double>& a_V,
                        Data<IRL2D::Moments>* a_liquid_moments,
                        Data<IRL2D::Moments>* a_gas_moments,
                        Data<IRL2D::Parabola>* a_interface) {
  const bool correct_fluxes = true;
  const BasicMesh& mesh = a_liquid_moments->getMesh();

  IRL2D::BezierList (*CreateFluxCell)(
      const IRL2D::Vec& P00, const IRL2D::Vec& P10, const double dt,
      const double time,
      const IRL2D::Vec (*vel)(const double t, const IRL2D::Vec& P),
      const IRL2D::Mat (*grad_vel)(const double t, const IRL2D::Vec& P),
      const bool correct_area, const double exact_area);

  if (a_advection_method == "SemiLagL") {
    CreateFluxCell = IRL2D::CreateLinearFluxCell;
  } else if (a_advection_method == "SemiLagQ") {
    CreateFluxCell = IRL2D::CreateFluxCell;
  }

  // for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
  //     const auto cell =
  //         IRL2D::RectangleFromBounds(IRL2D::Vec(mesh.x(i), mesh.y(j)),
  //                                    IRL2D::Vec(mesh.x(i + 1), mesh.y(j +
  //                                    1)));
  //     const auto cell_moments = IRL2D::ComputeMoments(cell);
  //     auto moments = IRL2D::ComputeMoments(cell, (*a_interface)(i, j));
  //     moments.m0() = (*a_liquid_moments)(i, j).m0();
  //     (*a_liquid_moments)(i, j) = moments;
  //     (*a_gas_moments)(i, j) = cell_moments - (*a_liquid_moments)(i, j);
  //   }
  // }
  // a_liquid_moments->updateBorder();
  // a_gas_moments->updateBorder();

  const IRL2D::Vec (*getExactVelocity2D)(const double t, const IRL2D::Vec& P);
  const IRL2D::Mat (*getExactGradient2D)(const double t, const IRL2D::Vec& P);

  if (a_simulation_type == "Rotation2D") {
    getExactVelocity2D = Rotation2D::getExactVelocity2D;
    getExactGradient2D = Rotation2D::getExactVelocityGradient2D;
  } else if (a_simulation_type == "Oscillation2D") {
    getExactVelocity2D = Oscillation2D::getExactVelocity2D;
    getExactGradient2D = Oscillation2D::getExactVelocityGradient2D;
  } else if (a_simulation_type == "Deformation2D") {
    getExactVelocity2D = Deformation2D::getExactVelocity2D;
    getExactGradient2D = Deformation2D::getExactVelocityGradient2D;
  }

  // Allocate storage for face fluxes
  Data<IRL2D::Moments> liq_flux[2] = {Data<IRL2D::Moments>(&mesh),
                                      Data<IRL2D::Moments>(&mesh)};
  Data<IRL2D::Moments> gas_flux[2] = {Data<IRL2D::Moments>(&mesh),
                                      Data<IRL2D::Moments>(&mesh)};

  // Calculate CFL
  double CFL = -DBL_MAX;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      CFL = std::fmax(CFL, std::fmax(a_U(i, j) * a_dt / mesh.dx(),
                                     a_V(i, j) * a_dt / mesh.dy()));
    }
  }

  // Generate band that is ceil(CFL) thick around interface
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (*a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }
  band.updateBorder();
  const int nlayers = 1 + static_cast<int>(std::ceil(CFL));
  for (int n = 0; n < nlayers; ++n) {
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        if (band(i, j) == 0) {
          for (int ii = -1; ii <= 1; ++ii) {
            for (int jj = -1; jj <= 1; ++jj) {
              if (band(i + ii, j + jj) == n + 1) {
                band(i, j) = n + 2;
              }
            }
          }
        }
      }
    }
    band.updateBorder();
  }

#ifdef USE_MPI
  int nmixed_global = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
        nmixed_global++;
      }
    }
  }

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nmixed_local = std::max(nmixed_global / size, 1);
  std::vector<int> proc_offset(size + 1);
  proc_offset[0] = 0;
  for (int r = 0; r < size; r++)
    proc_offset[r + 1] = proc_offset[r] + nmixed_local;
  proc_offset[size] = nmixed_global;
  for (int r = 1; r < size + 1; r++)
    proc_offset[r] = std::min(proc_offset[r], nmixed_global);
  nmixed_local = proc_offset[rank + 1] - proc_offset[rank];

  const int size_moments = 2 * 2 * 6;
  std::vector<double> flux_local(size_moments * nmixed_local);
  std::vector<double> flux_global(size_moments * nmixed_global);
#endif

  // Compute fluxes
  std::vector<IRL2D::BezierList> fluxes, clipped_fluxes;
  Data<IRL2D::BezierList> pre_images;

  int count = 0, count_local = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      // Initialize fluxes to zero
      for (int dim = 0; dim < 2; ++dim) {
        (liq_flux[dim])(i, j) = IRL2D::Moments();
        (gas_flux[dim])(i, j) = IRL2D::Moments();
      }
      // Only solve advection in narrow band near the interface
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
#ifdef USE_MPI
        if (count >= proc_offset[rank] && count < proc_offset[rank + 1]) {
#endif
          const auto x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
          const auto x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
          const auto cell = IRL2D::RectangleFromBounds(x0, x1);

          // Initialize face fluxes
          IRL2D::BezierList face_cell[2];
          face_cell[0] = CreateFluxCell(
              cell[3].first, cell[0].first, -a_dt, a_time + a_dt,
              getExactVelocity2D, getExactGradient2D, correct_fluxes,
              IRL2D::IntegrateFlux(cell[3].first, cell[0].first, a_dt, a_time,
                                   getExactVelocity2D));
          face_cell[1] = CreateFluxCell(
              cell[0].first, cell[1].first, -a_dt, a_time + a_dt,
              getExactVelocity2D, getExactGradient2D, correct_fluxes,
              IRL2D::IntegrateFlux(cell[0].first, cell[1].first, a_dt, a_time,
                                   getExactVelocity2D));
          // if ((i == 49 || i == 50) && (j == 31 || j == 32)) {
          // fluxes.push_back(face_cell[0]);
          // fluxes.push_back(face_cell[1]);
          // }

          for (int dim = 0; dim < 2; ++dim) {
            // Compute flux bounding boxes
            const auto bbox = IRL2D::BoundingBox(face_cell[dim]);
            int idmin[2], idmax[2];
            mesh.getIndices(bbox.first, idmin);
            mesh.getIndices(bbox.second, idmax);

            // Compute liquid fluxes by intersection on a n x n neighborhood
            // for (int ii = std::max(idmin[0], mesh.imino());
            //      ii <= std::min(idmax[0], mesh.imaxo()); ++ii) {
            //   for (int jj = std::max(idmin[1], mesh.imino());
            //        jj <= std::min(idmax[1], mesh.imaxo()); ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ii++) {
              for (int jj = j - nlayers; jj <= j + nlayers; jj++) {
                const auto xn0 = IRL2D::Vec(mesh.x(ii), mesh.y(jj));
                const auto xn1 = IRL2D::Vec(mesh.x(ii + 1), mesh.y(jj + 1));
                (liq_flux[dim])(i, j) += IRL2D::ComputeMoments(
                    face_cell[dim], xn0, xn1, (*a_interface)(ii, jj));
                // if ((i == 49 || i == 50) && (j == 31 || j == 32)) {
                // clipped_fluxes.push_back(ClipByRectangleAndParabola(
                //     face_cell[dim], xn0, xn1, (*a_interface)(ii, jj)));
                // }
              }
            }

#ifndef USE_MPI
            // Update gas fluxes
            (gas_flux[dim])(i, j) +=
                IRL2D::ComputeMoments(face_cell[dim]) - (liq_flux[dim])(i, j);
#else
          const auto fm = IRL2D::ComputeMoments(face_cell[dim]);
          flux_local[count_local++] = liq_flux[dim](i, j).m0();
          flux_local[count_local++] = liq_flux[dim](i, j).m1()[0];
          flux_local[count_local++] = liq_flux[dim](i, j).m1()[1];
          flux_local[count_local++] = liq_flux[dim](i, j).m2()[0][0];
          flux_local[count_local++] = liq_flux[dim](i, j).m2()[1][0];
          flux_local[count_local++] = liq_flux[dim](i, j).m2()[1][1];
          flux_local[count_local++] = fm.m0() - liq_flux[dim](i, j).m0();
          flux_local[count_local++] = fm.m1()[0] - liq_flux[dim](i, j).m1()[0];
          flux_local[count_local++] = fm.m1()[1] - liq_flux[dim](i, j).m1()[1];
          flux_local[count_local++] =
              fm.m2()[0][0] - liq_flux[dim](i, j).m2()[0][0];
          flux_local[count_local++] =
              fm.m2()[1][0] - liq_flux[dim](i, j).m2()[1][0];
          flux_local[count_local++] =
              fm.m2()[1][1] - liq_flux[dim](i, j).m2()[1][1];
#endif
          }
#ifdef USE_MPI
        }
        count++;
#endif
      }
    }
  }
  if (fluxes.size() > 0) {
    IRL2D::ToVTK(fluxes, "fluxes");
    IRL2D::ToVTK(clipped_fluxes, "clipped_fluxes");
  }

#ifdef USE_MPI
  std::vector<int> proc_count(size);
  for (int r = 0; r < size; r++) {
    proc_count[r] = size_moments * (proc_offset[r + 1] - proc_offset[r]);
    proc_offset[r] = size_moments * proc_offset[r];
  }
  MPI_Allgatherv(flux_local.data(), flux_local.size(), MPI_DOUBLE,
                 flux_global.data(), proc_count.data(), proc_offset.data(),
                 MPI_DOUBLE, MPI_COMM_WORLD);

  count_local = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
        for (int dim = 0; dim < 2; ++dim) {
          (liq_flux[dim])(i, j).m0() = flux_global[count_local++];
          (liq_flux[dim])(i, j).m1()[0] = flux_global[count_local++];
          (liq_flux[dim])(i, j).m1()[1] = flux_global[count_local++];
          (liq_flux[dim])(i, j).m2()[0][0] = flux_global[count_local++];
          (liq_flux[dim])(i, j).m2()[1][0] = flux_global[count_local++];
          (liq_flux[dim])(i, j).m2()[1][1] = flux_global[count_local++];
          (liq_flux[dim])(i, j).m2()[0][1] = (liq_flux[dim])(i, j).m2()[1][0];
          (gas_flux[dim])(i, j).m0() = flux_global[count_local++];
          (gas_flux[dim])(i, j).m1()[0] = flux_global[count_local++];
          (gas_flux[dim])(i, j).m1()[1] = flux_global[count_local++];
          (gas_flux[dim])(i, j).m2()[0][0] = flux_global[count_local++];
          (gas_flux[dim])(i, j).m2()[1][0] = flux_global[count_local++];
          (gas_flux[dim])(i, j).m2()[1][1] = flux_global[count_local++];
          (gas_flux[dim])(i, j).m2()[0][1] = (gas_flux[dim])(i, j).m2()[1][0];
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // liq_flux[0].updateBorder();
  // liq_flux[1].updateBorder();
  // gas_flux[0].updateBorder();
  // gas_flux[1].updateBorder();

  // Now calculate VOF from the face fluxes.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double cell_volume = mesh.cell_volume();
      const auto cell_centroid = IRL2D::Vec(mesh.xm(i), mesh.ym(j));
      const double liquid_volume_fraction =
          (*a_liquid_moments)(i, j).m0() / cell_volume;
      if (band(i, j) > 0) {
        // Compute moments in transported cell image
        (*a_liquid_moments)(i, j) =
            (*a_liquid_moments)(i, j) + (liq_flux[0])(i, j) -
            (liq_flux[0])(i + 1, j) + (liq_flux[1])(i, j) -
            (liq_flux[1])(i, j + 1);
        (*a_gas_moments)(i, j) = (*a_gas_moments)(i, j) + (gas_flux[0])(i, j) -
                                 (gas_flux[0])(i + 1, j) + (gas_flux[1])(i, j) -
                                 (gas_flux[1])(i, j + 1);

        const std::array<Data<IRL2D::Moments>*, 2> moment_list(
            {a_liquid_moments, a_gas_moments});

        for (int m = 0; m < 2; m++) {
          const auto moments = moment_list[m];
          const auto M0 = (*moments)(i, j).m0();
          const auto M1 = (*moments)(i, j).m1();
          const auto M2 = (*moments)(i, j).m2();

          // Correct moment 0
          double M0_final = M0;  // std::clamp(M0, 0.0, cell_volume);
          if (M0 < -1.0e-12 || M0 > 1.0 + 1.0e-12) {
            // std::cout << std::setprecision(15) << "M0(" << i << "," << j
            //           << ") = " << M0 << std::endl;
            // std::cout << "Parabola: " << (*a_interface)(i, j) << std::endl;
            // exit(0);
          }

          auto x0_k1 = M1 / IRL::safelyEpsilon(M0);
          const auto u0_k1 = getExactVelocity2D(a_time, x0_k1);
          const auto gradu0_k1 = getExactGradient2D(a_time, x0_k1);
          const auto x0_k2 = x0_k1 + a_dt * u0_k1 * 0.5;
          const auto u0_k2 = getExactVelocity2D(a_time + 0.5 * a_dt, x0_k2);
          const auto gradu0_k2 = getExactGradient2D(a_time + 0.5 * a_dt, x0_k2);
          const auto x0_k3 = x0_k1 + a_dt * u0_k2 * 0.5;
          const auto u0_k3 = getExactVelocity2D(a_time + 0.5 * a_dt, x0_k3);
          const auto gradu0_k3 = getExactGradient2D(a_time + 0.5 * a_dt, x0_k3);
          const auto x0_k4 = x0_k1 + a_dt * u0_k3;
          const auto u0_k4 = getExactVelocity2D(a_time + a_dt, x0_k4);
          const auto gradu0_k4 = getExactGradient2D(a_time + a_dt, x0_k4);

          const auto M1_final =
              M1 + a_dt * M0_final *
                       (u0_k1 + 2.0 * u0_k2 + 2.0 * u0_k3 + u0_k4) / 6.0;

          const auto M2t0_k1 = M0_final * IRL2D::outer_product(x0_k1, u0_k1);
          const auto M2t0_k2 = M0_final * IRL2D::outer_product(x0_k2, u0_k2);
          const auto M2t0_k3 = M0_final * IRL2D::outer_product(x0_k3, u0_k3);
          const auto M2t0_k4 = M0_final * IRL2D::outer_product(x0_k4, u0_k4);
          const auto M2t1_k1 = -M0_final * IRL2D::outer_product(x0_k1, x0_k1) *
                               gradu0_k1.transpose();
          const auto M2t1_k2 = -M0_final * IRL2D::outer_product(x0_k2, x0_k2) *
                               gradu0_k2.transpose();
          const auto M2t1_k3 = -M0_final * IRL2D::outer_product(x0_k3, x0_k3) *
                               gradu0_k3.transpose();
          const auto M2t1_k4 = -M0_final * IRL2D::outer_product(x0_k4, x0_k4) *
                               gradu0_k4.transpose();
          const auto M2t2_k1 = M2 * gradu0_k1.transpose();
          const auto M2t2_k2 =
              (M2 + a_dt * 0.5 *
                        (M2t0_k1 + M2t0_k1.transpose() + M2t1_k1 +
                         M2t1_k1.transpose() + M2t2_k1 + M2t2_k1.transpose())) *
              gradu0_k2.transpose();
          const auto M2t2_k3 =
              (M2 + a_dt * 0.5 *
                        (M2t0_k2 + M2t0_k2.transpose() + M2t1_k2 +
                         M2t1_k2.transpose() + M2t2_k2 + M2t2_k2.transpose())) *
              gradu0_k3.transpose();
          const auto M2t2_k4 = (M2 + a_dt * (M2t0_k3 + M2t0_k3.transpose() +
                                             M2t1_k3 + M2t1_k3.transpose() +
                                             M2t2_k3 + M2t2_k3.transpose())) *
                               gradu0_k4.transpose();
          const auto M2_rk4 =
              a_dt *
              ((M2t0_k1 + 2.0 * M2t0_k2 + 2.0 * M2t0_k3 + M2t0_k4) +
               (M2t1_k1 + 2.0 * M2t1_k2 + 2.0 * M2t1_k3 + M2t1_k4) +
               (M2t2_k1 + 2.0 * M2t2_k2 + 2.0 * M2t2_k3 + M2t2_k4)) *
              (1.0 / 6.0);
          const auto M2_final = M2 + M2_rk4 + M2_rk4.transpose();

          (*moments)(i, j).m0() = M0_final;
          (*moments)(i, j).m1() = M1_final;
          (*moments)(i, j).m2() = M2_final;
        }
      }
    }
  }

  // const double cell_volume = mesh.cell_volume();
  // std::vector<std::tuple<double, int, int>> sorted_volumes;
  // sorted_volumes.resize(mesh.getNx() * mesh.getNy());
  // count = 0;
  // double deficit = 0.0, extra = 0.0;
  // for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
  //   for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
  //     sorted_volumes[count++] =
  //         std::make_tuple((*a_liquid_moments)(i, j).m0(), i, j);
  //     if ((*a_liquid_moments)(i, j).m0() < 0.0) {
  //       deficit -= (*a_liquid_moments)(i, j).m0();
  //     } else if ((*a_liquid_moments)(i, j).m0() > cell_volume) {
  //       extra += (*a_liquid_moments)(i, j).m0() - cell_volume;
  //     }
  //   }
  // }
  // std::sort(sorted_volumes.begin(), sorted_volumes.end());

  // int nzero = -1, none = -1;
  // for (int n = 0; n < mesh.getNx() * mesh.getNy() - 1; n++) {
  //   const auto [volume0, i0, j0] = sorted_volumes[n];
  //   const auto [volume1, i1, j1] = sorted_volumes[n + 1];
  //   if (volume0 <= 0.0 && volume1 > 0.0) {
  //     nzero = n;
  //   } else if (volume0 <= cell_volume && volume1 > cell_volume) {
  //     none = n;
  //   }
  // }

  // if (deficit > 0.0) {
  //   for (int n = 0; n <= nzero; n++) {
  //     const auto [volume, i, j] = sorted_volumes[n];
  //     (*a_liquid_moments)(i, j).m0() = 0.0;
  //   }
  //   for (int n = nzero + 1; n <= none; n++) {
  //     const auto [volume, i, j] = sorted_volumes[n];
  //     if (volume < deficit) {
  //       (*a_liquid_moments)(i, j).m0() = 0.0;
  //       deficit -= volume;
  //     } else {
  //       (*a_liquid_moments)(i, j).m0() -= deficit;
  //       deficit = 0.0;
  //       break;
  //     }
  //   }
  // }

  // if (extra > 0.0) {
  //   for (int n = none + 1; n < mesh.getNx() * mesh.getNy(); n++) {
  //     const auto [volume, i, j] = sorted_volumes[n];
  //     (*a_liquid_moments)(i, j).m0() = cell_volume;
  //   }
  //   for (int n = none; n > nzero; n--) {
  //     const auto [volume, i, j] = sorted_volumes[n];
  //     if (cell_volume < volume + extra) {
  //       (*a_liquid_moments)(i, j).m0() = cell_volume;
  //       extra -= cell_volume - volume;
  //     } else {
  //       (*a_liquid_moments)(i, j).m0() += extra;
  //       extra = 0.0;
  //       break;
  //     }
  //   }
  // }

  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();
  correctMomentLocations(a_liquid_moments);
  correctMomentLocations(a_gas_moments);
}

void FullLag::advectVOF(const std::string& a_simulation_type,
                        const std::string& a_advection_method,
                        const std::string& a_reconstruction_method,
                        const double a_dt, const double a_time,
                        const Data<double>& a_U, const Data<double>& a_V,
                        Data<IRL2D::Moments>* a_liquid_moments,
                        Data<IRL2D::Moments>* a_gas_moments,
                        Data<IRL2D::Parabola>* a_interface) {
  const double tol = 100. * std::numeric_limits<double>::epsilon();
  const bool correct_fluxes = true;
  const BasicMesh& mesh = a_liquid_moments->getMesh();

  IRL2D::BezierList (*CreatePreImage)(
      const IRL2D::Vec& x0, const IRL2D::Vec& x1, const double dt,
      const double time,
      const IRL2D::Vec (*vel)(const double t, const IRL2D::Vec& P),
      const IRL2D::Mat (*grad_vel)(const double t, const IRL2D::Vec& P),
      const bool correct_area, const std::array<double, 4>& exact_area);

  if (a_advection_method == "FullLagL") {
    CreatePreImage = IRL2D::CreateLinearPreImage;
  } else if (a_advection_method == "FullLagQ") {
    CreatePreImage = IRL2D::CreatePreImage;
  }

  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      const auto cell =
          IRL2D::RectangleFromBounds(IRL2D::Vec(mesh.x(i), mesh.y(j)),
                                     IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
      const auto cell_moments = IRL2D::ComputeMoments(cell);
      auto moments = IRL2D::ComputeMoments(cell, (*a_interface)(i, j));
      moments.m0() = (*a_liquid_moments)(i, j).m0();
      (*a_liquid_moments)(i, j) = moments;
      (*a_gas_moments)(i, j) = cell_moments - (*a_liquid_moments)(i, j);
    }
  }
  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();

  const IRL2D::Vec (*getExactVelocity2D)(const double t, const IRL2D::Vec& P);
  const IRL2D::Mat (*getExactGradient2D)(const double t, const IRL2D::Vec& P);

  if (a_simulation_type == "Rotation2D") {
    getExactVelocity2D = Rotation2D::getExactVelocity2D;
    getExactGradient2D = Rotation2D::getExactVelocityGradient2D;
  } else if (a_simulation_type == "Oscillation2D") {
    getExactVelocity2D = Oscillation2D::getExactVelocity2D;
    getExactGradient2D = Oscillation2D::getExactVelocityGradient2D;
  } else if (a_simulation_type == "Deformation2D") {
    getExactVelocity2D = Deformation2D::getExactVelocity2D;
    getExactGradient2D = Deformation2D::getExactVelocityGradient2D;
  }

  // Calculate CFL
  double CFL = -DBL_MAX;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      CFL = std::fmax(CFL, std::fmax(a_U(i, j) * a_dt / mesh.dx(),
                                     a_V(i, j) * a_dt / mesh.dy()));
    }
  }

  // Generate band that is ceil(CFL) thick around interface
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (*a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }
  band.updateBorder();
  const int nlayers = 1 + static_cast<int>(std::ceil(CFL));
  for (int n = 0; n < nlayers; ++n) {
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        if (band(i, j) == 0) {
          for (int ii = -1; ii <= 1; ++ii) {
            for (int jj = -1; jj <= 1; ++jj) {
              if (band(i + ii, j + jj) == n + 1) {
                band(i, j) = n + 2;
              }
            }
          }
        }
      }
    }
    band.updateBorder();
  }

#ifdef USE_MPI
  int nmixed_global = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
        nmixed_global++;
      }
    }
  }

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nmixed_local = std::max(nmixed_global / size, 1);
  std::vector<int> proc_offset(size + 1);
  proc_offset[0] = 0;
  for (int r = 0; r < size; r++)
    proc_offset[r + 1] = proc_offset[r] + nmixed_local;
  proc_offset[size] = nmixed_global;
  for (int r = 1; r < size + 1; r++)
    proc_offset[r] = std::min(proc_offset[r], nmixed_global);
  nmixed_local = proc_offset[rank + 1] - proc_offset[rank];

  const int size_moments = 2 * 6;
  std::vector<double> mom_update_local(size_moments * nmixed_local);
  std::vector<double> mom_update_global(size_moments * nmixed_global);
#endif

  // Compute fluxes
  std::vector<IRL2D::BezierList> fluxes, clipped_fluxes;
  Data<IRL2D::BezierList> pre_images(&mesh);
  Data<IRL2D::Moments> liq_mom_update(&mesh);
  Data<IRL2D::Moments> gas_mom_update(&mesh);

  int count = 0, count_local = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      const double cell_volume = mesh.cell_volume();
      // Initialize fluxes to zero
      (liq_mom_update)(i, j) = IRL2D::Moments();
      (gas_mom_update)(i, j) = IRL2D::Moments();

      // Only solve advection in narrow band near the interface
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
#ifdef USE_MPI
        if (count >= proc_offset[rank] && count < proc_offset[rank + 1]) {
#endif
          const auto x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
          const auto x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
          const auto cell = IRL2D::RectangleFromBounds(x0, x1);

          // Initialize face fluxes
          std::array<double, 4> exact_fluxes;
          for (int n = 0; n < 4; n++) {
            exact_fluxes[n] =
                IRL2D::IntegrateFlux(cell[n].first, cell[(n + 1) % 4].first,
                                     a_dt, a_time, getExactVelocity2D);
          }
          const auto preimage =
              CreatePreImage(x0, x1, -a_dt, a_time + a_dt, getExactVelocity2D,
                             getExactGradient2D, correct_fluxes, exact_fluxes);

          fluxes.push_back(preimage);

          // Compute flux bounding boxes
          const auto bbox = IRL2D::BoundingBox(preimage);
          static int idmin[2], idmax[2];
          mesh.getIndices(bbox.first, idmin);
          mesh.getIndices(bbox.second, idmax);

          // Compute liquid fluxes by intersection on a n x n neighborhood
          // for (int ii = std::max(idmin[0], mesh.imino());
          //      ii <= std::min(idmax[0], mesh.imaxo()); ++ii) {
          //   for (int jj = std::max(idmin[1], mesh.imino());
          //        jj <= std::min(idmax[1], mesh.imaxo()); ++jj) {
          for (int ii = i - nlayers; ii <= i + nlayers; ii++) {
            for (int jj = j - nlayers; jj <= j + nlayers; jj++) {
              const auto xn0 = IRL2D::Vec(mesh.x(ii), mesh.y(jj));
              const auto xn1 = IRL2D::Vec(mesh.x(ii + 1), mesh.y(jj + 1));
              (liq_mom_update)(i, j) += IRL2D::ComputeMoments(
                  preimage, xn0, xn1, (*a_interface)(ii, jj));
              clipped_fluxes.push_back(ClipByRectangleAndParabola(
                  preimage, xn0, xn1, (*a_interface)(ii, jj)));
            }
          }

          // (liq_mom_update)(i, j).m0() =
          //     std::clamp((liq_mom_update)(i, j).m0(), 0.0, cell_volume);

#ifndef USE_MPI
          // Update gas fluxes
          (gas_mom_update)(i, j) +=
              IRL2D::ComputeMoments(preimage) - (liq_mom_update)(i, j);
#else
        const auto fm = IRL2D::ComputeMoments(preimage);
        mom_update_local[count_local++] = liq_mom_update(i, j).m0();
        mom_update_local[count_local++] = liq_mom_update(i, j).m1()[0];
        mom_update_local[count_local++] = liq_mom_update(i, j).m1()[1];
        mom_update_local[count_local++] = liq_mom_update(i, j).m2()[0][0];
        mom_update_local[count_local++] = liq_mom_update(i, j).m2()[1][0];
        mom_update_local[count_local++] = liq_mom_update(i, j).m2()[1][1];
        mom_update_local[count_local++] = fm.m0() - liq_mom_update(i, j).m0();
        mom_update_local[count_local++] =
            fm.m1()[0] - liq_mom_update(i, j).m1()[0];
        mom_update_local[count_local++] =
            fm.m1()[1] - liq_mom_update(i, j).m1()[1];
        mom_update_local[count_local++] =
            fm.m2()[0][0] - liq_mom_update(i, j).m2()[0][0];
        mom_update_local[count_local++] =
            fm.m2()[1][0] - liq_mom_update(i, j).m2()[1][0];
        mom_update_local[count_local++] =
            fm.m2()[1][1] - liq_mom_update(i, j).m2()[1][1];
#endif
#ifdef USE_MPI
        }
        count++;
#endif
      }
    }
  }
  // if (fluxes.size() > 0) {
  //   IRL2D::ToVTK(fluxes, "fluxes");
  //   IRL2D::ToVTK(clipped_fluxes, "clipped_fluxes");
  // }

#ifdef USE_MPI
  std::vector<int> proc_count(size);
  for (int r = 0; r < size; r++) {
    proc_count[r] = size_moments * (proc_offset[r + 1] - proc_offset[r]);
    proc_offset[r] = size_moments * proc_offset[r];
  }
  MPI_Allgatherv(mom_update_local.data(), mom_update_local.size(), MPI_DOUBLE,
                 mom_update_global.data(), proc_count.data(),
                 proc_offset.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  count_local = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
        (liq_mom_update)(i, j).m0() = mom_update_global[count_local++];
        (liq_mom_update)(i, j).m1()[0] = mom_update_global[count_local++];
        (liq_mom_update)(i, j).m1()[1] = mom_update_global[count_local++];
        (liq_mom_update)(i, j).m2()[0][0] = mom_update_global[count_local++];
        (liq_mom_update)(i, j).m2()[1][0] = mom_update_global[count_local++];
        (liq_mom_update)(i, j).m2()[1][1] = mom_update_global[count_local++];
        (liq_mom_update)(i, j).m2()[0][1] = (liq_mom_update)(i, j).m2()[1][0];
        (gas_mom_update)(i, j).m0() = mom_update_global[count_local++];
        (gas_mom_update)(i, j).m1()[0] = mom_update_global[count_local++];
        (gas_mom_update)(i, j).m1()[1] = mom_update_global[count_local++];
        (gas_mom_update)(i, j).m2()[0][0] = mom_update_global[count_local++];
        (gas_mom_update)(i, j).m2()[1][0] = mom_update_global[count_local++];
        (gas_mom_update)(i, j).m2()[1][1] = mom_update_global[count_local++];
        (gas_mom_update)(i, j).m2()[0][1] = (gas_mom_update)(i, j).m2()[1][0];
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  gas_mom_update.updateBorder();
  liq_mom_update.updateBorder();

  // Now calculate VOF from the face fluxes.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double cell_volume = mesh.cell_volume();
      const auto cell_centroid = IRL2D::Vec(mesh.xm(i), mesh.ym(j));
      const double liquid_volume_fraction =
          (*a_liquid_moments)(i, j).m0() / cell_volume;
      if (band(i, j) > 0) {
        // Compute moments in transported cell image
        (*a_liquid_moments)(i, j) = (liq_mom_update)(i, j);
        (*a_gas_moments)(i, j) = (gas_mom_update)(i, j);

        const std::array<Data<IRL2D::Moments>*, 2> moment_list(
            {a_liquid_moments, a_gas_moments});

        for (int m = 0; m < 2; m++) {
          const auto moments = moment_list[m];
          const auto M0 = (*moments)(i, j).m0();
          const auto M1 = (*moments)(i, j).m1();
          const auto M2 = (*moments)(i, j).m2();

          // Correct moment 0
          double M0_final = M0;  // std::clamp(M0, 0.0, cell_volume);

          auto x0_k1 = M1 / IRL::safelyEpsilon(M0);
          const auto u0_k1 = getExactVelocity2D(a_time, x0_k1);
          const auto gradu0_k1 = getExactGradient2D(a_time, x0_k1);
          const auto x0_k2 = x0_k1 + a_dt * u0_k1 * 0.5;
          const auto u0_k2 = getExactVelocity2D(a_time + 0.5 * a_dt, x0_k2);
          const auto gradu0_k2 = getExactGradient2D(a_time + 0.5 * a_dt, x0_k2);
          const auto x0_k3 = x0_k1 + a_dt * u0_k2 * 0.5;
          const auto u0_k3 = getExactVelocity2D(a_time + 0.5 * a_dt, x0_k3);
          const auto gradu0_k3 = getExactGradient2D(a_time + 0.5 * a_dt, x0_k3);
          const auto x0_k4 = x0_k1 + a_dt * u0_k3;
          const auto u0_k4 = getExactVelocity2D(a_time + a_dt, x0_k4);
          const auto gradu0_k4 = getExactGradient2D(a_time + a_dt, x0_k4);

          const auto M1_final =
              M1 + a_dt * M0_final *
                       (u0_k1 + 2.0 * u0_k2 + 2.0 * u0_k3 + u0_k4) / 6.0;

          const auto M2t0_k1 = M0_final * IRL2D::outer_product(x0_k1, u0_k1);
          const auto M2t0_k2 = M0_final * IRL2D::outer_product(x0_k2, u0_k2);
          const auto M2t0_k3 = M0_final * IRL2D::outer_product(x0_k3, u0_k3);
          const auto M2t0_k4 = M0_final * IRL2D::outer_product(x0_k4, u0_k4);
          const auto M2t1_k1 = -M0_final * IRL2D::outer_product(x0_k1, x0_k1) *
                               gradu0_k1.transpose();
          const auto M2t1_k2 = -M0_final * IRL2D::outer_product(x0_k2, x0_k2) *
                               gradu0_k2.transpose();
          const auto M2t1_k3 = -M0_final * IRL2D::outer_product(x0_k3, x0_k3) *
                               gradu0_k3.transpose();
          const auto M2t1_k4 = -M0_final * IRL2D::outer_product(x0_k4, x0_k4) *
                               gradu0_k4.transpose();
          const auto M2t2_k1 = M2 * gradu0_k1.transpose();
          const auto M2t2_k2 =
              (M2 + a_dt * 0.5 *
                        (M2t0_k1 + M2t0_k1.transpose() + M2t1_k1 +
                         M2t1_k1.transpose() + M2t2_k1 + M2t2_k1.transpose())) *
              gradu0_k2.transpose();
          const auto M2t2_k3 =
              (M2 + a_dt * 0.5 *
                        (M2t0_k2 + M2t0_k2.transpose() + M2t1_k2 +
                         M2t1_k2.transpose() + M2t2_k2 + M2t2_k2.transpose())) *
              gradu0_k3.transpose();
          const auto M2t2_k4 = (M2 + a_dt * (M2t0_k3 + M2t0_k3.transpose() +
                                             M2t1_k3 + M2t1_k3.transpose() +
                                             M2t2_k3 + M2t2_k3.transpose())) *
                               gradu0_k4.transpose();
          const auto M2_rk4 =
              a_dt *
              ((M2t0_k1 + 2.0 * M2t0_k2 + 2.0 * M2t0_k3 + M2t0_k4) +
               (M2t1_k1 + 2.0 * M2t1_k2 + 2.0 * M2t1_k3 + M2t1_k4) +
               (M2t2_k1 + 2.0 * M2t2_k2 + 2.0 * M2t2_k3 + M2t2_k4)) *
              (1.0 / 6.0);
          const auto M2_final = M2 + M2_rk4 + M2_rk4.transpose();

          (*moments)(i, j).m0() = M0_final;
          (*moments)(i, j).m1() = M1_final;
          (*moments)(i, j).m2() = M2_final;
        }
      }
    }
  }

  const double cell_volume = mesh.cell_volume();
  std::vector<std::tuple<double, int, int>> sorted_volumes;
  sorted_volumes.resize(mesh.getNx() * mesh.getNy());
  count = 0;
  double deficit = 0.0, extra = 0.0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      sorted_volumes[count++] =
          std::make_tuple((*a_liquid_moments)(i, j).m0(), i, j);
      if ((*a_liquid_moments)(i, j).m0() < 0.0) {
        deficit -= (*a_liquid_moments)(i, j).m0();
      } else if ((*a_liquid_moments)(i, j).m0() > cell_volume) {
        extra += (*a_liquid_moments)(i, j).m0() - cell_volume;
      }
    }
  }
  std::sort(sorted_volumes.begin(), sorted_volumes.end());

  int nzero = -1, none = -1;
  for (int n = 0; n < mesh.getNx() * mesh.getNy() - 1; n++) {
    const auto [volume0, i0, j0] = sorted_volumes[n];
    const auto [volume1, i1, j1] = sorted_volumes[n + 1];
    if (volume0 <= 0.0 && volume1 > 0.0) {
      nzero = n;
    } else if (volume0 <= cell_volume && volume1 > cell_volume) {
      none = n;
    }
  }

  if (deficit > 0.0) {
    for (int n = 0; n <= nzero; n++) {
      const auto [volume, i, j] = sorted_volumes[n];
      (*a_liquid_moments)(i, j).m0() = 0.0;
    }
    for (int n = nzero + 1; n <= none; n++) {
      const auto [volume, i, j] = sorted_volumes[n];
      if (volume < deficit) {
        (*a_liquid_moments)(i, j).m0() = 0.0;
        deficit -= volume;
      } else {
        (*a_liquid_moments)(i, j).m0() -= deficit;
        deficit = 0.0;
        break;
      }
    }
  }

  if (extra > 0.0) {
    for (int n = none + 1; n < mesh.getNx() * mesh.getNy(); n++) {
      const auto [volume, i, j] = sorted_volumes[n];
      (*a_liquid_moments)(i, j).m0() = cell_volume;
    }
    for (int n = none; n > nzero; n--) {
      const auto [volume, i, j] = sorted_volumes[n];
      if (cell_volume < volume + extra) {
        (*a_liquid_moments)(i, j).m0() = cell_volume;
        extra -= cell_volume - volume;
      } else {
        (*a_liquid_moments)(i, j).m0() += extra;
        extra = 0.0;
        break;
      }
    }
  }

  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();
  correctMomentLocations(a_liquid_moments);
  correctMomentLocations(a_gas_moments);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

void ReSyFullLagL::advectVOF(const std::string& a_simulation_type,
                        const std::string& a_advection_method,
                        const std::string& a_reconstruction_method,
                        const double a_dt, const double a_time,
                        const Data<double>& a_U, const Data<double>& a_V,
                        Data<IRL2D::Moments>* a_liquid_moments,
                        Data<IRL2D::Moments>* a_gas_moments,
                        Data<IRL2D::Parabola>* a_interface) {
  const double tol = 100. * std::numeric_limits<double>::epsilon();
  const bool correct_fluxes = true;
  const BasicMesh& mesh = a_liquid_moments->getMesh();

  IRL2D::BezierList (*CreatePreImage)(
      const IRL2D::Vec& x0, const IRL2D::Vec& x1, const double dt,
      const double time,
      const IRL2D::Vec (*vel)(const double t, const IRL2D::Vec& P),
      const IRL2D::Mat (*grad_vel)(const double t, const IRL2D::Vec& P),
      const bool correct_area, const std::array<double, 4>& exact_area);

  CreatePreImage = IRL2D::CreateLinearPreImage;

  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      const auto cell =
          IRL2D::RectangleFromBounds(IRL2D::Vec(mesh.x(i), mesh.y(j)),
                                     IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
      const auto cell_moments = IRL2D::ComputeMoments(cell);
      auto moments = IRL2D::ComputeMoments(cell, (*a_interface)(i, j));
      moments.m0() = (*a_liquid_moments)(i, j).m0();
      (*a_liquid_moments)(i, j) = moments;
      (*a_gas_moments)(i, j) = cell_moments - (*a_liquid_moments)(i, j);
    }
  }
  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();

  const IRL2D::Vec (*getExactVelocity2D)(const double t, const IRL2D::Vec& P);
  const IRL2D::Mat (*getExactGradient2D)(const double t, const IRL2D::Vec& P);

  if (a_simulation_type == "Rotation2D") {
    getExactVelocity2D = Rotation2D::getExactVelocity2D;
    getExactGradient2D = Rotation2D::getExactVelocityGradient2D;
  } else if (a_simulation_type == "Oscillation2D") {
    getExactVelocity2D = Oscillation2D::getExactVelocity2D;
    getExactGradient2D = Oscillation2D::getExactVelocityGradient2D;
  } else if (a_simulation_type == "Deformation2D") {
    getExactVelocity2D = Deformation2D::getExactVelocity2D;
    getExactGradient2D = Deformation2D::getExactVelocityGradient2D;
  }

  // Calculate CFL
  double CFL = -DBL_MAX;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      CFL = std::fmax(CFL, std::fmax(a_U(i, j) * a_dt / mesh.dx(),
                                     a_V(i, j) * a_dt / mesh.dy()));
    }
  }

  // Generate band that is ceil(CFL) thick around interface
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      band(i, j) = 0;
      const double liquid_volume_fraction =
          (*a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }
  band.updateBorder();
  const int nlayers = 1 + static_cast<int>(std::ceil(CFL));
  for (int n = 0; n < nlayers; ++n) {
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        if (band(i, j) == 0) {
          for (int ii = -1; ii <= 1; ++ii) {
            for (int jj = -1; jj <= 1; ++jj) {
              if (band(i + ii, j + jj) == n + 1) {
                band(i, j) = n + 2;
              }
            }
          }
        }
      }
    }
    band.updateBorder();
  }

#ifdef USE_MPI
  int nmixed_global = 0;
  int n_simplices = 1;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
        nmixed_global++;
      }
    }
  }

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nmixed_local = std::max(nmixed_global / size, 1);
  std::vector<int> proc_offset(size + 1);
  proc_offset[0] = 0;
  for (int r = 0; r < size; r++)
    proc_offset[r + 1] = proc_offset[r] + nmixed_local;
  proc_offset[size] = nmixed_global;
  for (int r = 1; r < size + 1; r++)
    proc_offset[r] = std::min(proc_offset[r], nmixed_global) * n_simplices;
  nmixed_local = proc_offset[rank + 1] - proc_offset[rank];

  const int size_moments = 2 * 6;
  std::vector<double> mom_update_local(size_moments * nmixed_local);
  std::vector<double> mom_update_global(size_moments * n_simplices * nmixed_global);
#endif

  // Compute fluxes
  Data<IRL2D::BezierList> pre_images(&mesh);
  Data<IRL2D::Moments> liq_mom_final(&mesh);
  Data<IRL2D::Moments> gas_mom_final(&mesh);
  Data<std::vector<IRL2D::Moments>> tri_liq_mom_update(&mesh);
  Data<std::vector<IRL2D::Moments>> tri_gas_mom_update(&mesh);
  Data<std::vector<IRL2D::Mat>> mapping_A(&mesh);
  Data<std::vector<IRL2D::Vec>> mapping_b(&mesh);
  Data<std::vector<IRL2D::Vec>> mapped_m1_triangle(&mesh);
  int num_triangles = 8;

  int count = 0, count_local = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      const double cell_volume = mesh.cell_volume();
      // initializing
      tri_liq_mom_update(i,j) = std::vector<IRL2D::Moments>(num_triangles, IRL2D::Moments());
      tri_gas_mom_update(i,j) = std::vector<IRL2D::Moments>(num_triangles, IRL2D::Moments());
      liq_mom_final(i,j) = IRL2D::Moments();
      gas_mom_final(i,j) = IRL2D::Moments();
      mapping_A(i,j) = std::vector<IRL2D::Mat>(num_triangles, IRL2D::Mat());
      mapping_b(i,j) = std::vector<IRL2D::Vec>(num_triangles, IRL2D::Vec());
      mapped_m1_triangle(i,j) = std::vector<IRL2D::Vec>(num_triangles, IRL2D::Vec());

      // Only solve advection in narrow band near the interface
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
#ifdef USE_MPI
        if (count >= proc_offset[rank] && count < proc_offset[rank + 1]) {
#endif
          const auto x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
          const auto x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
          const auto cell = IRL2D::RectangleFromBounds(x0, x1);

          // Initialize face fluxes
          std::array<double, 4> exact_fluxes;
          for (int n = 0; n < 4; n++) {
            exact_fluxes[n] =
                IRL2D::IntegrateFlux(cell[n].first, cell[(n + 1) % 4].first,
                                     a_dt, a_time, getExactVelocity2D);
          }
          const auto preimage =
              CreatePreImage(x0, x1, -a_dt, a_time + a_dt, getExactVelocity2D,
                             getExactGradient2D, correct_fluxes, exact_fluxes);

          // triangulation of cell and preimage
          std::vector<IRL2D::BezierList> triangulated_cell = IRL2D::TriangulateCell(cell, false);
          std::vector<IRL2D::BezierList> triangulated_preimage = IRL2D::TriangulateCell(preimage, true);

          for (int t = 0; t < num_triangles; ++t){
            // mapping from preimage to cell
            std::pair<IRL2D::Mat, IRL2D::Vec> mapping_MatVec 
            = IRL2D::MappingMatVec(triangulated_preimage[t], triangulated_cell[t]);
            mapping_A(i,j)[t] = mapping_MatVec.first; mapping_b(i,j)[t] = mapping_MatVec.second;
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii){
              for (int jj = j - nlayers; jj <= j + nlayers; ++jj){
                const auto xn0 = IRL2D::Vec(mesh.x(ii), mesh.y(jj));
                const auto xn1 = IRL2D::Vec(mesh.x(ii + 1), mesh.y(jj + 1));
                tri_liq_mom_update(i,j)[t] += IRL2D::ComputeMoments(triangulated_preimage[t], 
                                              xn0, xn1, (*a_interface)(ii, jj));
              }
            }
            tri_gas_mom_update(i,j)[t] = IRL2D::ComputeMoments(triangulated_preimage[t])
                                         - (tri_liq_mom_update)(i, j)[t];
            liq_mom_final(i,j) += IRL2D::ComputeMappedTriangleMoments(tri_liq_mom_update(i,j)[t],
                                  mapping_A(i,j)[t], mapping_b(i,j)[t]);
            gas_mom_final(i,j) += IRL2D::ComputeMappedTriangleMoments(tri_gas_mom_update(i,j)[t],
                                  mapping_A(i,j)[t], mapping_b(i,j)[t]);
          }

#ifdef USE_MPI
          const auto fm = IRL2D::ComputeMoments(preimage);
          mom_update_local[count_local++] = liq_mom_final(i,j).m0();
          mom_update_local[count_local++] = liq_mom_final(i,j).m1()[0];
          mom_update_local[count_local++] = liq_mom_final(i,j).m1()[1];
          mom_update_local[count_local++] = liq_mom_final(i,j).m2()[0][0];
          mom_update_local[count_local++] = liq_mom_final(i,j).m2()[1][0];
          mom_update_local[count_local++] = liq_mom_final(i,j).m2()[1][1];
          mom_update_local[count_local++] = gas_mom_final(i,j).m0();
          mom_update_local[count_local++] = gas_mom_final(i,j).m1()[0];
          mom_update_local[count_local++] = gas_mom_final(i,j).m1()[1];
          mom_update_local[count_local++] = gas_mom_final(i,j).m2()[0][0];
          mom_update_local[count_local++] = gas_mom_final(i,j).m2()[1][0];
          mom_update_local[count_local++] = gas_mom_final(i,j).m2()[1][1];
        }
        count++;
#endif
      }
    }
  } // end of calculating moments for each preimage

#ifdef USE_MPI
  std::vector<int> proc_count(size);
  for (int r = 0; r < size; r++) {
    proc_count[r] = size_moments * (proc_offset[r + 1] - proc_offset[r]);
    proc_offset[r] = size_moments * proc_offset[r];
  }
  MPI_Allgatherv(mom_update_local.data(), mom_update_local.size(), MPI_DOUBLE,
                 mom_update_global.data(), proc_count.data(),
                 proc_offset.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  count_local = 0;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
        (liq_mom_final)(i, j).m0() = mom_update_global[count_local++];
        (liq_mom_final)(i, j).m1()[0] = mom_update_global[count_local++];
        (liq_mom_final)(i, j).m1()[1] = mom_update_global[count_local++];
        (liq_mom_final)(i, j).m2()[0][0] = mom_update_global[count_local++];
        (liq_mom_final)(i, j).m2()[1][0] = mom_update_global[count_local++];
        (liq_mom_final)(i, j).m2()[1][1] = mom_update_global[count_local++];
        (liq_mom_final)(i, j).m2()[0][1] = (liq_mom_final)(i, j).m2()[1][0];
        (gas_mom_final)(i, j).m0() = mom_update_global[count_local++];
        (gas_mom_final)(i, j).m1()[0] = mom_update_global[count_local++];
        (gas_mom_final)(i, j).m1()[1] = mom_update_global[count_local++];
        (gas_mom_final)(i, j).m2()[0][0] = mom_update_global[count_local++];
        (gas_mom_final)(i, j).m2()[1][0] = mom_update_global[count_local++];
        (gas_mom_final)(i, j).m2()[1][1] = mom_update_global[count_local++];
        (gas_mom_final)(i, j).m2()[0][1] = (gas_mom_final)(i, j).m2()[1][0];
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  liq_mom_final.updateBorder();
  gas_mom_final.updateBorder();

  // Updating moments in the narrow band
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      if (band(i, j) > 0) {
        (*a_liquid_moments)(i,j) = liq_mom_final(i,j);
        (*a_gas_moments)(i,j) = gas_mom_final(i,j);
      }
    }
  }

  const double cell_volume = mesh.cell_volume();
  std::vector<std::tuple<double, int, int>> sorted_volumes;
  sorted_volumes.resize(mesh.getNx() * mesh.getNy());
  count = 0;
  double deficit = 0.0, extra = 0.0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      sorted_volumes[count++] =
          std::make_tuple((*a_liquid_moments)(i, j).m0(), i, j);
      if ((*a_liquid_moments)(i, j).m0() < 0.0) {
        deficit -= (*a_liquid_moments)(i, j).m0();
      } else if ((*a_liquid_moments)(i, j).m0() > cell_volume) {
        extra += (*a_liquid_moments)(i, j).m0() - cell_volume;
      }
    }
  }
  std::sort(sorted_volumes.begin(), sorted_volumes.end());

  int nzero = -1, none = -1;
  for (int n = 0; n < mesh.getNx() * mesh.getNy() - 1; n++) {
    const auto [volume0, i0, j0] = sorted_volumes[n];
    const auto [volume1, i1, j1] = sorted_volumes[n + 1];
    if (volume0 <= 0.0 && volume1 > 0.0) {
      nzero = n;
    } else if (volume0 <= cell_volume && volume1 > cell_volume) {
      none = n;
    }
  }

  if (deficit > 0.0) {
    for (int n = 0; n <= nzero; n++) {
      const auto [volume, i, j] = sorted_volumes[n];
      (*a_liquid_moments)(i, j).m0() = 0.0;
    }
    for (int n = nzero + 1; n <= none; n++) {
      const auto [volume, i, j] = sorted_volumes[n];
      if (volume < deficit) {
        (*a_liquid_moments)(i, j).m0() = 0.0;
        deficit -= volume;
      } else {
        (*a_liquid_moments)(i, j).m0() -= deficit;
        deficit = 0.0;
        break;
      }
    }
  }

  if (extra > 0.0) {
    for (int n = none + 1; n < mesh.getNx() * mesh.getNy(); n++) {
      const auto [volume, i, j] = sorted_volumes[n];
      (*a_liquid_moments)(i, j).m0() = cell_volume;
    }
    for (int n = none; n > nzero; n--) {
      const auto [volume, i, j] = sorted_volumes[n];
      if (cell_volume < volume + extra) {
        (*a_liquid_moments)(i, j).m0() = cell_volume;
        extra -= cell_volume - volume;
      } else {
        (*a_liquid_moments)(i, j).m0() += extra;
        extra = 0.0;
        break;
      }
    }
  }

  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();
  correctMomentLocations(a_liquid_moments);
  correctMomentLocations(a_gas_moments);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

void correctMomentLocations(Data<IRL2D::Moments>* a_liquid_moments) {
  const BasicMesh& mesh = (*a_liquid_moments).getMesh();
  // Fix distance to recreate volume fraction

  // Extract moment in tensor notation
  Data<IRL2D::Vec> Shift(&mesh);
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      Shift(i, j) = IRL2D::Vec();
    }
  }

  // x- boundary
  for (int i = mesh.imino(); i < mesh.imin(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      Shift(i, j)[0] -= mesh.lx();
    }
  }

  // x+ boundary
  for (int i = mesh.imax() + 1; i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      Shift(i, j)[0] += mesh.lx();
    }
  }

  // y- boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j < mesh.jmin(); ++j) {
      Shift(i, j)[1] -= mesh.ly();
    }
  }

  // y+ boundary
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmax() + 1; j <= mesh.jmaxo(); ++j) {
      Shift(i, j)[1] += mesh.ly();
    }
  }

  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      const auto M0 = (*a_liquid_moments)(i, j).m0();
      const auto M1 = (*a_liquid_moments)(i, j).m1();
      const auto M2 = (*a_liquid_moments)(i, j).m2();
      auto M1_final = M1 + M0 * Shift(i, j);
      auto M2_final = M2 + IRL2D::outer_product(M1, Shift(i, j)) +
                      IRL2D::outer_product(Shift(i, j), M1) +
                      M0 * IRL2D::outer_product(Shift(i, j), Shift(i, j));
      (*a_liquid_moments)(i, j).m1() = M1_final;
      (*a_liquid_moments)(i, j).m2() = M2_final;
    }
  }
}