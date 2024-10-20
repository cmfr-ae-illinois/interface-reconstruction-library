// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <mpi.h>
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
  const bool correct_fluxes = true;
  const BasicMesh& mesh = a_liquid_moments->getMesh();

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

  // setPhaseQuantities((*a_interface), a_liquid_moments, a_gas_moments);

  // Allocate storage for face fluxes
  Data<IRL2D::Moments> face_liquid_flux[2] = {Data<IRL2D::Moments>(&mesh),
                                              Data<IRL2D::Moments>(&mesh)};
  Data<IRL2D::Moments> face_gas_flux[2] = {Data<IRL2D::Moments>(&mesh),
                                           Data<IRL2D::Moments>(&mesh)};

  // Calculate CFL
  double CFL = -DBL_MAX;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      CFL = std::fmax(CFL, std::fmax(a_U(i, j) * a_dt / mesh.dx(),
                                     a_V(i, j) * a_dt / mesh.dy()));
    }
  }

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

  std::vector<IRL2D::BezierList> fluxes;
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
      // Initialize fluxes to zero
      for (int dim = 0; dim < 2; ++dim) {
        (face_liquid_flux[dim])(i, j) = IRL2D::Moments();
        (face_gas_flux[dim])(i, j) = IRL2D::Moments();
      }
      // Only solve advection in narrow band near the interface
      if (band(i, j) > 0 || band(i - 1, j) > 0 || band(i, j - 1) > 0) {
        const auto x0 = IRL2D::Vec(mesh.x(i), mesh.y(j));
        const auto x1 = IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1));
        const auto cell = IRL2D::RectangleFromBounds(x0, x1);

        // Initialize face fluxes
        IRL2D::BezierList face_cell[2];
        face_cell[0] = IRL2D::CreateFluxCell(
            cell[3].first, cell[0].first, -a_dt, a_time, getExactVelocity2D,
            getExactGradient2D, correct_fluxes,
            IntegrateFlux(cell[3].first, cell[0].first, a_dt, a_time,
                          getExactVelocity2D));
        face_cell[1] = IRL2D::CreateFluxCell(
            cell[0].first, cell[1].first, -a_dt, a_time, getExactVelocity2D,
            getExactGradient2D, correct_fluxes,
            IntegrateFlux(cell[0].first, cell[1].first, a_dt, a_time,
                          getExactVelocity2D));
        fluxes.push_back(face_cell[0]);
        fluxes.push_back(face_cell[1]);

        // Compute liquid fluxes by intersection on a n x n neighborhood
        for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            const auto xn0 = IRL2D::Vec(mesh.x(ii), mesh.y(jj));
            const auto xn1 = IRL2D::Vec(mesh.x(ii + 1), mesh.y(jj + 1));
            for (int dim = 0; dim < 2; ++dim) {
              (face_liquid_flux[dim])(i, j) += IRL2D::ComputeMoments(
                  face_cell[dim], xn0, xn1, (*a_interface)(ii, jj));
            }
          }
        }

        // Update gas fluxes
        for (int dim = 0; dim < 2; ++dim) {
          (face_gas_flux[dim])(i, j) += IRL2D::ComputeMoments(face_cell[dim]) -
                                        (face_liquid_flux[dim])(i, j);
        }
      }
    }
  }

  IRL2D::ToVTK(fluxes, "fluxes");

  face_liquid_flux[0].updateBorder();
  face_liquid_flux[1].updateBorder();
  face_gas_flux[0].updateBorder();
  face_gas_flux[1].updateBorder();

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
            (*a_liquid_moments)(i, j) + (face_liquid_flux[0])(i, j) -
            (face_liquid_flux[0])(i + 1, j) + (face_liquid_flux[1])(i, j) -
            (face_liquid_flux[1])(i, j + 1);
        (*a_gas_moments)(i, j) =
            (*a_gas_moments)(i, j) + (face_gas_flux[0])(i, j) -
            (face_gas_flux[0])(i + 1, j) + (face_gas_flux[1])(i, j) -
            (face_gas_flux[1])(i, j + 1);

        const std::array<Data<IRL2D::Moments>*, 2> moment_list(
            {a_liquid_moments, a_gas_moments});

        for (int m = 0; m < 2; m++) {
          const auto moments = moment_list[m];
          const auto M0 = (*moments)(i, j).m0();
          const auto M1 = (*moments)(i, j).m1();
          const auto M2 = (*moments)(i, j).m2();

          // Correct moment 0
          double M0_final = M0;
          if (M0 < 0.0) {
            M0_final = 0.0;
          } else if (M0 > cell_volume) {
            M0_final = cell_volume;
          }
          (*moments)(i, j).m0() = M0_final;

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
      } else if (liquid_volume_fraction < IRL::global_constants::VF_LOW ||
                 liquid_volume_fraction > IRL::global_constants::VF_HIGH) {
        const double clipped_vfrac =
            std::max(0.0, std::min(1.0, liquid_volume_fraction));
        const auto cell_moments =
            IRL2D::ComputeMoments(IRL2D::RectangleFromBounds(
                IRL2D::Vec(mesh.x(i), mesh.y(j)),
                IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1))));
        (*a_liquid_moments)(i, j) = clipped_vfrac * cell_moments;
        (*a_gas_moments)(i, j) = (1.0 - clipped_vfrac) * cell_moments;
      }
    }
  }
  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();
  correctMomentLocations(a_liquid_moments);
  correctMomentLocations(a_gas_moments);
}

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