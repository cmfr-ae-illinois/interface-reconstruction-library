// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/2d_advector/solver.h"

#include <stdio.h>
#include <cstdio>
#include <string>

#include "examples/2d_advector/basic_mesh.h"
#include "examples/2d_advector/data.h"
#include "examples/2d_advector/rotation_2d.h"

void setPhaseQuantities(const Data<IRL2D::Parabola>& a_interface,
                        Data<IRL2D::Moments>* a_liquid_moments,
                        Data<IRL2D::Moments>* a_gas_moments,
                        Data<double>* a_vfrac) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      const auto cell =
          IRL2D::RectangleFromBounds(IRL2D::Vec(mesh.x(i), mesh.y(j)),
                                     IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
      const auto cell_moments = IRL2D::ComputeMoments(cell);
      (*a_liquid_moments)(i, j) =
          IRL2D::ComputeMoments(cell, a_interface(i, j));
      (*a_gas_moments)(i, j) = cell_moments - (*a_liquid_moments)(i, j);
      (*a_vfrac)(i, j) = (*a_liquid_moments)(i, j).m0() / mesh.cell_volume();
    }
  }
  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();
  a_vfrac->updateBorder();
}

void writeDiagnosticsHeader(void) {
  printf("%10s %20s %12s %20s %20s %20s %20s %20s %20s %20s\n", "Iteration",
         "Time", "CFL", "liquidVFSum", "liquidVolSum", "ChangeLiquidVFSum",
         "ChangeLiquidVolSum", "AdvectionDuration", "ReconDuration",
         "OutputDuration", "InterfaceCells");
}

void writeOutDiagnostics(const int a_iteration, const double a_dt,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<IRL2D::Moments>& a_liquid_moments,
                         const Data<IRL2D::Parabola>& a_interface,
                         std::chrono::duration<double> a_VOF_duration,
                         std::chrono::duration<double> a_recon_duration,
                         std::chrono::duration<double> a_write_duration) {
  const BasicMesh& mesh = a_U.getMesh();
  static double initial_liquid_volume_fraction_sum;
  static double initial_liquid_volume_sum;
  // Calculate CFL
  double CFL = -DBL_MAX;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      CFL = std::fmax(CFL, std::fmax(a_U(i, j) * a_dt / mesh.dx(),
                                     a_V(i, j) * a_dt / mesh.dy()));
    }
  }

  // Calculate sum of volume fraction and sum of liquid volume
  double liquid_volume_fraction_sum = 0.0;
  double liquid_volume_sum = 0.0;
  int number_of_interface_cells = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          a_liquid_moments(i, j).m0() / mesh.cell_volume();
      liquid_volume_fraction_sum += liquid_volume_fraction;
      liquid_volume_sum += a_liquid_moments(i, j).m0();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        ++number_of_interface_cells;
      }
    }
  }
  // Save initial values to compare against.
  if (a_iteration == 0) {
    initial_liquid_volume_fraction_sum = liquid_volume_fraction_sum;
    initial_liquid_volume_sum = liquid_volume_sum;
  }
  printf(
      "%10d %20.4E %12.3F %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20d"
      "\n",
      a_iteration, a_simulation_time, CFL, liquid_volume_fraction_sum,
      liquid_volume_sum,
      liquid_volume_fraction_sum - initial_liquid_volume_fraction_sum,
      liquid_volume_sum - initial_liquid_volume_sum, a_VOF_duration.count(),
      a_recon_duration.count(), a_write_duration.count(),
      number_of_interface_cells);
}

// void printError(const BasicMesh& mesh,
//                 const Data<IRL::GeneralMoments3D<2>>& liquid_moments,
//                 const Data<IRL::GeneralMoments3D<2>>&
//                 starting_liquid_moments) {
//   double linf_error_m0 = 0.0;
//   double linf_error_m1 = 0.0;
//   double linf_error_m2 = 0.0;
//   double l1_error_m0 = 0.0;
//   double l1_error_m1 = 0.0;
//   double l1_error_m2 = 0.0;
//   double l2_error_m0 = 0.0;
//   double l2_error_m1 = 0.0;
//   double l2_error_m2 = 0.0;
//   double scale_m0 = 1.0 / std::pow(mesh.dx(), 3.0);
//   double scale_m1 = 1.0 / std::pow(mesh.dx(), 4.0);
//   double scale_m2 = 1.0 / std::pow(mesh.dx(), 5.0);
//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
//         const double liquid_volume_fraction =
//             liquid_moments(i, j, k)[0] / mesh.cell_volume();
//         if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
//             liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
//           auto mom_err =
//               (liquid_moments(i, j, k) - starting_liquid_moments(i, j, k));
//           linf_error_m0 = std::max(linf_error_m0, std::abs(mom_err[0]));
//           linf_error_m1 = std::max(linf_error_m1, std::abs(mom_err[1]));
//           linf_error_m1 = std::max(linf_error_m1, std::abs(mom_err[2]));
//           linf_error_m1 = std::max(linf_error_m1, std::abs(mom_err[3]));
//           linf_error_m2 = std::max(linf_error_m2, std::abs(mom_err[4]));
//           linf_error_m2 = std::max(linf_error_m2, std::abs(2.0 *
//           mom_err[5])); linf_error_m2 = std::max(linf_error_m2, std::abs(2.0
//           * mom_err[6])); linf_error_m2 = std::max(linf_error_m2,
//           std::abs(mom_err[7])); linf_error_m2 = std::max(linf_error_m2,
//           std::abs(2.0 * mom_err[8])); linf_error_m2 =
//           std::max(linf_error_m2, std::abs(mom_err[9])); l1_error_m0 +=
//           std::abs(mom_err[0]); l1_error_m1 += std::abs(mom_err[1]);
//           l1_error_m1 += std::abs(mom_err[2]);
//           l1_error_m1 += std::abs(mom_err[3]);
//           l1_error_m2 += std::abs(mom_err[4]);
//           l1_error_m2 += std::abs(2.0 * mom_err[5]);
//           l1_error_m2 += std::abs(2.0 * mom_err[6]);
//           l1_error_m2 += std::abs(mom_err[7]);
//           l1_error_m2 += std::abs(2.0 * mom_err[8]);
//           l1_error_m2 += std::abs(mom_err[9]);
//           l2_error_m0 += mom_err[0] * mom_err[0];
//           l2_error_m1 += mom_err[1] * mom_err[1];
//           l2_error_m1 += mom_err[2] * mom_err[2];
//           l2_error_m1 += mom_err[3] * mom_err[3];
//           l2_error_m2 += mom_err[4] * mom_err[4];
//           l2_error_m2 += 4.0 * mom_err[5] * mom_err[5];
//           l2_error_m2 += 4.0 * mom_err[6] * mom_err[6];
//           l2_error_m2 += mom_err[7] * mom_err[7];
//           l2_error_m2 += 4.0 * mom_err[8] * mom_err[8];
//           l2_error_m2 += mom_err[9] * mom_err[9];
//         }
//       }
//     }
//   }
//   l1_error_m0 /=
//       (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
//   l1_error_m1 /=
//       (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
//   l1_error_m2 /=
//       (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
//   l2_error_m0 /=
//       (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
//   l2_error_m1 /=
//       (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
//   l2_error_m2 /=
//       (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
//   linf_error_m0 *= scale_m0;
//   linf_error_m1 *= scale_m1;
//   linf_error_m2 *= scale_m2;
//   l1_error_m0 *= scale_m0;
//   l1_error_m1 *= scale_m1;
//   l1_error_m2 *= scale_m2;
//   l2_error_m0 = std::sqrt(l2_error_m0) * scale_m0;
//   l2_error_m1 = std::sqrt(l2_error_m1) * scale_m1;
//   l2_error_m2 = std::sqrt(l2_error_m2) * scale_m2;
//   std::cout << std::scientific << std::setprecision(2)
//             << "Linf M0 = " << linf_error_m0 << std::endl;
//   std::cout << "Linf M1 = " << linf_error_m1 << std::endl;
//   std::cout << "Linf M2 = " << linf_error_m2 << std::endl;
//   std::cout << "L1   M0 = " << l1_error_m0 << std::endl;
//   std::cout << "L1   M1 = " << l1_error_m1 << std::endl;
//   std::cout << "L1   M2 = " << l1_error_m2 << std::endl;
//   std::cout << "L2   M0 = " << l2_error_m0 << std::endl;
//   std::cout << "L2   M1 = " << l2_error_m1 << std::endl;
//   std::cout << "L2   M2 = " << l2_error_m2 << std::endl;
// }
