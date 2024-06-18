// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/paraboloid_advector/solver.h"

#include <stdio.h>
#include <cstdio>
#include <string>
#include "irl/paraboloid_reconstruction/parametrized_surface.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/planar_reconstruction/planar_localizer.h"

#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/deformation_3d.h"
#include "examples/paraboloid_advector/rotation_3d.h"
#include "examples/paraboloid_advector/translation_3d.h"

// Convert and store the mesh cells into localizers.
void initializeLocalizers(Data<IRL::PlanarLocalizer>* a_localizers) {
  // For each cell in the domain, construct the cell as a RectangularCuboid
  // and get the localizer.
  const BasicMesh& mesh = a_localizers->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt lower_corner(mesh.x(i), mesh.y(j), mesh.z(k));
        IRL::Pt upper_corner(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
        (*a_localizers)(i, j, k) =
            IRL::RectangularCuboid::fromBoundingPts(lower_corner, upper_corner)
                .getLocalizer();
      }
    }
  }
}

void initializeLocalizedParaboloids(
    const Data<IRL::PlanarLocalizer>& a_cell_localizers,
    const Data<IRL::Paraboloid>& a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>*
        a_linked_localized_paraboloids) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_linked_localized_paraboloids)(i, j, k) =
            IRL::LocalizedParaboloidLink<double>(&a_cell_localizers(i, j, k),
                                                 &a_interface(i, j, k));
      }
    }
  }
}

void setPhaseQuantities(const Data<IRL::Paraboloid>& a_interface,
                        Data<IRL::GeneralMoments3D<2>>* a_liquid_moments,
                        Data<IRL::GeneralMoments3D<2>>* a_gas_moments) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const auto moments = IRL::getVolumeMoments<
            IRL::SeparatedMoments<IRL::GeneralMoments3D<2>>>(
            cell, a_interface(i, j, k));
        (*a_liquid_moments)(i, j, k) = moments[0];
        (*a_gas_moments)(i, j, k) = moments[1];
      }
    }
  }
  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();
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
                         const Data<double>& a_W,
                         const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
                         const Data<IRL::Paraboloid>& a_interface,
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
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        CFL = std::fmax(CFL,
                        std::fmax(a_U(i, j, k) * a_dt / mesh.dx(),
                                  std::fmax(a_V(i, j, k) * a_dt / mesh.dy(),
                                            a_W(i, j, k) * a_dt / mesh.dz())));
      }
    }
  }
  // Calculate sum of volume fraction and sum of liquid volume
  double liquid_volume_fraction_sum = 0.0;
  double liquid_volume_sum = 0.0;
  int number_of_interface_cells = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        liquid_volume_fraction_sum += liquid_volume_fraction;
        liquid_volume_sum += a_liquid_moments(i, j, k)[0];
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          ++number_of_interface_cells;
        }
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

void writeInterfaceToFile(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::Paraboloid>& a_liquid_gas_interface, const double a_time,
    VTKOutput* a_output, const bool print) {
  const BasicMesh& mesh = a_liquid_moments.getMesh();

  std::vector<IRL::ParametrizedSurfaceOutput> surfaces;
  double radius = 0.25;
  double total_surface = 0.0, avg_mean_curv = 0.0;
  double exact_curv = 2.0 / radius;
  double max_mean_curv_error = 0.0;
  double l2_mean_curv_error = 0.0, l1_mean_curv_error = 0.0, l2_counter = 0.0;

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int split_proc = static_cast<int>(std::cbrt(static_cast<double>(size)));
  int imin = 0, nx = mesh.getNx(), jmin = 0, ny = mesh.getNy(), kmin = 0,
      nz = mesh.getNz();

  if (size > 1) {
    if (size == split_proc * split_proc * split_proc) {
      for (int i = 0; i < split_proc; i++) {
        for (int j = 0; j < split_proc; j++) {
          for (int k = 0; k < split_proc; k++) {
            if (i + split_proc * j + split_proc * split_proc * k == rank) {
              imin = i * (mesh.getNx() / split_proc);
              nx = std::min((i + 1) * (mesh.getNx() / split_proc),
                            mesh.getNx()) -
                   imin;
              jmin = j * (mesh.getNy() / split_proc);
              ny = std::min((j + 1) * (mesh.getNy() / split_proc),
                            mesh.getNy()) -
                   jmin;
              kmin = k * (mesh.getNz() / split_proc);
              nz = std::min((k + 1) * (mesh.getNz() / split_proc),
                            mesh.getNz()) -
                   kmin;
            }
          }
        }
      }
    } else {
      imin = rank * (mesh.getNx() / size);
      nx = std::min((rank + 1) * (mesh.getNx() / size), mesh.getNx()) - imin;
      jmin = 0;
      ny = mesh.getNy();
      kmin = 0;
      nz = mesh.getNz();
    }
  }

  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
        const double liquid_volume_fraction =
            a_liquid_moments(i, j, k)[0] / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          auto volume_and_surface = IRL::getVolumeMoments<IRL::AddSurfaceOutput<
              IRL::Volume, IRL::ParametrizedSurfaceOutput>>(
              cell, a_liquid_gas_interface(i, j, k));
          if (volume_and_surface.getMoments() > -DBL_MAX) {
            auto surface = volume_and_surface.getSurface();
            double length_scale = std::min(
                std::pow(cell.calculateVolume(), 1.0 / 3.0) / 3.0, 1.0e-2);
            surface.setLengthScale(length_scale);
            if (surface.getSurfaceArea() >
                1.0e-6 * length_scale * length_scale) {
              surfaces.push_back(surface);
            }
            total_surface += surface.getSurfaceArea();
            avg_mean_curv += surface.getMeanCurvatureIntegral();
            max_mean_curv_error = std::max(
                max_mean_curv_error,
                std::abs(surface.getAverageMeanCurvature() - exact_curv));
            // }
            l1_mean_curv_error +=
                std::abs(surface.getAverageMeanCurvature() - exact_curv);
            l2_mean_curv_error +=
                (surface.getAverageMeanCurvature() - exact_curv) *
                (surface.getAverageMeanCurvature() - exact_curv);
            l2_counter += 1.0;
          }
        }
      }
    }
  }

  double max_mean_curv_error_global = 0.0;
  double l2_mean_curv_error_global = 0.0, l1_mean_curv_error_global = 1.0,
         l2_counter_global = 0.0;
  double total_surface_global = 0.0, avg_mean_curv_global = 0.0;

  MPI_Allreduce(&l2_counter, &l2_counter_global, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&l2_mean_curv_error, &l2_mean_curv_error_global, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&l1_mean_curv_error, &l1_mean_curv_error_global, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&avg_mean_curv, &avg_mean_curv_global, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&total_surface, &total_surface_global, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&max_mean_curv_error, &max_mean_curv_error_global, 1,
                MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  l1_mean_curv_error_global /= l2_counter_global;
  l2_mean_curv_error_global /= l2_counter_global;
  l2_mean_curv_error_global = std::sqrt(l2_mean_curv_error_global);
  avg_mean_curv_global /= total_surface_global;
  total_surface_global = std::sqrt(total_surface_global);

  if (rank == 0 && print) {
    std::cout << "Linf  K = " << std::scientific << std::setprecision(2)
              << max_mean_curv_error / exact_curv << std::endl;
    std::cout << "L1    K = " << l1_mean_curv_error_global / exact_curv
              << std::endl;
    std::cout << "L2    K = " << l2_mean_curv_error_global / exact_curv
              << std::endl;
  }
  a_output->writeParametrizedInterface(a_time, surfaces);
}

void printError(const BasicMesh& mesh,
                const Data<IRL::GeneralMoments3D<2>>& liquid_moments,
                const Data<IRL::GeneralMoments3D<2>>& starting_liquid_moments) {
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    double linf_error_m0 = 0.0;
    double linf_error_m1 = 0.0;
    double linf_error_m2 = 0.0;
    double l1_error_m0 = 0.0;
    double l1_error_m1 = 0.0;
    double l1_error_m2 = 0.0;
    double l2_error_m0 = 0.0;
    double l2_error_m1 = 0.0;
    double l2_error_m2 = 0.0;
    double scale_m0 = 1.0 / std::pow(mesh.dx(), 3.0);
    double scale_m1 = 1.0 / std::pow(mesh.dx(), 4.0);
    double scale_m2 = 1.0 / std::pow(mesh.dx(), 5.0);
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
          const double liquid_volume_fraction =
              liquid_moments(i, j, k)[0] / mesh.cell_volume();
          if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
              liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
            auto mom_err =
                (liquid_moments(i, j, k) - starting_liquid_moments(i, j, k));
            linf_error_m0 = std::max(linf_error_m0, std::abs(mom_err[0]));
            linf_error_m1 = std::max(linf_error_m1, std::abs(mom_err[1]));
            linf_error_m1 = std::max(linf_error_m1, std::abs(mom_err[2]));
            linf_error_m1 = std::max(linf_error_m1, std::abs(mom_err[3]));
            linf_error_m2 = std::max(linf_error_m2, std::abs(mom_err[4]));
            linf_error_m2 = std::max(linf_error_m2, std::abs(2.0 * mom_err[5]));
            linf_error_m2 = std::max(linf_error_m2, std::abs(2.0 * mom_err[6]));
            linf_error_m2 = std::max(linf_error_m2, std::abs(mom_err[7]));
            linf_error_m2 = std::max(linf_error_m2, std::abs(2.0 * mom_err[8]));
            linf_error_m2 = std::max(linf_error_m2, std::abs(mom_err[9]));
            l1_error_m0 += std::abs(mom_err[0]);
            l1_error_m1 += std::abs(mom_err[1]);
            l1_error_m1 += std::abs(mom_err[2]);
            l1_error_m1 += std::abs(mom_err[3]);
            l1_error_m2 += std::abs(mom_err[4]);
            l1_error_m2 += std::abs(2.0 * mom_err[5]);
            l1_error_m2 += std::abs(2.0 * mom_err[6]);
            l1_error_m2 += std::abs(mom_err[7]);
            l1_error_m2 += std::abs(2.0 * mom_err[8]);
            l1_error_m2 += std::abs(mom_err[9]);
            l2_error_m0 += mom_err[0] * mom_err[0];
            l2_error_m1 += mom_err[1] * mom_err[1];
            l2_error_m1 += mom_err[2] * mom_err[2];
            l2_error_m1 += mom_err[3] * mom_err[3];
            l2_error_m2 += mom_err[4] * mom_err[4];
            l2_error_m2 += 4.0 * mom_err[5] * mom_err[5];
            l2_error_m2 += 4.0 * mom_err[6] * mom_err[6];
            l2_error_m2 += mom_err[7] * mom_err[7];
            l2_error_m2 += 4.0 * mom_err[8] * mom_err[8];
            l2_error_m2 += mom_err[9] * mom_err[9];
          }
        }
      }
    }
    l1_error_m0 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    l1_error_m1 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    l1_error_m2 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    l2_error_m0 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    l2_error_m1 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    l2_error_m2 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    linf_error_m0 *= scale_m0;
    linf_error_m1 *= scale_m1;
    linf_error_m2 *= scale_m2;
    l1_error_m0 *= scale_m0;
    l1_error_m1 *= scale_m1;
    l1_error_m2 *= scale_m2;
    l2_error_m0 = std::sqrt(l2_error_m0) * scale_m0;
    l2_error_m1 = std::sqrt(l2_error_m1) * scale_m1;
    l2_error_m2 = std::sqrt(l2_error_m2) * scale_m2;
    std::cout << std::scientific << std::setprecision(2)
              << "Linf M0 = " << linf_error_m0 << std::endl;
    std::cout << "Linf M1 = " << linf_error_m1 << std::endl;
    std::cout << "Linf M2 = " << linf_error_m2 << std::endl;
    std::cout << "L1   M0 = " << l1_error_m0 << std::endl;
    std::cout << "L1   M1 = " << l1_error_m1 << std::endl;
    std::cout << "L1   M2 = " << l1_error_m2 << std::endl;
    std::cout << "L2   M0 = " << l2_error_m0 << std::endl;
    std::cout << "L2   M1 = " << l2_error_m1 << std::endl;
    std::cout << "L2   M2 = " << l2_error_m2 << std::endl;
  }
}
