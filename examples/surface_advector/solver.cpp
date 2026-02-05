// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/surface_advector/solver.h"
#include "irl/geometry/polygons/polygon.h"

// Convert and store the mesh cells into localizers.
void initializeLocalizers(Data<IRL::PlanarLocalizer>* a_localizers) {
  // For each cell, construct cell as a RectangularCuboid and obtain localizer.
  const BasicMesh& mesh = a_localizers->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt lower_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        IRL::Pt upper_pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
        (*a_localizers)(i, j, k) =
            IRL::RectangularCuboid::fromBoundingPts(lower_pt, upper_pt)
                .getLocalizer();
      }
    }
  }
}

void initializeLocalizedInterfaces(
    const Data<IRL::PlanarLocalizer>& a_cell_localizers,
    const Data<IRL::SeparatorVariant>& a_interface,
    Data<IRL::LocalizedSeparatorVariantLink>* a_linked_localized_interfaces) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_linked_localized_interfaces)(i, j, k) =
            IRL::LocalizedSeparatorVariantLink(&a_cell_localizers(i, j, k),
                                               &a_interface(i, j, k));
      }
    }
  }
}

void setPhaseQuantities(const Data<IRL::SeparatorVariant>& a_interface,
                        Data<IRL::VolumeMoments>* a_liq_moments,
                        Data<IRL::VolumeMoments>* a_gas_moments) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        const auto moments =
            IRL::getVolumeMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
                cell, a_interface(i, j, k));
        (*a_liq_moments)(i, j, k) = moments[0];
        (*a_gas_moments)(i, j, k) = moments[1];
      }
    }
  }
  a_liq_moments->updateBorder();
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
                         const Data<IRL::VolumeMoments>& a_liq_moments,
                         const Data<IRL::SeparatorVariant>& a_interface,
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
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        liquid_volume_fraction_sum += liquid_volume_fraction;
        liquid_volume_sum += a_liq_moments(i, j, k).volume();
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
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::SeparatorVariant>& a_liquid_gas_interface,
    const Data<double>& a_surfactant_concentration, const double a_time,
    VTKOutput* a_output, const bool print) {
  using VolumeAndParaboloid =
      IRL::AddSurfaceOutput<IRL::Volume,
                            IRL::ParaboloidParametrizedSurfaceOutput>;
  using VolumeAndCylinder =
      IRL::AddSurfaceOutput<IRL::Volume,
                            IRL::CylinderParametrizedSurfaceOutput>;

  const BasicMesh& mesh = a_liq_moments.getMesh();

  std::vector<IRL::Polygon> polygons;
  std::vector<IRL::ParaboloidParametrizedSurfaceOutput> paraboloids;
  std::vector<IRL::CylinderParametrizedSurfaceOutput> cylinders;
  std::vector<double> polygon_data, paraboloid_data, cylinder_data;

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const double liquid_volume_fraction =
            a_liq_moments(i, j, k).volume() / mesh.cell_volume();
        if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
            liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          if (const auto ptr = std::get_if<IRL::PlanarSeparator>(
                  &a_liquid_gas_interface(i, j, k))) {
            const auto polygon =
                IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(cell, *ptr,
                                                                     (*ptr)[0]);
            polygons.push_back(polygon);
            polygon_data.push_back(a_surfactant_concentration(i, j, k));
          } else if (const auto ptr = std::get_if<IRL::Paraboloid>(
                         &a_liquid_gas_interface(i, j, k))) {
            auto volume_and_surface =
                IRL::getVolumeMoments<VolumeAndParaboloid>(cell, *ptr);
            auto surface = volume_and_surface.getSurface();
            double length_scale = std::min(0.25 * mesh.dx(), 1.0e-2);
            surface.setLengthScale(length_scale);
            if (surface.getSurfaceArea() >
                1.0e-6 * length_scale * length_scale) {
              paraboloids.push_back(surface);
              paraboloid_data.push_back(a_surfactant_concentration(i, j, k));
            }
          } else if (const auto ptr = std::get_if<IRL::Cylinder>(
                         &a_liquid_gas_interface(i, j, k))) {
            auto volume_and_surface =
                IRL::getVolumeMoments<VolumeAndCylinder>(cell, *ptr);
            auto surface = volume_and_surface.getSurface();
            double length_scale = std::min(0.25 * mesh.dx(), 1.0e-2);
            surface.setLengthScale(length_scale);
            if (surface.getSurfaceArea() >
                1.0e-6 * length_scale * length_scale) {
              cylinders.push_back(surface);
              cylinder_data.push_back(a_surfactant_concentration(i, j, k));
            }
          }
        }
      }
    }
  }

  a_output->writeVTKInterface(a_time, polygons, paraboloids, cylinders,
                              polygon_data, paraboloid_data, cylinder_data);
}

void printError(const BasicMesh& mesh,
                const Data<IRL::VolumeMoments>& liq_moments,
                const Data<IRL::VolumeMoments>& starting_liq_moments) {
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    double linf_error_m0 = 0.0, linf_error_m1 = 0.0;
    double l1_error_m0 = 0.0, l1_error_m1 = 0.0;
    double l2_error_m0 = 0.0, l2_error_m1 = 0.0;
    double scale_m0 = 1.0 / std::pow(mesh.dx(), 3.0);
    double scale_m1 = 1.0 / std::pow(mesh.dx(), 4.0);
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
          const double liquid_volume_fraction =
              liq_moments(i, j, k).volume() / mesh.cell_volume();
          if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
              liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
            auto mom_err =
                (liq_moments(i, j, k) - starting_liq_moments(i, j, k));
            linf_error_m0 = std::max(linf_error_m0, std::abs(mom_err.volume()));
            linf_error_m1 =
                std::max(linf_error_m1, std::abs(mom_err.centroid()[0]));
            linf_error_m1 =
                std::max(linf_error_m1, std::abs(mom_err.centroid()[1]));
            linf_error_m1 =
                std::max(linf_error_m1, std::abs(mom_err.centroid()[2]));
            l1_error_m0 += std::abs(mom_err.volume());
            l1_error_m1 += std::abs(mom_err.centroid()[0]);
            l1_error_m1 += std::abs(mom_err.centroid()[1]);
            l1_error_m1 += std::abs(mom_err.centroid()[2]);
            l2_error_m0 += mom_err.volume() * mom_err.volume();
            l2_error_m1 += mom_err.centroid()[0] * mom_err.centroid()[0];
            l2_error_m1 += mom_err.centroid()[1] * mom_err.centroid()[1];
            l2_error_m1 += mom_err.centroid()[2] * mom_err.centroid()[2];
          }
        }
      }
    }
    l1_error_m0 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    l1_error_m1 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    l2_error_m0 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    l2_error_m1 /=
        (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
    linf_error_m0 *= scale_m0;
    linf_error_m1 *= scale_m1;
    l1_error_m0 *= scale_m0;
    l1_error_m1 *= scale_m1;
    l2_error_m0 = std::sqrt(l2_error_m0) * scale_m0;
    l2_error_m1 = std::sqrt(l2_error_m1) * scale_m1;
    std::cout << std::scientific << std::setprecision(2)
              << "Linf M0 = " << linf_error_m0 << std::endl;
    std::cout << "Linf M1 = " << linf_error_m1 << std::endl;
    std::cout << "L1   M0 = " << l1_error_m0 << std::endl;
    std::cout << "L1   M1 = " << l1_error_m1 << std::endl;
    std::cout << "L2   M0 = " << l2_error_m0 << std::endl;
    std::cout << "L2   M1 = " << l2_error_m1 << std::endl;
  }
}

void connectMesh(
    const BasicMesh& a_mesh,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface) {
  IRL::LocalizedSeparatorVariantLink* neighbor_ptr;
  // Provide mesh connectivity information.
  int unique_id = 0;

  for (int i = a_mesh.imino(); i <= a_mesh.imaxo(); ++i) {
    for (int j = a_mesh.jmino(); j <= a_mesh.jmaxo(); ++j) {
      for (int k = a_mesh.kmino(); k <= a_mesh.kmaxo(); ++k) {
        (*a_link_localized_interface)(i, j, k).setId(unique_id);
        neighbor_ptr = i - 1 < a_mesh.imino()
                           ? nullptr
                           : &(*a_link_localized_interface)(i - 1, j, k);
        (*a_link_localized_interface)(i, j, k).setEdgeConnectivity(
            0, neighbor_ptr);
        neighbor_ptr = i + 1 > a_mesh.imaxo()
                           ? nullptr
                           : &(*a_link_localized_interface)(i + 1, j, k);
        (*a_link_localized_interface)(i, j, k).setEdgeConnectivity(
            1, neighbor_ptr);
        neighbor_ptr = j - 1 < a_mesh.jmino()
                           ? nullptr
                           : &(*a_link_localized_interface)(i, j - 1, k);
        (*a_link_localized_interface)(i, j, k).setEdgeConnectivity(
            2, neighbor_ptr);
        neighbor_ptr = j + 1 > a_mesh.jmaxo()
                           ? nullptr
                           : &(*a_link_localized_interface)(i, j + 1, k);
        (*a_link_localized_interface)(i, j, k).setEdgeConnectivity(
            3, neighbor_ptr);
        neighbor_ptr = k - 1 < a_mesh.kmino()
                           ? nullptr
                           : &(*a_link_localized_interface)(i, j, k - 1);
        (*a_link_localized_interface)(i, j, k).setEdgeConnectivity(
            4, neighbor_ptr);
        neighbor_ptr = k + 1 > a_mesh.kmaxo()
                           ? nullptr
                           : &(*a_link_localized_interface)(i, j, k + 1);
        (*a_link_localized_interface)(i, j, k).setEdgeConnectivity(
            5, neighbor_ptr);
        ++unique_id;
      }
    }
  }
}
