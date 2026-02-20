// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/variant_advector/solver.h"
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
  printf("%10s %20s %12s %20s %20s %20s %20s %20s %20s %20s %20s\n",
         "Iteration", "Time", "CFL", "liquidVFSum", "liquidVolSum",
         "ChangeLiquidVFSum", "ChangeLiquidVolSum", "AdvectionDuration",
         "ReconDuration", "OutputDuration", "InterfaceCells");
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
    const double a_time, VTKOutput* a_output, const bool print) {
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
            }
          }
        }
      }
    }
  }

  a_output->writeVTKInterface(a_time, polygons, paraboloids, cylinders);
}

void writeInterfaceWithScalarToFile(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::SeparatorVariant>& a_liquid_gas_interface,
    const std::vector<InterfaceScalarField>& a_scalar_fields,
    const double a_time, VTKOutput* a_output, const bool print) {
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
            }
          }
        }
      }
    }
  }

  a_output->writeVTKInterfaceWithScalar(a_time, polygons, paraboloids,
                                        cylinders, a_scalar_fields);
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
          // if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          //     liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
          auto mom_err = (liq_moments(i, j, k) - starting_liq_moments(i, j, k));
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
          // }
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

// calculate moments of the polygon
std::vector<double> calculatePolygonSurfaceMoments(
    const IRL::Polygon& polygon) {
  std::vector<double> moments(10, 0.0);

  const int n = polygon.getNumberOfVertices();

  if (n < 3) return moments;

  // common point for all triangles
  const IRL::Pt& v0 = polygon[0];

  // finding moments for each triangle
  for (int i = 1; i < n - 1; i++) {
    const IRL::Pt& v1 = polygon[i];
    const IRL::Pt& v2 = polygon[i + 1];

    // forming triangle
    IRL::Polygon triangle;
    triangle.addVertex(v0);
    triangle.addVertex(v1);
    triangle.addVertex(v2);
    triangle.calculateAndSetPlaneOfExistence();

    // M0
    double triangle_area = triangle.calculateVolume();
    moments[0] += triangle_area;

    // M1
    IRL::Pt triangle_centroid = (v0 + v1 + v2) * (1.0 / 3.0);
    moments[1] += triangle_centroid[0] * triangle_area;
    moments[2] += triangle_centroid[1] * triangle_area;
    moments[3] += triangle_centroid[2] * triangle_area;

    // M2
    const auto a = v1 - v0;
    const auto b = v2 - v0;
    const double factor = 2.0 * triangle_area;

    auto accumulate = [&](const IRL::Pt& u, const IRL::Pt& v, double scale) {
      moments[4] += scale * u[0] * v[0];  // M00
      moments[5] += scale * u[0] * v[1];  // M01
      moments[6] += scale * u[0] * v[2];  // M02
      moments[7] += scale * u[1] * v[1];  // M11
      moments[8] += scale * u[1] * v[2];  // M12
      moments[9] += scale * u[2] * v[2];  // M22
    };
    accumulate(v0, v0, factor * 0.5);
    accumulate(v0, a, factor * (1.0 / 6.0));
    accumulate(a, v0, factor * (1.0 / 6.0));
    accumulate(v0, b, factor * (1.0 / 6.0));
    accumulate(b, v0, factor * (1.0 / 6.0));
    accumulate(a, a, factor * (1.0 / 12.0));
    accumulate(b, b, factor * (1.0 / 12.0));
    accumulate(a, b, factor * (1.0 / 24.0));
    accumulate(b, a, factor * (1.0 / 24.0));
  }

  return moments;
}

// recentering moments
std::vector<double> recenterMoments(const std::vector<double>& moments,
                                    const IRL::Pt& xc) {
  std::vector<double> centered_moments(10, 0.0);

  Eigen::Vector3d x_ref(xc[0], xc[1], xc[2]);

  double M0 = moments[0];
  Eigen::Vector3d M1(moments[1], moments[2], moments[3]);
  Eigen::Matrix3d M2 = Eigen::Matrix3d::Zero();
  M2(0, 0) = moments[4];
  M2(0, 1) = moments[5];
  M2(0, 2) = moments[6];
  M2(1, 0) = M2(0, 1);
  M2(1, 1) = moments[7];
  M2(1, 2) = moments[8];
  M2(2, 0) = M2(0, 2);
  M2(2, 1) = M2(1, 2);
  M2(2, 2) = moments[9];

  Eigen::Vector3d M1_recentered = M1 - (M0 * x_ref);
  Eigen::Matrix3d M2_recentered = M2 - (x_ref * M1.transpose()) -
                                  (M1 * x_ref.transpose()) +
                                  (M0 * x_ref * x_ref.transpose());

  centered_moments[0] = M0;
  centered_moments[1] = M1_recentered[0];
  centered_moments[2] = M1_recentered[1];
  centered_moments[3] = M1_recentered[2];
  centered_moments[4] = M2_recentered(0, 0);
  centered_moments[5] = M2_recentered(0, 1);
  centered_moments[6] = M2_recentered(0, 2);
  centered_moments[7] = M2_recentered(1, 1);
  centered_moments[8] = M2_recentered(1, 2);
  centered_moments[9] = M2_recentered(2, 2);

  return centered_moments;
}

// getting starting and ending moments
void writeMomentsToBinary(const Data<IRL::VolumeMoments>& liq_moments_1,
                          const Data<IRL::VolumeMoments>& liq_moments_2,
                          const Data<IRL::SeparatorVariant>& interfaces_1,
                          const Data<IRL::SeparatorVariant>& interfaces_2,
                          const std::string& output_dir,
                          const std::string& case_name,
                          const std::string& reconstruction_method) {
  const BasicMesh& mesh = interfaces_1.getMesh();
  const int Nx = mesh.getNx();
  const int Ny = mesh.getNy();
  const int Nz = mesh.getNz();

  const int mom_count = 10;

  // storing moments in vectors
  Data<std::vector<double>> volume_moments_1_centered(&mesh);
  Data<std::vector<double>> volume_moments_2_centered(&mesh);
  Data<std::vector<double>> surface_moments_1_centered(&mesh);
  Data<std::vector<double>> surface_moments_2_centered(&mesh);

  // default everything to zero
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        volume_moments_1_centered(i, j, k).resize(mom_count, 0.0);
        volume_moments_2_centered(i, j, k).resize(mom_count, 0.0);
        surface_moments_1_centered(i, j, k).resize(mom_count, 0.0);
        surface_moments_2_centered(i, j, k).resize(mom_count, 0.0);
      }
    }
  }

  using VolumeMomentsAndSurface =
      IRL::AddSurfaceOutput<IRL::VolumeMoments,
                            IRL::ParaboloidParametrizedSurfaceOutput>;

  // --- starting ---
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        const double vf_1 =
            liq_moments_1(i, j, k).volume() / mesh.cell_volume();
        if (vf_1 > IRL::global_constants::VF_LOW &&
            vf_1 < IRL::global_constants::VF_HIGH) {
          // mixed
          IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
          IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
          IRL::Pt xc = 0.5 * (x0 + x1);
          IRL::RectangularCuboid cell =
              IRL::RectangularCuboid::fromBoundingPts(x0, x1);

          // volume moments
          IRL::GeneralMoments3D<2> vm1 = std::visit(
              [&](const auto& interface) -> IRL::GeneralMoments3D<2> {
                return IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(
                    cell, interface);
              },
              interfaces_1(i, j, k));

          // storing volume moments
          std::vector<double> temp_vm1(10, 0.0);
          for (int m = 0; m < mom_count; m++) {
            temp_vm1[m] = vm1[m];
          }
          volume_moments_1_centered(i, j, k) = recenterMoments(temp_vm1, xc);

          // surface moments (starting is always a paraboloid)
          const IRL::Paraboloid paraboloid =
              std::get<IRL::Paraboloid>(interfaces_1(i, j, k));
          auto surface =
              IRL::getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid)
                  .getSurface();
          IRL::GeneralSurfaceMoments3D<2> sm1 = surface.getSurfaceMoments<2>();

          // storing surface moments
          std::vector<double> temp_sm1(10, 0.0);
          for (int m = 0; m < mom_count; m++) {
            temp_sm1[m] = sm1[m];
          }
          surface_moments_1_centered(i, j, k) = recenterMoments(temp_sm1, xc);
        }
      }
    }
  }

  // --- ending ---
  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        const double vf_2 =
            liq_moments_2(i, j, k).volume() / mesh.cell_volume();
        if (vf_2 > IRL::global_constants::VF_LOW &&
            vf_2 < IRL::global_constants::VF_HIGH) {
          // mixed
          IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
          IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
          IRL::Pt xc = 0.5 * (x0 + x1);
          IRL::RectangularCuboid cell =
              IRL::RectangularCuboid::fromBoundingPts(x0, x1);

          // volume moments
          IRL::GeneralMoments3D<2> vm2 = std::visit(
              [&](const auto& interface) -> IRL::GeneralMoments3D<2> {
                return IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(
                    cell, interface);
              },
              interfaces_2(i, j, k));

          // storing volume moments
          std::vector<double> temp_vm2(10, 0.0);
          for (int m = 0; m < mom_count; m++) {
            temp_vm2[m] = vm2[m];
          }
          volume_moments_2_centered(i, j, k) = recenterMoments(temp_vm2, xc);

          // calculating and storing surface moments
          std::vector<double> temp_sm2(10, 0.0);
          std::visit(
              [&](auto const& separator) {
                using T = std::decay_t<decltype(separator)>;

                if constexpr (std::is_same_v<T, IRL::Paraboloid>) {
                  const IRL::Paraboloid& paraboloid = separator;
                  auto surface = IRL::getVolumeMoments<VolumeMomentsAndSurface>(
                                     cell, paraboloid)
                                     .getSurface();
                  IRL::GeneralSurfaceMoments3D<2> sm2 =
                      surface.template getSurfaceMoments<2>();
                  for (int m = 0; m < mom_count; m++) {
                    temp_sm2[m] = sm2[m];
                  }
                } else if constexpr (std::is_same_v<T, IRL::PlanarSeparator>) {
                  const IRL::PlanarSeparator& planar_separator = separator;
                  IRL::Polygon polygon =
                      IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                          cell, planar_separator, planar_separator[0]);
                  if (polygon.getNumberOfVertices() > 2) {
                    polygon.calculateAndSetPlaneOfExistence();
                    temp_sm2 = calculatePolygonSurfaceMoments(polygon);
                  }
                }
              },
              interfaces_2(i, j, k));
          surface_moments_2_centered(i, j, k) = recenterMoments(temp_sm2, xc);
        }
      }
    }
  }

  // --- writing data to binary ---
  std::string filename = output_dir + "/moments_" + case_name + "_" +
                         reconstruction_method + "_Nx" + std::to_string(Nx) +
                         ".bin";

  std::ofstream out(filename, std::ios::binary);
  if (!out) {
    throw std::runtime_error("Cannot open " + filename + " for writing.");
  }

  // simple header: magic, version, Nx, Ny, Nz
  int32_t magic = 0x4D4F4D33;  // 'MOM3' in hex
  int32_t version = 1;
  int32_t dims[3] = {Nx, Ny, Nz};

  out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
  out.write(reinterpret_cast<const char*>(&version), sizeof(version));
  out.write(reinterpret_cast<const char*>(dims), sizeof(dims));

  double buf[40];

  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        const auto& m1v = volume_moments_1_centered(i, j, k);
        const auto& m2v = volume_moments_2_centered(i, j, k);
        const auto& m1s = surface_moments_1_centered(i, j, k);
        const auto& m2s = surface_moments_2_centered(i, j, k);

        // volume start
        buf[0] = m1v[0];
        buf[1] = m1v[1];
        buf[2] = m1v[2];
        buf[3] = m1v[3];
        buf[4] = m1v[4];
        buf[5] = m1v[5];
        buf[6] = m1v[6];
        buf[7] = m1v[7];
        buf[8] = m1v[8];
        buf[9] = m1v[9];

        // volume end
        buf[10] = m2v[0];
        buf[11] = m2v[1];
        buf[12] = m2v[2];
        buf[13] = m2v[3];
        buf[14] = m2v[4];
        buf[15] = m2v[5];
        buf[16] = m2v[6];
        buf[17] = m2v[7];
        buf[18] = m2v[8];
        buf[19] = m2v[9];

        // surface start
        buf[20] = m1s[0];
        buf[21] = m1s[1];
        buf[22] = m1s[2];
        buf[23] = m1s[3];
        buf[24] = m1s[4];
        buf[25] = m1s[5];
        buf[26] = m1s[6];
        buf[27] = m1s[7];
        buf[28] = m1s[8];
        buf[29] = m1s[9];

        // surface end
        buf[30] = m2s[0];
        buf[31] = m2s[1];
        buf[32] = m2s[2];
        buf[33] = m2s[3];
        buf[34] = m2s[4];
        buf[35] = m2s[5];
        buf[36] = m2s[6];
        buf[37] = m2s[7];
        buf[38] = m2s[8];
        buf[39] = m2s[9];

        out.write(reinterpret_cast<const char*>(buf), sizeof(buf));
      }
    }
  }

  out.close();
}