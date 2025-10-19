// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/implicit_surface_reconstruction/solver.h"
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

void reconstructSurface(const std::string& a_reconstruction_method,
                        const int a_nx) {
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Set mesh
  const int GC = 5;
  // const IRL::Pt lower_domain(-0.5, -0.5, -0.5);
  // const IRL::Pt upper_domain(0.5, 0.5, 0.5);
  const IRL::Pt lower_domain(-2.5, -2.5, -2.5);
  const IRL::Pt upper_domain(2.5, 2.5, 2.5);
  BasicMesh cc_mesh(a_nx, a_nx, a_nx, GC);
  cc_mesh.setCellBoundaries(lower_domain, upper_domain);

  // Allocate local data
  Data<double> velU(&cc_mesh), velV(&cc_mesh), velW(&cc_mesh);
  Data<double> vf(&cc_mesh);
  Data<IRL::VolumeMoments> liq_moments(&cc_mesh), gas_moments(&cc_mesh);

  // Allocate interfaces and localizers
  // Data<IRL::PlanarLocalizer> cell_localizers(&cc_mesh);
  Data<IRL::SeparatorVariant> interface(&cc_mesh), interface_2(&cc_mesh);
  // Data<IRL::LocalizedSeparatorVariantLink>
  // link_localized_interface(&cc_mesh);

  // // Generate localizers from cells
  // initializeLocalizers(&cell_localizers);
  // // Link interfaces to localizers
  // initializeLocalizedInterfaces(cell_localizers, interface,
  //                               &link_localized_interface);
  // // Generate mesh connectivity
  // connectMesh(cc_mesh, &link_localized_interface);

  // Set constants in IRL
  IRL::setMinimumVolumeToTrack(10.0 * DBL_EPSILON * cc_mesh.dx() *
                               cc_mesh.dy() * cc_mesh.dz());
  IRL::setVolumeFractionBounds(1.0e-13);
  IRL::setVolumeFractionTolerance(1.0e-12);

  if (rank == 0) {
    std::cout << "Selected volume fraction bounds = "
              << IRL::global_constants::VF_LOW << std::endl;
  }

  // initializing implicit surface
  constexpr std::size_t max_refine_level = 4;
  // IRL::Sphere<double, max_refine_level> sphere(0.0, 0.0, 0.0, 0.15);
  IRL::Orthocircle<double, max_refine_level> orthocircle;

  // initiailizing volume fraction
  for (int i = cc_mesh.imin(); i <= cc_mesh.imax(); i++) {
    for (int j = cc_mesh.jmin(); j <= cc_mesh.jmax(); j++) {
      for (int k = cc_mesh.kmin(); k <= cc_mesh.kmax(); k++) {
        IRL::Pt x0(cc_mesh.x(i), cc_mesh.y(j), cc_mesh.z(k));
        IRL::Pt x1(cc_mesh.x(i + 1), cc_mesh.y(j + 1), cc_mesh.z(k + 1));
        IRL::RectangularCuboid cell =
            IRL::RectangularCuboid::fromBoundingPts(x0, x1);
        // IRL::ImplicitSurfaceCutter<IRL::Sphere<double, max_refine_level>,
        //                            IRL::VolumeMoments>
        //     cutter(sphere, cell);
        IRL::ImplicitSurfaceCutter<IRL::Orthocircle<double, max_refine_level>,
                                   IRL::VolumeMoments>
            cutter(orthocircle, cell);
        liq_moments(i, j, k) = cutter.computeVolumeMoments();
        gas_moments(i, j, k) = cell.calculateMoments() - liq_moments(i, j, k);
      }
    }
  }

  VTKOutput vtk_io("viz_out", "viz", cc_mesh);
  VTKOutput vtk_io2("viz_out2", "viz2", cc_mesh);
  double simulation_time = 0.0;

  if (rank == 0) {
    vtk_io.writeVTKFile(simulation_time);
  }

  // construct interface
  std::string method_1 = "Jibben";
  std::string method_2 = "PU";
  // getReconstruction(a_reconstruction_method, liq_moments, gas_moments, 0.0,
  //                   velU, velV, velW, &interface);

  getReconstruction(method_1, liq_moments, gas_moments, 0.0, velU, velV, velW,
                    &interface);
  getReconstruction(method_2, liq_moments, gas_moments, 0.0, velU, velV, velW,
                    &interface_2);

  writeInterfaceToFile(liq_moments, interface, simulation_time, &vtk_io, true);
  writeInterfaceToFile(liq_moments, interface_2, simulation_time, &vtk_io2,
                       true);
}
