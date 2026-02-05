// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_SURFACE_ADVECTOR_SOLVER_H_
#define EXAMPLES_SURFACE_ADVECTOR_SOLVER_H_

#include <mpi.h>
#include <stdio.h>
#include <sys/stat.h>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>

#include "irl/cylinder_reconstruction/cylinder_parametrized_surface.h"
#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid_parametrized_surface.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include "examples/surface_advector/basic_mesh.h"
#include "examples/surface_advector/data.h"
#include "examples/surface_advector/reconstruction_types.h"
#include "examples/surface_advector/vof_advection.h"
#include "examples/surface_advector/vtk.h"

#include "examples/surface_advector/block_i.h"
#include "examples/surface_advector/translation_3d.h"

/// \brief Handles running and advancing the solution according to provided
/// static functions in structs.
template <class SimulationType>
int runSimulation(const std::string& a_case_name,
                  const std::string& a_advection_method,
                  const std::string& a_reconstruction_method,
                  const double a_max_cfl, const double a_end_time,
                  const int a_viz_frequency, const int a_nx);

// \brief Convert and store the mesh cells into localizers.
void initializeLocalizers(Data<IRL::PlanarLocalizer>* a_localizers);

/// \brief Initialize the linked localized interfaces used during advection.
void initializeLocalizedInterfaces(
    const Data<IRL::PlanarLocalizer>& a_cell_localizers,
    const Data<IRL::SeparatorVariant>& a_interface,
    Data<IRL::LocalizedSeparatorVariantLink>* a_linked_localized_interface);

/// \brief Set phase quantities according to the given
void setPhaseQuantities(const Data<IRL::SeparatorVariant>& a_interface,
                        Data<IRL::VolumeMoments>* a_liq_moments,
                        Data<IRL::VolumeMoments>* a_gas_moments);

/// \brief Write out the header for the diagnostics.
void writeDiagnosticsHeader(void);

/// \brief Write out diagnostics to the screen.
void writeOutDiagnostics(const int a_iteration, const double a_timestep,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<double>& a_W,
                         const Data<IRL::VolumeMoments>& a_liq_moments,
                         const Data<IRL::SeparatorVariant>& a_interface,
                         std::chrono::duration<double> a_VOF_duration,
                         std::chrono::duration<double> a_recon_duration,
                         std::chrono::duration<double> a_write_duration);

/// \brief Generates triangulated surface and writes to provided VTK file
void writeInterfaceToFile(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::SeparatorVariant>& a_liquid_gas_interface,
    const Data<double>& a_surfactant_concentration, const double a_time,
    VTKOutput* a_output, const bool print);

void printError(const BasicMesh& mesh,
                const Data<IRL::VolumeMoments>& liq_moments,
                const Data<IRL::VolumeMoments>& starting_liq_moments);

void connectMesh(
    const BasicMesh& a_mesh,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface);

//******************************************************************* //
//     Template function definitions placed below this.
//******************************************************************* //
template <class SimulationType>
int runSimulation(const std::string& a_case_name,
                  const std::string& a_advection_method,
                  const std::string& a_reconstruction_method,
                  const double a_max_cfl, const double a_end_time,
                  const int a_viz_frequency, const int a_nx) {
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Set mesh
  BasicMesh cc_mesh = SimulationType::setMesh(a_nx);
  const double timestep = SimulationType::getTimeStep(cc_mesh, a_max_cfl);

  // Allocate local data
  Data<double> velU(&cc_mesh), velV(&cc_mesh), velW(&cc_mesh);
  Data<IRL::VolumeMoments> liq_moments(&cc_mesh), gas_moments(&cc_mesh);
  Data<double> surfactant_mass(&cc_mesh), surfactant_concentration(&cc_mesh);

  // Allocate interfaces and localizers
  Data<IRL::PlanarLocalizer> cell_localizers(&cc_mesh);
  Data<IRL::SeparatorVariant> interface(&cc_mesh);
  Data<IRL::LocalizedSeparatorVariantLink> link_localized_interface(&cc_mesh);

  // Generate localizers from cells
  initializeLocalizers(&cell_localizers);
  // Link interfaces to localizers
  initializeLocalizedInterfaces(cell_localizers, interface,
                                &link_localized_interface);
  // Generate mesh connectivity
  connectMesh(cc_mesh, &link_localized_interface);

  // Set constants in IRL
  IRL::setMinimumVolumeToTrack(10.0 * DBL_EPSILON * cc_mesh.dx() *
                               cc_mesh.dy() * cc_mesh.dz());
  IRL::setVolumeFractionBounds(1.0e-13);
  IRL::setVolumeFractionTolerance(1.0e-12);

  if (rank == 0) {
    std::cout << "Selected volume fraction bounds = "
              << IRL::global_constants::VF_LOW << std::endl;
  }

  // Initialize case data
  SimulationType::initialize(&velU, &velV, &velW, &interface, 0.0, a_end_time);
  setPhaseQuantities(interface, &liq_moments, &gas_moments);
  const auto starting_liq_moments = liq_moments;

  // Initialize block I polyhedron
  IRL::BlockI blockI = IRL::BlockI({IRL::Pt(-0.5 * 0.69, -0.5, 0.75),
                                    IRL::Pt(0.5 * 0.69, -0.5, 0.75),
                                    IRL::Pt(0.5 * 0.69, -0.23, 0.75),
                                    IRL::Pt(0.2753623188 * 0.69, -0.23, 0.75),
                                    IRL::Pt(0.2753623188 * 0.69, 0.23, 0.75),
                                    IRL::Pt(0.5 * 0.69, 0.23, 0.75),
                                    IRL::Pt(0.5 * 0.69, 0.5, 0.75),
                                    IRL::Pt(-0.5 * 0.69, 0.5, 0.75),
                                    IRL::Pt(-0.5 * 0.69, 0.23, 0.75),
                                    IRL::Pt(-0.2753623188 * 0.69, 0.23, 0.75),
                                    IRL::Pt(-0.2753623188 * 0.69, -0.23, 0.75),
                                    IRL::Pt(-0.5 * 0.69, -0.23, 0.75),
                                    IRL::Pt(-0.5 * 0.69, -0.5, -0.75),
                                    IRL::Pt(0.5 * 0.69, -0.5, -0.75),
                                    IRL::Pt(0.5 * 0.69, -0.23, -0.75),
                                    IRL::Pt(0.2753623188 * 0.69, -0.23, -0.75),
                                    IRL::Pt(0.2753623188 * 0.69, 0.23, -0.75),
                                    IRL::Pt(0.5 * 0.69, 0.23, -0.75),
                                    IRL::Pt(0.5 * 0.69, 0.5, -0.75),
                                    IRL::Pt(-0.5 * 0.69, 0.5, -0.75),
                                    IRL::Pt(-0.5 * 0.69, 0.23, -0.75),
                                    IRL::Pt(-0.2753623188 * 0.69, 0.23, -0.75),
                                    IRL::Pt(-0.2753623188 * 0.69, -0.23, -0.75),
                                    IRL::Pt(-0.5 * 0.69, -0.23, -0.75)});

  const double normalization_factor = 0.25;
  for (auto& vertex : blockI) {
    for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
      vertex[d] *= normalization_factor;
    }
  }
  const auto centroid = IRL::Pt(0.5, 0.5, 0.75);
  for (auto& vertex : blockI) {
    for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
      vertex[d] += centroid[d];
    }
  }

  // // Print block I polyhedron
  // IRL::HalfEdgePolyhedron<IRL::Pt> half_edge;
  // blockI.setHalfEdgeVersion(&half_edge);
  // auto seg_half_edge = half_edge.generateSegmentedPolyhedron();
  // std::ofstream myfile("blockI.vtu");
  // if (myfile.is_open()) {
  //   myfile << seg_half_edge;
  //   myfile.close();
  // }

  VTKOutput vtk_io("viz_out", "viz", cc_mesh);
  vtk_io.addData("VelocityX", velU);
  vtk_io.addData("VelocityY", velV);
  vtk_io.addData("VelocityZ", velW);
  double simulation_time = 0.0;
  int iteration = 0;

  if (rank == 0) {
    vtk_io.writeVTKFile(simulation_time);
  }
  getReconstruction(a_reconstruction_method, liq_moments, gas_moments, 0.0,
                    velU, velV, velW, &interface);
  resetMoments(link_localized_interface, &liq_moments, &gas_moments);

  // Initialize surfactant concentration
  for (int i = cc_mesh.imin(); i <= cc_mesh.imax(); ++i) {
    for (int j = cc_mesh.jmin(); j <= cc_mesh.jmax(); ++j) {
      for (int k = cc_mesh.kmin(); k <= cc_mesh.kmax(); ++k) {
        surfactant_concentration(i, j, k) = 0.0;
        surfactant_mass(i, j, k) = 0.0;
        const double vfrac =
            liq_moments(i, j, k).volume() / cc_mesh.cell_volume();
        if (vfrac >= IRL::global_constants::VF_LOW &&
            vfrac <= IRL::global_constants::VF_HIGH) {
          const IRL::Pt lower_pt(cc_mesh.x(i), cc_mesh.y(j), cc_mesh.z(k));
          const IRL::Pt upper_pt(cc_mesh.x(i + 1), cc_mesh.y(j + 1),
                                 cc_mesh.z(k + 1));
          const auto cuboid_cell =
              IRL::RectangularCuboid::fromBoundingPts(lower_pt, upper_pt);

          IRL::PlanarLocalizer localizer = cuboid_cell.getLocalizer();
          IRL::HalfEdgePolyhedron<IRL::Pt> half_edge;
          blockI.setHalfEdgeVersion(&half_edge);
          auto seg_half_edge = half_edge.generateSegmentedPolyhedron();
          decltype(seg_half_edge) dummy_clipped_polytope;
          for (int n = 0; n < 6; n++) {
            splitHalfEdgePolytope(&seg_half_edge, &dummy_clipped_polytope,
                                  &half_edge, localizer[n]);
          }

          if (const IRL::PlanarSeparator* separator =
                  std::get_if<IRL::PlanarSeparator>(
                      &(interface(i, j, k)))) {  // If plane
                                                 // const auto polygon =
            //     IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            //         cuboid_cell, *separator, (*separator)[0]);
            // auto polygon_clipped =
            //     IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            //         blockI, *separator, (*separator)[0]);
            // for (int n = 0; n < 6; n++) {
            //   polygon_clipped =
            //       cutPolygonByPlane(polygon_clipped, localizer[n], 1.0);
            // }
            // auto surface_area = polygon.calculateVolume();
            // auto surface_area_clipped = polygon_clipped.calculateVolume();
            // surfactant_concentration(i, j, k) =
            //     surface_area_clipped / surface_area;
            // surfactant_mass(i, j, k) = surface_area_clipped;
          } else if (const IRL::Paraboloid* paraboloid =
                         std::get_if<IRL::Paraboloid>(
                             &(interface(i, j, k)))) {  // If paraboloid
            using VolumeAndSurface =
                IRL::AddSurfaceOutput<IRL::Volume,
                                      IRL::ParaboloidParametrizedSurfaceOutput>;
            auto surface_and_volume = IRL::getVolumeMoments<VolumeAndSurface>(
                cuboid_cell, *paraboloid);
            auto surface = surface_and_volume.getSurface();
            auto surface_area = surface.getSurfaceArea();

            auto surface_and_volume_clipped =
                IRL::intersectPolyhedronWithParaboloid<VolumeAndSurface>(
                    &seg_half_edge, &half_edge, *paraboloid);
            auto surface_clipped = surface_and_volume_clipped.getSurface();
            auto surface_area_clipped = surface_clipped.getSurfaceArea();
            surfactant_concentration(i, j, k) =
                surface_area_clipped / surface_area;
            surfactant_mass(i, j, k) = surface_area_clipped;
          }
        }
      }
    }
  }

  if (rank == 0) {
    writeInterfaceToFile(liq_moments, interface, surfactant_concentration,
                         simulation_time, &vtk_io, true);
    writeDiagnosticsHeader();
  }

  std::string output_folder = "viz";
  const int dir_err = mkdir(output_folder.c_str(), 0777);
  std::chrono::duration<double> advect_VOF_time(0.0);
  std::chrono::duration<double> recon_time(0.0);
  std::chrono::duration<double> write_time(0.0);
  if (rank == 0) {
    writeOutDiagnostics(iteration, timestep, simulation_time, velU, velV, velW,
                        liq_moments, interface, advect_VOF_time, recon_time,
                        write_time);
    // printError(cc_mesh, liq_moments, starting_liq_moments);
  }
  while (simulation_time < a_end_time) {
    const double time_step_to_use =
        std::fmin(timestep, a_end_time - simulation_time);
    SimulationType::setVelocity(simulation_time + 0.5 * time_step_to_use, &velU,
                                &velV, &velW);

    auto start = std::chrono::system_clock::now();
    advectVOF(a_advection_method, a_reconstruction_method, time_step_to_use,
              simulation_time, velU, velV, velW, &interface,
              &link_localized_interface, &liq_moments, &gas_moments,
              &surfactant_mass);

    if (rank == 0) {
      writeOutDiagnostics(iteration + 1, time_step_to_use,
                          simulation_time + time_step_to_use, velU, velV, velW,
                          liq_moments, interface, advect_VOF_time, recon_time,
                          write_time);
      if (simulation_time + time_step_to_use >= a_end_time) {
        Data<IRL::SeparatorVariant> ref_interface(&cc_mesh);
        Data<IRL::VolumeMoments> ref_liq_moments(&cc_mesh);
        Data<IRL::VolumeMoments> ref_gas_moments(&cc_mesh);
        SimulationType::initialize(&velU, &velV, &velW, &ref_interface,
                                   a_end_time, a_end_time);
        setPhaseQuantities(ref_interface, &ref_liq_moments, &ref_gas_moments);
        printError(cc_mesh, liq_moments, ref_liq_moments);
      }
    }

    auto advect_end = std::chrono::system_clock::now();
    advect_VOF_time = advect_end - start;
    getReconstruction(a_reconstruction_method, liq_moments, gas_moments,
                      time_step_to_use, velU, velV, velW, &interface);
    auto recon_end = std::chrono::system_clock::now();
    recon_time = recon_end - advect_end;

    for (int i = cc_mesh.imin(); i <= cc_mesh.imax(); ++i) {
      for (int j = cc_mesh.jmin(); j <= cc_mesh.jmax(); ++j) {
        for (int k = cc_mesh.kmin(); k <= cc_mesh.kmax(); ++k) {
          const double vfrac =
              liq_moments(i, j, k).volume() / cc_mesh.cell_volume();
          if (vfrac >= IRL::global_constants::VF_LOW &&
              vfrac <= IRL::global_constants::VF_HIGH) {
            const IRL::Pt lower_pt(cc_mesh.x(i), cc_mesh.y(j), cc_mesh.z(k));
            const IRL::Pt upper_pt(cc_mesh.x(i + 1), cc_mesh.y(j + 1),
                                   cc_mesh.z(k + 1));
            const auto cuboid_cell =
                IRL::RectangularCuboid::fromBoundingPts(lower_pt, upper_pt);

            if (const IRL::PlanarSeparator* separator =
                    std::get_if<IRL::PlanarSeparator>(
                        &(interface(i, j, k)))) {  // If plane

            } else if (const IRL::Paraboloid* paraboloid =
                           std::get_if<IRL::Paraboloid>(
                               &(interface(i, j, k)))) {  // If paraboloid
              using VolumeAndSurface = IRL::AddSurfaceOutput<
                  IRL::Volume, IRL::ParaboloidParametrizedSurfaceOutput>;
              auto surface_and_volume = IRL::getVolumeMoments<VolumeAndSurface>(
                  cuboid_cell, *paraboloid);
              auto surface = surface_and_volume.getSurface();
              auto surface_area = surface.getSurfaceArea();
              surfactant_concentration(i, j, k) =
                  surfactant_mass(i, j, k) / IRL::safelyEpsilon(surface_area);
            }
          }
        }
      }
    }

    if (a_viz_frequency > 0 && iteration % a_viz_frequency == 0) {
      if (rank == 0) {
        vtk_io.writeVTKFile(simulation_time);
        writeInterfaceToFile(liq_moments, interface, surfactant_concentration,
                             simulation_time, &vtk_io,
                             simulation_time + time_step_to_use >= a_end_time);
      }
    }
    auto write_end = std::chrono::system_clock::now();
    write_time = write_end - recon_end;

    simulation_time += time_step_to_use;
    ++iteration;
  }

  return 0;
}

#endif  // EXAMPLES_SURFACE_ADVECTOR_SOLVER_H_
