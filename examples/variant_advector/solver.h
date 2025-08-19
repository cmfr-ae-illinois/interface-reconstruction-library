// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_VARIANT_ADVECTOR_SOLVER_H_
#define EXAMPLES_VARIANT_ADVECTOR_SOLVER_H_

#include <mpi.h>
#include <stdio.h>
#include <sys/stat.h>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include "examples/variant_advector/basic_mesh.h"
#include "examples/variant_advector/data.h"
#include "examples/variant_advector/reconstruction_types.h"
#include "examples/variant_advector/vof_advection.h"
#include "examples/variant_advector/vtk.h"

#include "examples/variant_advector/translation_3d.h"

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
    const double a_time, VTKOutput* a_output, const bool print);

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

  writeInterfaceToFile(liq_moments, interface, simulation_time, &vtk_io, true);
  if (rank == 0) {
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
              simulation_time, velU, velV, velW, &link_localized_interface,
              &liq_moments, &gas_moments);

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

    if (a_viz_frequency > 0 && iteration % a_viz_frequency == 0) {
      if (rank == 0) {
        vtk_io.writeVTKFile(simulation_time);
      }
      writeInterfaceToFile(liq_moments, interface, simulation_time, &vtk_io,
                           simulation_time + time_step_to_use >= a_end_time);
    }
    auto write_end = std::chrono::system_clock::now();
    write_time = write_end - recon_end;

    simulation_time += time_step_to_use;
    ++iteration;
  }

  return 0;
}

#endif  // EXAMPLES_VARIANT_ADVECTOR_SOLVER_H_
