// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_SOLVER_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_SOLVER_H_

#include <mpi.h>
#include <sys/stat.h>
#include <chrono>
#include <iostream>
#include <string>

#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/planar_reconstruction/planar_localizer.h"

#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/deformation_3d.h"
#include "examples/paraboloid_advector/reconstruction_types.h"
#include "examples/paraboloid_advector/rotation_3d.h"
#include "examples/paraboloid_advector/translation_3d.h"
#include "examples/paraboloid_advector/vof_advection.h"
#include "examples/paraboloid_advector/vtk.h"

/// \brief Handles running and advancing the solution according to provided
/// static functions in structs.
template <class SimulationType>
int runSimulation(const std::string& a_simulation_type,
                  const std::string& a_advection_method,
                  const std::string& a_reconstruction_method, const double a_dt,
                  const double a_end_time, const int a_visualization_frequency);

// \brief Convert and store the mesh cells into localizers.
void initializeLocalizers(Data<IRL::PlanarLocalizer>* a_localizers);

/// \brief Initialize the linked localized paraboloids used during advection.
void initializeLocalizedParaboloids(
    const Data<IRL::PlanarLocalizer>& a_cell_localizers,
    const Data<IRL::Paraboloid>& a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>* a_linked_localized_paraboloid);

/// \brief Set phase quantities according to the given
void setPhaseQuantities(const Data<IRL::Paraboloid>& a_interface,
                        Data<IRL::GeneralMoments3D<2>>* a_liquid_moments,
                        Data<IRL::GeneralMoments3D<2>>* a_gas_moments);

/// \brief Write out the header for the diagnostics.
void writeDiagnosticsHeader(void);

/// \brief Write out diagnostics to the screen.
void writeOutDiagnostics(const int a_iteration, const double a_dt,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<double>& a_W,
                         const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
                         const Data<IRL::Paraboloid>& a_interface,
                         std::chrono::duration<double> a_VOF_duration,
                         std::chrono::duration<double> a_recon_duration,
                         std::chrono::duration<double> a_write_duration);

/// \brief Generates triangulated surface and writes to provided VTK file
void writeInterfaceToFile(
    const Data<IRL::GeneralMoments3D<2>>& a_liquid_moments,
    const Data<IRL::Paraboloid>& a_liquid_gas_interface, const double a_time,
    VTKOutput* a_output, const bool print);

void printError(const BasicMesh& mesh,
                const Data<IRL::GeneralMoments3D<2>>& liquid_moments,
                const Data<IRL::GeneralMoments3D<2>>& starting_liquid_moments);

//******************************************************************* //
//     Template function definitions placed below this.
//******************************************************************* //
template <class SimulationType>
int runSimulation(const std::string& a_simulation_type,
                  const std::string& a_advection_method,
                  const std::string& a_reconstruction_method, const double a_dt,
                  const double a_end_time, const int a_visualization_frequency,
                  const int a_nx) {
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Set mesh
  BasicMesh cc_mesh = SimulationType::setMesh(a_nx);

  // Allocate local data
  Data<double> velU(&cc_mesh);
  Data<double> velV(&cc_mesh);
  Data<double> velW(&cc_mesh);
  Data<IRL::GeneralMoments3D<2>> liquid_moments(&cc_mesh);
  Data<IRL::GeneralMoments3D<2>> gas_moments(&cc_mesh);

  // Allocate
  Data<IRL::PlanarLocalizer> cell_localizers(&cc_mesh);
  initializeLocalizers(&cell_localizers);
  Data<IRL::Paraboloid> interface(&cc_mesh);
  Data<IRL::LocalizedParaboloidLink<double>> link_localized_paraboloids(
      &cc_mesh);
  initializeLocalizedParaboloids(cell_localizers, interface,
                                 &link_localized_paraboloids);
  connectMesh(cc_mesh, &link_localized_paraboloids);
  // Set constants in IRL
  IRL::setMinimumVolumeToTrack(10.0 * DBL_EPSILON * cc_mesh.dx() *
                               cc_mesh.dy() * cc_mesh.dz());
  IRL::setVolumeFractionBounds(1.0e-13);
  IRL::setVolumeFractionTolerance(1.0e-12);

  if (rank == 0) {
    std::cout << "VolumeFractionBounds = " << IRL::global_constants::VF_LOW
              << std::endl;
  }

  // Initialize data
  SimulationType::initialize(&velU, &velV, &velW, &interface, 0.0);
  setPhaseQuantities(interface, &liquid_moments, &gas_moments);
  const auto starting_liquid_moments = liquid_moments;

  VTKOutput vtk_io("viz_out", "viz", cc_mesh);
  vtk_io.addData("VelocityX", velU);
  vtk_io.addData("VelocityY", velV);
  vtk_io.addData("VelocityZ", velW);
  double simulation_time = 0.0;
  int iteration = 0;

  if (rank == 0) {
    vtk_io.writeVTKFile(simulation_time);
  }
  getReconstruction(a_reconstruction_method, liquid_moments, gas_moments,
                    &link_localized_paraboloids, 0.0, velU, velV, velW,
                    &interface);
  resetMoments(link_localized_paraboloids, &liquid_moments, &gas_moments);

  writeInterfaceToFile(liquid_moments, interface, simulation_time, &vtk_io,
                       true);
  if (rank == 0) {
    writeDiagnosticsHeader();
  }
  std::string output_folder = "viz";
  const int dir_err = mkdir(output_folder.c_str(), 0777);
  std::chrono::duration<double> advect_VOF_time(0.0);
  std::chrono::duration<double> recon_time(0.0);
  std::chrono::duration<double> write_time(0.0);
  if (rank == 0) {
    writeOutDiagnostics(iteration, a_dt, simulation_time, velU, velV, velW,
                        liquid_moments, interface, advect_VOF_time, recon_time,
                        write_time);
    printError(cc_mesh, liquid_moments, starting_liquid_moments);
  }
  while (simulation_time < a_end_time) {
    const double time_step_to_use =
        std::fmin(a_dt, a_end_time - simulation_time);
    SimulationType::setVelocity(simulation_time + 0.5 * time_step_to_use, &velU,
                                &velV, &velW);

    auto start = std::chrono::system_clock::now();
    advectVOF(a_simulation_type, a_advection_method, a_reconstruction_method,
              time_step_to_use, simulation_time, velU, velV, velW,
              &link_localized_paraboloids, &liquid_moments, &gas_moments,
              &interface);

    if (rank == 0) {
      writeOutDiagnostics(iteration + 1, time_step_to_use,
                          simulation_time + time_step_to_use, velU, velV, velW,
                          liquid_moments, interface, advect_VOF_time,
                          recon_time, write_time);
      if (simulation_time + time_step_to_use >= a_end_time) {
        Data<IRL::Paraboloid> ref_interface(&cc_mesh);
        Data<IRL::GeneralMoments3D<2>> ref_liquid_moments(&cc_mesh);
        Data<IRL::GeneralMoments3D<2>> ref_gas_moments(&cc_mesh);
        SimulationType::initialize(&velU, &velV, &velW, &ref_interface,
                                   simulation_time + time_step_to_use);
        setPhaseQuantities(ref_interface, &ref_liquid_moments,
                           &ref_gas_moments);
        printError(cc_mesh, liquid_moments, ref_liquid_moments);
      }
    }

    auto advect_end = std::chrono::system_clock::now();
    advect_VOF_time = advect_end - start;
    getReconstruction(a_reconstruction_method, liquid_moments, gas_moments,
                      &link_localized_paraboloids, time_step_to_use, velU, velV,
                      velW, &interface);
    auto recon_end = std::chrono::system_clock::now();
    recon_time = recon_end - advect_end;

    if (a_visualization_frequency > 0 &&
        iteration % a_visualization_frequency == 0) {
      if (rank == 0) {
        vtk_io.writeVTKFile(simulation_time);
      }
      writeInterfaceToFile(liquid_moments, interface, simulation_time, &vtk_io,
                           simulation_time + time_step_to_use >= a_end_time);
    }
    auto write_end = std::chrono::system_clock::now();
    write_time = write_end - recon_end;

    simulation_time += time_step_to_use;
    ++iteration;
    // }
  }

  return 0;
}

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_SOLVER_H_
