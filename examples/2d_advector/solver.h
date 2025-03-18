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

#include <sys/stat.h>
#include <chrono>
#include <iostream>
#include <string>

#include "examples/2d_advector/basic_mesh.h"
#include "examples/2d_advector/data.h"
#include "examples/2d_advector/irl2d.h"
#include "examples/2d_advector/reconstruction_types.h"
#include "examples/2d_advector/rotation_2d.h"
#include "examples/2d_advector/vof_advection.h"
#include "examples/2d_advector/vtk.h"

/// \brief Handles running and advancing the solution according to provided
/// static functions in structs.
template <class SimulationType>
int runSimulation(const std::string& a_simulation_type,
                  const std::string& a_advection_method,
                  const std::string& a_reconstruction_method, const double a_dt,
                  const double a_end_time, const int a_visualization_frequency);

/// \brief Set phase quantities according to the given
void setPhaseQuantities(const Data<IRL2D::Parabola>& a_interface,
                        Data<IRL2D::Moments>* a_liquid_moments,
                        Data<IRL2D::Moments>* a_gas_moments,
                        Data<double>* a_vfrac);

/// \brief Write out the header for the diagnostics.
void writeDiagnosticsHeader(void);

/// \brief Write out diagnostics to the screen.
void writeOutDiagnostics(const int a_iteration, const double a_dt,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<IRL2D::Moments>& a_liquid_moments,
                         const Data<IRL2D::Parabola>& a_interface,
                         std::chrono::duration<double> a_VOF_duration,
                         std::chrono::duration<double> a_recon_duration,
                         std::chrono::duration<double> a_write_duration);

/// \brief Generates triangulated surface and writes to provided VTK file
void writeInterfaceToFile(const Data<IRL2D::Moments>& a_liquid_moments,
                          const Data<IRL2D::Parabola>& a_liquid_gas_interface,
                          const double a_time, VTKOutput* a_output,
                          const bool print);

void printError(const BasicMesh& mesh,
                const Data<IRL2D::Moments>& liquid_moments,
                const Data<IRL2D::Moments>& starting_liquid_moments);

//******************************************************************* //
//     Template function definitions placed below this.
//******************************************************************* //
template <class SimulationType>
int runSimulation(const std::string& a_simulation_type,
                  const std::string& a_advection_method,
                  const std::string& a_reconstruction_method, const double a_dt,
                  const double a_end_time, const int a_visualization_frequency,
                  const int a_nx) {
#ifdef USE_MPI
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  // Set mesh
  BasicMesh cc_mesh = SimulationType::setMesh(a_nx);

  // Allocate local data
  Data<double> velU(&cc_mesh);
  Data<double> velV(&cc_mesh);
  Data<double> vfrac(&cc_mesh);
  Data<IRL2D::Moments> liquid_moments(&cc_mesh);
  Data<IRL2D::Moments> gas_moments(&cc_mesh);
  Data<IRL2D::Parabola> interface(&cc_mesh);

  // Set constants in IRL
  IRL::setMinimumVolumeToTrack(10.0 * DBL_EPSILON * cc_mesh.dx() *
                               cc_mesh.dy());
  IRL::setVolumeFractionBounds(1.0e-9);
  IRL::setVolumeFractionTolerance(1.0e-9);

#ifdef USE_MPI
  if (!rank) {
#endif
    std::cout << "VolumeFractionBounds = " << IRL::global_constants::VF_LOW
              << std::endl;
#ifdef USE_MPI
  }
#endif

  // Initialize data
  SimulationType::initialize(&velU, &velV, &interface, 0.0);
  setPhaseQuantities(interface, &liquid_moments, &gas_moments, &vfrac);
  const auto starting_liquid_moments = liquid_moments;

  // cell coordinates
  std::cout << "CELL COORDINATES -----------------------" << std::endl;
  std::cout << "x0 = " << IRL2D::Vec(cc_mesh.x(23), cc_mesh.y(51)) << std::endl;
  std::cout << "x1 = " << IRL2D::Vec(cc_mesh.x(24), cc_mesh.y(52)) << std::endl;

  // exact interface parameters
  std::cout << "INTERFACE PARAMETERS -----------------------" << std::endl;
  std::cout << "Exact coefficient = " << interface(23,51).coeff() << std::endl;
  std::cout << "Exact datum = " << interface(23,51).datum() << std::endl;
  std::cout << "Exact reference frame = " << interface(23,51).frame() << std::endl;

  // exact moments of a cell
  std::cout << "LIQUID MOMENTS -----------------------" << std::endl;
  std::cout << "m0 = " << liquid_moments(23,51).m0() << std::endl;
  std::cout << "m1 = " << liquid_moments(23,51).m1() << std::endl;
  std::cout << "m2 = " << liquid_moments(23,51).m2() << std::endl;
  std::cout << "GAS MOMENTS -----------------------" << std::endl;
  std::cout << "m0 = " << gas_moments(23,51).m0() << std::endl;
  std::cout << "m1 = " << gas_moments(23,51).m1() << std::endl;
  std::cout << "m2 = " << gas_moments(23,51).m2() << std::endl;
  
  // exact volume fraction and centroid
  std::cout << "VOLFRAC & CENTROID -----------------------" << std::endl;
  std::cout << "liquid vol frac = " << liquid_moments(23,51).m0()/cc_mesh.cell_volume() << std::endl;
  std::cout << "gas vol frac = " << gas_moments(23,51).m0()/cc_mesh.cell_volume() << std::endl;
  std::cout << "liquid centroid = " << liquid_moments(23,51).m1()/liquid_moments(23,51).m0() << std::endl;
  std::cout << "gas centroid = " << gas_moments(23,51).m1()/gas_moments(23,51).m0() << std::endl;

  // std::cout << "------------------------------------------------------------------------------" << std::endl;

  VTKOutput vtk_io("viz_out", "viz", cc_mesh);
  vtk_io.addData("VelocityX", velU);
  vtk_io.addData("VelocityY", velV);
  vtk_io.addData("VFrac", vfrac);
  double simulation_time = 0.0;
  int iteration = 0;

  vtk_io.writeVTKFile(simulation_time);
  getReconstruction(a_reconstruction_method, liquid_moments, gas_moments, 0.0,
                    velU, velV, &interface);
  
  // for (int i = cc_mesh.imin(); i <= cc_mesh.imax(); ++i){
  //   for (int j = cc_mesh.jmin(); j <= cc_mesh.jmax(); ++j){
  //     if (i != 25 && j != 51){
  //       interface(i,j) = IRL2D::Parabola();
  //     }
  //   }
  // }

    // // checking how many planes are there and where
  // int plane_count = 0;
  // for (int i = cc_mesh.imin(); i <= cc_mesh.imax(); ++i){
  //   for (int j = cc_mesh.jmin(); j <= cc_mesh.jmax(); ++j){
  //     double liq_vol_frac = liquid_moments(i,j).m0() / cc_mesh.cell_volume();
  //     if (liq_vol_frac >= IRL::global_constants::VF_LOW &&
  //         liq_vol_frac <= IRL::global_constants::VF_HIGH){
  //           double alpha = interface(i,j).coeff();
  //           //if (std::abs(alpha) < 1.0e-1){
  //             plane_count++;
  //             std::cout << "(i,j) = (" << i << "," << j << "), alpha = " << alpha << std::endl;
  //           //}
  //         }
  //   }
  // }
  // std::cout << "Cells with plane = " << plane_count << std::endl;

  setPhaseQuantities(interface, &liquid_moments, &gas_moments, &vfrac);
  
  std::cout << "---------------------------------------------------------------" << std::endl;
  // moments after reconstruction
  // exact moments of a cell
  // std::cout << "new liquid m0 = " << liquid_moments(23,51).m0() << std::endl;
  // std::cout << "new liquid m1 = " << liquid_moments(23,51).m1() << std::endl;
  // std::cout << "new liquid m2 = " << liquid_moments(23,51).m2() << std::endl;
  // std::cout << "new gas m0 = " << gas_moments(23,51).m0() << std::endl;
  // std::cout << "new gas m1 = " << gas_moments(23,51).m1() << std::endl;
  // std::cout << "new gas m2 = " << gas_moments(23,51).m2() << std::endl;

  // // new volume fraction and centroid
  // std::cout << "new liquid vol frac = " << liquid_moments(23,51).m0()/cc_mesh.cell_volume() << std::endl;
  // std::cout << "new gas vol frac = " << gas_moments(23,51).m0()/cc_mesh.cell_volume() << std::endl;
  // std::cout << "new liquid centroid = " << liquid_moments(23,51).m1()/liquid_moments(23,51).m0() << std::endl;
  // std::cout << "new gas centroid = " << gas_moments(23,51).m1()/gas_moments(23,51).m0() << std::endl;

  // // interface parameters after reconstruction of circle
  // std::cout << "new coefficient = " << interface(23,51).coeff() << std::endl;
  // std::cout << "new datum = " << interface(23,51).datum() << std::endl;
  // std::cout << "new reference frame = " << interface(23,51).frame() << std::endl;

  std::cout << "INTERFACE PARAMETERS -----------------------" << std::endl;
  std::cout << "coefficient = " << interface(23,51).coeff() << std::endl;
  std::cout << "datum = " << interface(23,51).datum() << std::endl;
  std::cout << "reference frame = " << interface(23,51).frame() << std::endl;

  // exact moments of a cell
  std::cout << "LIQUID MOMENTS -----------------------" << std::endl;
  std::cout << "m0 = " << liquid_moments(23,51).m0() << std::endl;
  std::cout << "m1 = " << liquid_moments(23,51).m1() << std::endl;
  std::cout << "m2 = " << liquid_moments(23,51).m2() << std::endl;
  std::cout << "GAS MOMENTS -----------------------" << std::endl;
  std::cout << "m0 = " << gas_moments(23,51).m0() << std::endl;
  std::cout << "m1 = " << gas_moments(23,51).m1() << std::endl;
  std::cout << "m2 = " << gas_moments(23,51).m2() << std::endl;
  
  // exact volume fraction and centroid
  std::cout << "VOLFRAC & CENTROID -----------------------" << std::endl;
  std::cout << "liquid vol frac = " << liquid_moments(23,51).m0()/cc_mesh.cell_volume() << std::endl;
  std::cout << "gas vol frac = " << gas_moments(23,51).m0()/cc_mesh.cell_volume() << std::endl;
  std::cout << "liquid centroid = " << liquid_moments(23,51).m1()/liquid_moments(23,51).m0() << std::endl;
  std::cout << "gas centroid = " << gas_moments(23,51).m1()/gas_moments(23,51).m0() << std::endl;

  std::cout << "--------------------------------------------------------------------------" << std::endl;

  vtk_io.writeVTKInterface(simulation_time, interface);
  std::string output_folder = "viz";
  const int dir_err = mkdir(output_folder.c_str(), 0777);
  std::chrono::duration<double> advect_VOF_time(0.0);
  std::chrono::duration<double> recon_time(0.0);
  std::chrono::duration<double> write_time(0.0);
#ifdef USE_MPI
  if (!rank) {
#endif
    writeDiagnosticsHeader();
    writeOutDiagnostics(iteration, a_dt, simulation_time, velU, velV,
                        liquid_moments, interface, advect_VOF_time, recon_time,
                        write_time);
#ifdef USE_MPI
  }
#endif
  // printError(cc_mesh, liquid_moments, starting_liquid_moments);
  while (simulation_time < a_end_time) {
    const double time_step_to_use =
        std::fmin(a_dt, a_end_time - simulation_time);
    SimulationType::setVelocity(simulation_time + 0.5 * time_step_to_use, &velU,
                                &velV);

    auto start = std::chrono::system_clock::now();
    advectVOF(a_simulation_type, a_advection_method, a_reconstruction_method,
              time_step_to_use, simulation_time, velU, velV, &liquid_moments,
              &gas_moments, &interface);

#ifdef USE_MPI
    if (!rank) {
#endif
      writeOutDiagnostics(iteration + 1, time_step_to_use,
                          simulation_time + time_step_to_use, velU, velV,
                          liquid_moments, interface, advect_VOF_time,
                          recon_time, write_time);
#ifdef USE_MPI
    }
#endif

    auto advect_end = std::chrono::system_clock::now();
    advect_VOF_time = advect_end - start;
    getReconstruction(a_reconstruction_method, liquid_moments, gas_moments,
                      time_step_to_use, velU, velV, &interface);
    auto recon_end = std::chrono::system_clock::now();
    recon_time = recon_end - advect_end;

    if (a_visualization_frequency > 0 &&
        ((iteration + 1) % a_visualization_frequency == 0 ||
         time_step_to_use == a_end_time - simulation_time)) {
      for (int i = cc_mesh.imin(); i <= cc_mesh.imax(); ++i) {
        for (int j = cc_mesh.jmin(); j <= cc_mesh.jmax(); ++j) {
          vfrac(i, j) = liquid_moments(i, j).m0() / cc_mesh.cell_volume();
        }
      }
      vtk_io.writeVTKFile(simulation_time);
      vtk_io.writeVTKInterface(simulation_time, interface);
    }

    // setPhaseQuantities(interface, &liquid_moments, &gas_moments, &vfrac);
    auto write_end = std::chrono::system_clock::now();
    write_time = write_end - recon_end;

    simulation_time += time_step_to_use;
    ++iteration;
  }

  vtk_io.writeVTKFile(simulation_time);
  vtk_io.writeVTKInterface(simulation_time, interface);

  Data<double> ref_vfrac(&cc_mesh);
  Data<IRL2D::Parabola> ref_interface(&cc_mesh);
  Data<IRL2D::Moments> ref_liquid_moments(&cc_mesh);
  Data<IRL2D::Moments> ref_gas_moments(&cc_mesh);
  SimulationType::initialize(&velU, &velV, &ref_interface, a_end_time);
  setPhaseQuantities(ref_interface, &ref_liquid_moments, &ref_gas_moments,
                     &ref_vfrac);
  getReconstruction(a_reconstruction_method, ref_liquid_moments,
                    ref_gas_moments, a_end_time, velU, velV, &ref_interface);
  setPhaseQuantities(ref_interface, &ref_liquid_moments, &ref_gas_moments,
                     &ref_vfrac);
#ifdef USE_MPI
  if (!rank) {
#endif
    printError(cc_mesh, liquid_moments, ref_liquid_moments);
#ifdef USE_MPI
  }
#endif

  return 0;
}

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_SOLVER_H_
