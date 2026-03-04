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
#include <Eigen/Dense>
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

// #include "irl/interface_reconstruction_methods/jibben_neighborhood.h"

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

void writeInterfaceWithScalarToFile(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::SeparatorVariant>& a_liquid_gas_interface,
    std::vector<InterfaceScalarField>* a_scalar_fields, const double a_time,
    VTKOutput* a_output, const bool print);

void printError(const BasicMesh& mesh,
                const Data<IRL::VolumeMoments>& liq_moments,
                const Data<IRL::VolumeMoments>& starting_liq_moments);

void connectMesh(
    const BasicMesh& a_mesh,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface);

std::vector<double> calculatePolygonSurfaceMoments(const IRL::Polygon& polygon);

void writeMomentsToBinary(const Data<IRL::VolumeMoments>& liq_moments_1,
                          const Data<IRL::VolumeMoments>& liq_moments_2,
                          const Data<IRL::SeparatorVariant>& interfaces_1,
                          const Data<IRL::SeparatorVariant>& interfaces_2,
                          const std::string& output_dir,
                          const std::string& case_name,
                          const std::string& reconstruction_method);

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

  // storing initial interface
  Data<IRL::SeparatorVariant> starting_interface(&cc_mesh);
  starting_interface = interface;

  // outputting interfaces in respective directory
  // std::string output_dir =
  //     "/home/parinht2/Desktop/ppic "
  //     "paper/advection_convergence/run_2/interfaces";
  // std::string vtk_output_dir = output_dir + "/" + a_case_name + "_" +
  //                              a_reconstruction_method + "_" +
  //                              std::to_string(a_nx);

  // VTKOutput vtk_io(vtk_output_dir, "viz", cc_mesh);
  VTKOutput vtk_io("viz_out", "viz", cc_mesh);
  vtk_io.addData("VelocityX", velU);
  vtk_io.addData("VelocityY", velV);
  vtk_io.addData("VelocityZ", velW);
  double simulation_time = 0.0;
  int iteration = 0;

  if (rank == 0) {
    vtk_io.writeVTKFile(simulation_time);
  }

  std::vector<InterfaceScalarField> scalar_fields;
  getReconstruction(a_reconstruction_method, liq_moments, gas_moments, 0.0,
                    velU, velV, velW, &interface, &scalar_fields);

  resetMoments(link_localized_interface, &liq_moments, &gas_moments);

  if (a_reconstruction_method == "JibbenM" ||
      a_reconstruction_method == "Hybrid") {
    writeInterfaceWithScalarToFile(liq_moments, interface, &scalar_fields,
                                   simulation_time, &vtk_io, true);
  } else {
    writeInterfaceToFile(liq_moments, interface, simulation_time, &vtk_io,
                         true);
  }

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

  // jibben interface output
  // std::string output_dir =
  // "/home/parinht2/Desktop/jibben_pu_hybrid/debugging/"; int output_count = 1;

  // outputting plane neighborhood
  // if (a_reconstruction_method == "LVIRA") {
  //   int count = 0;
  //   std::vector<IRL::Polygon> polygons;
  //   std::vector<IRL::Pt> plic_centroids;
  //   std::vector<IRL::Normal> plic_normals;
  //   std::string plane_output_dir =
  //   "/home/parinht2/Desktop/plane_neighborhood/"; for (int i =
  //   cc_mesh.imin(); i <= cc_mesh.imax(); ++i) {
  //     for (int j = cc_mesh.jmin(); j <= cc_mesh.jmax(); ++j) {
  //       for (int k = cc_mesh.kmin(); k <= cc_mesh.kmax(); ++k) {
  //         const double vf =
  //             liq_moments(i, j, k).volume() / cc_mesh.cell_volume();
  //         if (vf < IRL::global_constants::VF_LOW ||
  //             vf > IRL::global_constants::VF_HIGH) {
  //           continue;
  //         }
  //         // only want to do this for one mixed cell
  //         if (count > 0) {
  //           continue;
  //         }
  //         const int nlayers = 2;
  //         for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
  //           for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
  //             for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
  //               if (const auto ptr = std::get_if<IRL::PlanarSeparator>(
  //                       &interface(ii, jj, kk))) {
  //                 const auto polygon =
  //                     IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
  //                         IRL::RectangularCuboid::fromBoundingPts(
  //                             IRL::Pt(cc_mesh.x(ii), cc_mesh.y(jj),
  //                                     cc_mesh.z(kk)),
  //                             IRL::Pt(cc_mesh.x(ii + 1), cc_mesh.y(jj + 1),
  //                                     cc_mesh.z(kk + 1))),
  //                         *ptr, (*ptr)[0]);
  //                 IRL::Pt centroid = polygon.calculateCentroid();
  //                 IRL::Normal normal =
  //                 polygon.getPlaneOfExistence().normal();
  //                 polygons.push_back(polygon);
  //                 plic_centroids.push_back(centroid);
  //                 plic_normals.push_back(normal);
  //               }
  //             }
  //           }
  //         }
  //         count++;
  //       }
  //     }
  //   }
  //   writePolygonsVTK(plane_output_dir + "plane_neighborhood.vtk", polygons);
  //   writeScatterVTK(plic_centroids,
  //                   plane_output_dir + "plane_neighborhood_scatter.vtk");
  //   writeVectorsVTK(plane_output_dir + "plane_neighborhood_normals.vtk",
  //                   plic_centroids, plic_normals);
  // }

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
        // storing final interface
        // Data<IRL::SeparatorVariant> ending_interface(&cc_mesh);
        // const auto ending_liq_moments = liq_moments;
        // getReconstruction(a_reconstruction_method, liq_moments, gas_moments,
        //                   time_step_to_use, velU, velV, velW,
        //                   &ending_interface);
        // // writing moments to file
        // std::string moments_output_dir =
        //     "/home/parinht2/Desktop/ppic "
        //     "paper/advection_convergence/run_2/binary_files";
        // writeMomentsToBinary(starting_liq_moments, ending_liq_moments,
        //                      starting_interface, ending_interface,
        //                      moments_output_dir, a_case_name,
        //                      a_reconstruction_method);

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

    std::vector<InterfaceScalarField> scalar_fields;
    getReconstruction(a_reconstruction_method, liq_moments, gas_moments,
                      time_step_to_use, velU, velV, velW, &interface,
                      &scalar_fields);

    // outputting interfaces based on Jibben metrics

    auto recon_end = std::chrono::system_clock::now();
    recon_time = recon_end - advect_end;

    if (a_viz_frequency > 0 && iteration % a_viz_frequency == 0) {
      if (rank == 0) {
        vtk_io.writeVTKFile(simulation_time);
      }
      if (a_reconstruction_method == "JibbenM" ||
          a_reconstruction_method == "Hybrid") {
        writeInterfaceWithScalarToFile(
            liq_moments, interface, &scalar_fields, simulation_time, &vtk_io,
            simulation_time + time_step_to_use >= a_end_time);
      } else {
        writeInterfaceToFile(liq_moments, interface, simulation_time, &vtk_io,
                             simulation_time + time_step_to_use >= a_end_time);
      }
    }
    auto write_end = std::chrono::system_clock::now();
    write_time = write_end - recon_end;

    simulation_time += time_step_to_use;
    ++iteration;

    // if (iteration == 113) {
    //   std::vector<IRL::Polygon> polygons;
    //   std::vector<IRL::ParaboloidParametrizedSurfaceOutput> paraboloids;
    //   std::vector<IRL::CylinderParametrizedSurfaceOutput> cylinders;
    //   if (const auto ptr =
    //           std::get_if<IRL::Paraboloid>(&(interface(35, 34, 28)))) {
    //     std::cout << "Paraboloid detected" << std::endl;
    //     IRL::Paraboloid paraboloid = *ptr;
    //     using VolumeAndParaboloid =
    //         IRL::AddSurfaceOutput<IRL::Volume,
    //                               IRL::ParaboloidParametrizedSurfaceOutput>;
    //     auto cell = IRL::RectangularCuboid::fromBoundingPts(
    //         IRL::Pt(cc_mesh.x(34), cc_mesh.y(33), cc_mesh.z(27)),
    //         IRL::Pt(cc_mesh.x(36 + 1), cc_mesh.y(35 + 1), cc_mesh.z(29 +
    //         1)));
    //     auto volume_and_surface =
    //         IRL::getVolumeMoments<VolumeAndParaboloid>(cell, *ptr);
    //     auto surface = volume_and_surface.getSurface();
    //     paraboloids.push_back(surface);
    //   }
    //   // outputting paraboloid
    //   std::cout << "Number of polygons = " << polygons.size() << std::endl;
    //   std::cout << "Number of paraboloids = " << paraboloids.size()
    //             << std::endl;
    //   std::cout << "Number of cylinders = " << cylinders.size() << std::endl;
    //   std::string output =
    //       "/home/parinht2/Desktop/jibben_pu_hybrid/debugging/paraboloid";
    //   VTKOutput vtk_parab(output, "viz", cc_mesh);
    //   vtk_parab.writeVTKFile(simulation_time);
    //   vtk_parab.writeVTKInterface(
    //       simulation_time, polygons, paraboloids, cylinders,
    //       simulation_time + time_step_to_use >= a_end_time);
    // }

    // if (iteration == 113) {
    //   break;
    // }
  }

  return 0;
}

#endif  // EXAMPLES_VARIANT_ADVECTOR_SOLVER_H_
