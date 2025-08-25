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
#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/implicit_surfaces.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include "examples/sphere_initialization/basic_mesh.h"
#include "examples/sphere_initialization/data.h"
#include "examples/sphere_initialization/reconstruction_types.h"
#include "examples/sphere_initialization/vtk.h"

void reconstructSurface(const std::string& a_reconstruction_method,
                        const int a_nx);

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

/// \brief Generates triangulated surface and writes to provided VTK file
void writeInterfaceToFile(
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::SeparatorVariant>& a_liquid_gas_interface,
    const double a_time, VTKOutput* a_output, const bool print);

void connectMesh(
    const BasicMesh& a_mesh,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface);

//******************************************************************* //
//     Function definitions placed below this.
//******************************************************************* //
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

#endif  // EXAMPLES_VARIANT_ADVECTOR_SOLVER_H_
