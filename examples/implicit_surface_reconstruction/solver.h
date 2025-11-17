// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_SOLVER_H_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_SOLVER_H_

#include <mpi.h>
#include <stdio.h>
#include <sys/stat.h>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>
#include "Eigen/Dense"

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/implicit_surfaces.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include "examples/implicit_surface_reconstruction/basic_mesh.h"
#include "examples/implicit_surface_reconstruction/data.h"
#include "examples/implicit_surface_reconstruction/reconstruction_types.h"
#include "examples/implicit_surface_reconstruction/vtk.h"

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

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_SOLVER_H_
