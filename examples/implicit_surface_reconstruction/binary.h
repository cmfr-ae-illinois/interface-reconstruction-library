// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_BINARY_H_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_BINARY_H_

#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.h"

#include "examples/implicit_surface_reconstruction/basic_mesh.h"
#include "examples/implicit_surface_reconstruction/data.h"

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void writeMomentsToBinary(
    const Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                         IRL::GeneralSurfaceMoments3D<SM_ORDER>>>& moments,
    const std::string& filename);

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void readMomentsFromBinary(
    const std::string& filename,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* moments);

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void coarsenMomentsFromBinary(
    const std::string& fine_filename, int factor,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* coarse_moments);

#include "examples/implicit_surface_reconstruction/binary.tpp"

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_BINARY_H_
