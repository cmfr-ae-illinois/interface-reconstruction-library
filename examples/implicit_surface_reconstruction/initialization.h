// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_H_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_H_

#include <mpi.h>
#include <stdio.h>
#include <sys/stat.h>
#include <Eigen/Dense>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>

#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/implicit_surfaces.h"

#include "examples/implicit_surface_reconstruction/basic_mesh.h"
#include "examples/implicit_surface_reconstruction/data.h"

std::vector<std::tuple<int, int, int>> getCellStatus(Data<int>* cell_status);

struct Range {
  int begin;
  int end;
};

inline Range block_partition(int N, int rank, int size) {
  const int q = N / size;
  const int r = N % size;
  const int start = rank * q + std::min(rank, r);
  const int count = q + (rank < r ? 1 : 0);
  return {start, start + count};
}

template <std::size_t MAX_REFINE, std::size_t VM_ORDER, std::size_t SM_ORDER>
void getInitializedField(
    const Data<int>& cell_status,
    std::vector<std::tuple<int, int, int>> mixed_cells_list_root,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* moments);

#include "examples/implicit_surface_reconstruction/initialization.tpp"

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_H_
