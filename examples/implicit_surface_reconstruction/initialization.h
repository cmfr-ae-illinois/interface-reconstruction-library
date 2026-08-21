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
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/implicit_surfaces.h"

#include "examples/implicit_surface_reconstruction/basic_mesh.h"
#include "examples/implicit_surface_reconstruction/data.h"
#include "examples/implicit_surface_reconstruction/sparse_moments.h"
#include "examples/implicit_surface_reconstruction/surface_select.h"

// for parallelization
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

struct CellStatusStats {
  std::uint64_t cells = 0;
  std::uint64_t mixed = 0;
  std::uint64_t inside = 0;
  std::uint64_t outside = 0;
  double time = 0.0;
};

template <class SurfaceType>
CellStatusStats getCellStatus(const BasicMesh& mesh, const SurfaceType& surface,
                              InsideCellMask* inside_cells,
                              std::vector<std::uint32_t>* mixed_cell_indices);

template <class SurfaceType, std::size_t VM_ORDER, std::size_t SM_ORDER>
std::vector<SparseMixedCellMoments<VM_ORDER, SM_ORDER>> getInitializedField(
    const BasicMesh& mesh,
    const std::vector<std::uint32_t>& mixed_cell_indices_root,
    const SurfaceType& surface);

template <class SurfaceType, std::size_t VM_ORDER, std::size_t SM_ORDER>
std::size_t initializeMomentsAndWriteBin(const BasicMesh& mesh,
                                         const SurfaceType& surface,
                                         const std::string& bin_path);

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void run_initialization(const std::string& shape, int Nx,
                        const std::string& output_dir);

#include "examples/implicit_surface_reconstruction/initialization.tpp"

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_H_
