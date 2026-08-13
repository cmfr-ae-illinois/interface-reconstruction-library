// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_TPP_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_TPP_

#include "examples/implicit_surface_reconstruction/initialization.h"

// wrapper to set refine level to be 0 for a surface
template <class S>
struct ZeroRefineSurface : S {
  using S::S;
  ZeroRefineSurface(const S& s) : S(s) {}
  static constexpr std::size_t getMaxRefineLevel() { return 0; }
};

struct SizeRange {
  std::size_t begin;
  std::size_t end;
};

inline SizeRange blockPartitionSize(const std::size_t count, const int rank,
                                    const int size) {
  const std::size_t size_value = static_cast<std::size_t>(size);
  const std::size_t rank_value = static_cast<std::size_t>(rank);
  const std::size_t quotient = count / size_value;
  const std::size_t remainder = count % size_value;
  const std::size_t begin =
      rank_value * quotient + std::min(rank_value, remainder);
  const std::size_t local_count = quotient + (rank_value < remainder ? 1 : 0);
  return {begin, begin + local_count};
}

// finding mixed cells ----------------------------------------------
template <class SurfaceType>
void getCellStatus(const BasicMesh& mesh, const SurfaceType& surface,
                   InsideCellMask* inside_cells,
                   std::vector<std::uint32_t>* mixed_cell_indices) {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  using S0 = ZeroRefineSurface<std::decay_t<SurfaceType>>;
  S0 surface0(surface);

  if (inside_cells == nullptr || mixed_cell_indices == nullptr) {
    throw std::invalid_argument("Sparse cell-status outputs cannot be null");
  }

  const std::size_t cell_count =
      static_cast<std::size_t>(mesh.getNx()) * mesh.getNy() * mesh.getNz();
  if (cell_count == 0 ||
      cell_count - 1 >
          static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw std::overflow_error("Cell index exceeds sparse file format");
  }

  // Partition complete 64-bit mask words so local masks can be concatenated
  // directly on rank 0 without allocating the global mask on every rank.
  const std::size_t total_word_count = (cell_count + 63) / 64;
  const SizeRange local_word_range =
      blockPartitionSize(total_word_count, rank, size);
  const std::size_t local_word_count =
      local_word_range.end - local_word_range.begin;
  const std::size_t first_linear_index = local_word_range.begin * 64;
  const std::size_t end_linear_index =
      std::min(cell_count, local_word_range.end * 64);

  // Allocate whole local words, including harmless padding after the final
  // physical cell on the last rank.
  InsideCellMask local_inside(local_word_count * 64);
  std::vector<std::uint32_t> local_mixed_cell_indices;

  for (std::size_t index = first_linear_index; index < end_linear_index;
       ++index) {
    int i, j, k;
    getCellIndicesFromLinearIndex(index, mesh.getNy(), mesh.getNz(), &i, &j,
                                  &k);
    i += mesh.imin();
    j += mesh.jmin();
    k += mesh.kmin();

    const IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
    const IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
    const IRL::RectangularCuboid cell =
        IRL::RectangularCuboid::fromBoundingPts(x0, x1);
    IRL::ImplicitSurfaceCutter<S0, IRL::Volume> cutter(surface0, cell);

    const int status = cutter.getBaseCellStatus();
    if (status == 0) {
      local_mixed_cell_indices.push_back(static_cast<std::uint32_t>(index));
    } else if (status == -1) {
      local_inside.set(index - first_linear_index);
    }
  }

  if (local_word_count >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("Local inside mask exceeds MPI count limit");
  }

  std::vector<int> mask_counts, mask_displacements;
  if (rank == 0) {
    *inside_cells = InsideCellMask(cell_count);
    mask_counts.resize(size);
    mask_displacements.resize(size);
    for (int r = 0; r < size; ++r) {
      const SizeRange range = blockPartitionSize(total_word_count, r, size);
      const std::size_t count = range.end - range.begin;
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
          range.begin >
              static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("Global inside mask exceeds MPI limits");
      }
      mask_counts[r] = static_cast<int>(count);
      mask_displacements[r] = static_cast<int>(range.begin);
    }
  }

  MPI_Gatherv(local_inside.words().data(), static_cast<int>(local_word_count),
              MPI_UINT64_T, rank == 0 ? inside_cells->words().data() : nullptr,
              rank == 0 ? mask_counts.data() : nullptr,
              rank == 0 ? mask_displacements.data() : nullptr, MPI_UINT64_T, 0,
              MPI_COMM_WORLD);

  if (local_mixed_cell_indices.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("Local mixed-cell list exceeds MPI count limit");
  }
  const int local_mixed_count =
      static_cast<int>(local_mixed_cell_indices.size());
  std::vector<int> mixed_counts, mixed_displacements;
  if (rank == 0) mixed_counts.resize(size);
  MPI_Gather(&local_mixed_count, 1, MPI_INT,
             rank == 0 ? mixed_counts.data() : nullptr, 1, MPI_INT, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    mixed_displacements.resize(size);
    std::size_t total_mixed_count = 0;
    for (int r = 0; r < size; ++r) {
      if (total_mixed_count >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("Global mixed-cell list exceeds MPI limits");
      }
      mixed_displacements[r] = static_cast<int>(total_mixed_count);
      total_mixed_count += static_cast<std::size_t>(mixed_counts[r]);
    }
    if (total_mixed_count >
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error("Global mixed-cell list exceeds MPI limits");
    }
    mixed_cell_indices->resize(total_mixed_count);
  } else {
    mixed_cell_indices->clear();
  }

  MPI_Gatherv(local_mixed_cell_indices.data(), local_mixed_count, MPI_UINT32_T,
              rank == 0 ? mixed_cell_indices->data() : nullptr,
              rank == 0 ? mixed_counts.data() : nullptr,
              rank == 0 ? mixed_displacements.data() : nullptr, MPI_UINT32_T, 0,
              MPI_COMM_WORLD);
}

// finding implicit surface moments -----------------------------------
template <class SurfaceType, std::size_t VM_ORDER, std::size_t SM_ORDER>
std::vector<SparseMixedCellMoments<VM_ORDER, SM_ORDER>> getInitializedField(
    const BasicMesh& mesh,
    const std::vector<std::uint32_t>& mixed_cell_indices_root,
    const SurfaceType& surface) {
  int rank = 0, size = 1;
  (void)MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  (void)MPI_Comm_size(MPI_COMM_WORLD, &size);

  int mixed_count =
      rank == 0 ? static_cast<int>(mixed_cell_indices_root.size()) : 0;
  MPI_Bcast(&mixed_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> counts, displacements;
  if (rank == 0) {
    counts.resize(size);
    displacements.resize(size);
    for (int r = 0; r < size; ++r) {
      const Range range = block_partition(mixed_count, r, size);
      counts[r] = range.end - range.begin;
      displacements[r] = range.begin;
    }
  }

  const Range local_range = block_partition(mixed_count, rank, size);
  const int local_count = local_range.end - local_range.begin;
  std::vector<std::uint32_t> local_indices(local_count);
  MPI_Scatterv(rank == 0 ? mixed_cell_indices_root.data() : nullptr,
               rank == 0 ? counts.data() : nullptr,
               rank == 0 ? displacements.data() : nullptr, MPI_UINT32_T,
               local_indices.data(), local_count, MPI_UINT32_T, 0,
               MPI_COMM_WORLD);

  using Record = SparseMixedCellMoments<VM_ORDER, SM_ORDER>;
  std::vector<Record> local;
  local.reserve(local_indices.size());

  for (const std::uint32_t linear_index : local_indices) {
    int i, j, k;
    getCellIndicesFromLinearIndex(linear_index, mesh.getNy(), mesh.getNz(), &i,
                                  &j, &k);
    i += mesh.imin();
    j += mesh.jmin();
    k += mesh.kmin();

    IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
    IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
    IRL::RectangularCuboid cell =
        IRL::RectangularCuboid::fromBoundingPts(x0, x1);

    IRL::ImplicitSurfaceCutter<SurfaceType, IRL::GeneralMoments3D<VM_ORDER>>
        cutter(surface, cell);

    auto vol = cutter.computeVolumeMoments();
    auto surf = cutter.template computeSurfaceMoments<SM_ORDER>(
        false, Eigen::Integrator<double, 2>::GaussKronrod15, 5);

    Record rec{};
    rec.linear_index = linear_index;
    for (std::size_t t = 0; t < Record::NV; ++t) rec.volume[t] = vol[t];
    for (std::size_t t = 0; t < Record::NS; ++t) rec.surface[t] = surf[t];
    local.push_back(rec);
  }

  const int record_bytes = static_cast<int>(sizeof(Record));
  std::vector<int> byte_counts, byte_displacements;
  if (rank == 0) {
    byte_counts.resize(size);
    byte_displacements.resize(size);
    for (int r = 0; r < size; ++r) {
      byte_counts[r] = counts[r] * record_bytes;
      byte_displacements[r] = displacements[r] * record_bytes;
    }
  }

  std::vector<Record> gathered;
  if (rank == 0) gathered.resize(mixed_count);
  MPI_Gatherv(local.data(), local_count * record_bytes, MPI_BYTE,
              rank == 0 ? gathered.data() : nullptr,
              rank == 0 ? byte_counts.data() : nullptr,
              rank == 0 ? byte_displacements.data() : nullptr, MPI_BYTE, 0,
              MPI_COMM_WORLD);
  return gathered;
}

// initializing moments and writing to binary ------------------------------
template <class SurfaceType, std::size_t VM_ORDER, std::size_t SM_ORDER>
void initializeMomentsAndWriteBin(const BasicMesh& mesh,
                                  const SurfaceType& surface,
                                  const std::string& bin_path) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  InsideCellMask inside_cells;
  std::vector<std::uint32_t> mixed_cell_indices;
  getCellStatus(mesh, surface, &inside_cells, &mixed_cell_indices);

  auto mixed_moments = getInitializedField<SurfaceType, VM_ORDER, SM_ORDER>(
      mesh, mixed_cell_indices, surface);

  if (rank == 0) {
    writeMomentsToBinary<VM_ORDER, SM_ORDER>(mesh, inside_cells, mixed_moments,
                                             bin_path);
  }
}

// running initialization for generic shape ---------------------------
template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void run_initialization(const std::string& shape, int Nx,
                        const std::string& output_dir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  BasicMesh mesh(Nx, Nx, Nx, 1);
  SurfaceVariant surf = makeSurface(shape, mesh);
  std::string bin_path = binary_filename(output_dir, shape, Nx);

  std::visit(
      [&](const auto& surface) {
        using S = std::decay_t<decltype(surface)>;

        auto t0 = std::chrono::steady_clock::now();
        initializeMomentsAndWriteBin<S, VM_ORDER, SM_ORDER>(mesh, surface,
                                                            bin_path);
        auto t1 = std::chrono::steady_clock::now();
        if (rank == 0) {
          const double dt = std::chrono::duration<double>(t1 - t0).count();
          std::cout << "✅ Initialization complete for " << shape << " in "
                    << dt << " s\n";
        }
      },
      surf);
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_TPP_
