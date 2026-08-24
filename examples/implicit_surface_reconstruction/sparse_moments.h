// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_2_SPARSE_MOMENTS_H_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_2_SPARSE_MOMENTS_H_

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "irl/moments/general_moments.h"
#include "irl/moments/general_surface_moments.h"

// instead of boolean vector, we use a bitset to store the inside/outside cell
// mask
class InsideCellMask {
 public:
  InsideCellMask() = default;
  explicit InsideCellMask(const std::size_t cell_count)
      : cell_count_m(cell_count), words_m((cell_count + 63) / 64, 0) {}

  // setting the bit corresponding to the cell index to 1 (inside)
  void set(const std::size_t index) {
    if (index >= cell_count_m) throw std::out_of_range("InsideCellMask index");
    words_m[index / 64] |= std::uint64_t{1} << (index % 64);
  }

  // getting the bit corresponding to the cell index (1 = inside, 0 = outside)
  bool get(const std::size_t index) const {
    if (index >= cell_count_m) throw std::out_of_range("InsideCellMask index");
    return ((words_m[index / 64] >> (index % 64)) & std::uint64_t{1}) != 0;
  }

  // number of inside cells (i.e., number of bits set to 1)
  std::size_t cellCount() const { return cell_count_m; }

  // access to bit storage
  const std::vector<std::uint64_t>& words() const { return words_m; }
  std::vector<std::uint64_t>& words() { return words_m; }

 private:
  std::size_t cell_count_m = 0;
  std::vector<std::uint64_t> words_m;
};

// converting a 3D cell index (i,j,k) to a linear index (row major flattening)
inline std::size_t getLinearCellIndex(const int i, const int j, const int k,
                                      const int ny, const int nz) {
  return (static_cast<std::size_t>(i) * static_cast<std::size_t>(ny) +
          static_cast<std::size_t>(j)) *
             static_cast<std::size_t>(nz) +
         static_cast<std::size_t>(k);
}

inline void getCellIndicesFromLinearIndex(const std::size_t index, const int ny,
                                          const int nz, int* i, int* j,
                                          int* k) {
  const std::size_t plane = static_cast<std::size_t>(ny) * nz;
  *i = static_cast<int>(index / plane);
  const std::size_t remainder = index % plane;
  *j = static_cast<int>(remainder / static_cast<std::size_t>(nz));
  *k = static_cast<int>(remainder % static_cast<std::size_t>(nz));
}

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
struct SparseMixedCellMoments {
  static constexpr std::size_t NV =
      (VM_ORDER + 1) * (VM_ORDER + 2) * (VM_ORDER + 3) / 6;
  static constexpr std::size_t NS =
      (SM_ORDER + 1) * (SM_ORDER + 2) * (SM_ORDER + 3) / 6;

  std::uint32_t linear_index = 0;
  double volume[NV]{};
  double surface[NS]{};
};

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_2_SPARSE_MOMENTS_H_
