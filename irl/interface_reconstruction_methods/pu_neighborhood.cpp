// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/interface_reconstruction_methods/pu_neighborhood.h"

namespace IRL {

PUNeighborhood::PUNeighborhood(void) {}

void PUNeighborhood::addMember(const SeparatorVariant& a_separator,
                               const Pt& a_centroid, const double a_weight) {
  separators_m.push_back(a_separator);
  centroids_m.push_back(a_centroid);
  weights_m.push_back(a_weight);
}

void PUNeighborhood::setMember(const UnsignedIndex_t a_index,
                               const SeparatorVariant& a_separator,
                               const Pt& a_centroid, const double a_weight) {
  separators_m[a_index] = a_separator;
  centroids_m[a_index] = a_centroid;
  weights_m[a_index] = a_weight;
}

void PUNeighborhood::emptyNeighborhood(void) {
  separators_m.resize(0);
  centroids_m.resize(0);
  weights_m.resize(0);
}
void PUNeighborhood::resize(const UnsignedIndex_t a_size) {
  separators_m.resize(a_size);
  centroids_m.resize(a_size);
  weights_m.resize(a_size);
}

void PUNeighborhood::reserve(const UnsignedIndex_t a_size) {
  separators_m.reserve(a_size);
  centroids_m.reserve(a_size);
  weights_m.reserve(a_size);
}

void PUNeighborhood::setCenterOfStencil(const UnsignedIndex_t a_index) {
  center_cell_index_m = a_index;
}

UnsignedIndex_t PUNeighborhood::size(void) const {
  assert(separators_m.size() == weights_m.size());
  return static_cast<UnsignedIndex_t>(separators_m.size());
}

const SeparatorVariant& PUNeighborhood::getSeparator(
    const UnsignedIndex_t a_index) const {
  return separators_m[a_index];
}

const double& PUNeighborhood::getWeight(const UnsignedIndex_t a_index) const {
  return weights_m[a_index];
}

const Pt& PUNeighborhood::getCentroid(const UnsignedIndex_t a_index) const {
  return centroids_m[a_index];
}

const UnsignedIndex_t& PUNeighborhood::getCenterOfStencil(void) const {
  return center_cell_index_m;
}

}  // namespace IRL
