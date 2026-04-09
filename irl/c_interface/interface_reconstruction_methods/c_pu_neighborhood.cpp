// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood.h"

#include <cassert>

extern "C" {

void c_PUNeigh_new(c_PUNeigh* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::PUNeighborhood;
}

void c_PUNeigh_delete(c_PUNeigh* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_PUNeigh_reserve(c_PUNeigh* a_self, const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->reserve(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_PUNeigh_setSize(c_PUNeigh* a_self, const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->resize(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_PUNeigh_addMember(c_PUNeigh* a_self,
                         const c_SeparatorVariant* a_separator,
                         const double* __restrict__ a_centroid,
                         const double* a_weight) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_weight != nullptr);
  assert(a_centroid != nullptr);
  IRL::Pt centroid = IRL::Pt::fromRawDoublePointer(a_centroid);
  a_self->obj_ptr->addMember(*a_separator->obj_ptr, centroid, *a_weight);
}

void c_PUNeigh_setMember(c_PUNeigh* a_self, const int* a_index,
                         const c_SeparatorVariant* a_separator,
                         const double* __restrict__ a_centroid,
                         const double* a_weight) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_weight != nullptr);
  assert(a_centroid != nullptr);
  assert(*a_index >= 0);
  IRL::Pt centroid = IRL::Pt::fromRawDoublePointer(a_centroid);
  a_self->obj_ptr->setMember(static_cast<IRL::UnsignedIndex_t>(*a_index),
                             *a_separator->obj_ptr, centroid, *a_weight);
}

void c_PUNeigh_emptyNeighborhood(c_PUNeigh* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->emptyNeighborhood();
}

void c_PUNeigh_setCenterOfStencil(c_PUNeigh* a_self, const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  a_self->obj_ptr->setCenterOfStencil(
      static_cast<IRL::UnsignedIndex_t>(*a_index));
}

// add set member function here

}  // end extern "C"
