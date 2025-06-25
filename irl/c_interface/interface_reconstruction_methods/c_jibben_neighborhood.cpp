// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/interface_reconstruction_methods/c_jibben_neighborhood.h"

#include <cassert>

extern "C" {

void c_JibbenNeigh_new(c_JibbenNeigh* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::JibbenNeighborhood;
}

void c_JibbenNeigh_delete(c_JibbenNeigh* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_JibbenNeigh_reserve(c_JibbenNeigh* a_self, const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->reserve(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_JibbenNeigh_localize(c_JibbenNeigh* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->localize();
}

void c_JibbenNeigh_setDelta(c_JibbenNeigh* a_self, const double* a_delta) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->setDelta(*a_delta);
}

void c_JibbenNeigh_addMember(c_JibbenNeigh* a_self, const c_Poly* a_polygon,
                             const double* a_weights) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_polygon != nullptr);
  assert(a_polygon->obj_ptr != nullptr);
  assert(a_weights != nullptr);
  a_self->obj_ptr->addMember(*a_polygon->obj_ptr, *a_weights);
}

void c_JibbenNeigh_emptyNeighborhood(c_JibbenNeigh* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->emptyNeighborhood();
}

void c_JibbenNeigh_setCenterOfStencil(c_JibbenNeigh* a_self,
                                      const int* a_center_cell_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_center_cell_index >= 0);
  assert(*a_center_cell_index < static_cast<int>(a_self->obj_ptr->size()));
  a_self->obj_ptr->setCenterOfStencil(
      static_cast<IRL::UnsignedIndex_t>(*a_center_cell_index));
}

}  // end extern "C"
