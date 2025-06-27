// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_JIBBEN_NEIGHBORHOOD_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_JIBBEN_NEIGHBORHOOD_H_

#include "irl/c_interface/geometry/polygons/c_polygon.h"
#include "irl/interface_reconstruction_methods/jibben_neighborhood.h"

extern "C" {

struct c_JibbenNeigh {
  IRL::JibbenNeighborhood* obj_ptr = nullptr;
};

void c_JibbenNeigh_new(c_JibbenNeigh* a_self);

void c_JibbenNeigh_delete(c_JibbenNeigh* a_self);

void c_JibbenNeigh_reserve(c_JibbenNeigh* a_self, const int* a_size);

void c_JibbenNeigh_setSize(c_JibbenNeigh* a_self, const int* a_size);

void c_JibbenNeigh_localize(c_JibbenNeigh* a_self);

void c_JibbenNeigh_setDelta(c_JibbenNeigh* a_self, const double* a_delta);

void c_JibbenNeigh_addMember(c_JibbenNeigh* a_self, const c_Poly* a_polygon,
                             const double* a_weights);

void c_JibbenNeigh_emptyNeighborhood(c_JibbenNeigh* a_self);

void c_JibbenNeigh_setCenterOfStencil(c_JibbenNeigh* a_self,
                                      const int* a_center_cell_index);
}

#endif  // IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_JIBBEN_NEIGHBORHOOD_H_
