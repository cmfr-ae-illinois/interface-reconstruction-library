// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_NEIGHBORHOOD_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_NEIGHBORHOOD_H_

#include "irl/c_interface/variant_reconstruction/c_separator_variant.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"

extern "C" {

struct c_PUNeigh {
  IRL::PUNeighborhood* obj_ptr = nullptr;
};

void c_PUNeigh_new(c_PUNeigh* a_self);

void c_PUNeigh_delete(c_PUNeigh* a_self);

void c_PUNeigh_reserve(c_PUNeigh* a_self, const int* a_size);

void c_PUNeigh_setSize(c_PUNeigh* a_self, const int* a_size);

void c_PUNeigh_addMember(c_PUNeigh* a_self,
                         const c_SeparatorVariant* a_separator,
                         const double* __restrict__ a_centroid,
                         const double* a_weight);

void c_PUNeigh_setMember(c_PUNeigh* a_self, const int* a_index,
                         const c_SeparatorVariant* a_separator,
                         const double* __restrict__ a_centroid,
                         const double* a_weight);

void c_PUNeigh_emptyNeighborhood(c_PUNeigh* a_self);

void c_PUNeigh_setCenterOfStencil(c_PUNeigh* a_self, const int* a_index);
}

#endif  // IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_PU_NEIGHBORHOOD_H_
