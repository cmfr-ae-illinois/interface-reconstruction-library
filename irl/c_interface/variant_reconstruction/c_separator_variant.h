// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_VARIANT_RECONSTRUCTION_C_SEPARATOR_VARIANT_H_
#define IRL_C_INTERFACE_VARIANT_RECONSTRUCTION_C_SEPARATOR_VARIANT_H_

#include "irl/c_interface/data_structures/c_object_allocation_server_separator_variant.h"
#include "irl/data_structures/object_allocation_server.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/variant_reconstruction/separator_variant.h"

extern "C" {

struct c_SeparatorVariant {
  IRL::SeparatorVariant* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_SeparatorVariant_new(c_SeparatorVariant* a_self);

void c_SeparatorVariant_newFromObjectAllocationServer(
    c_SeparatorVariant* a_self,
    c_ObjServer_SeparatorVariant* a_object_allocation_server);

void c_SeparatorVariant_delete(c_SeparatorVariant* a_self);

void c_SeparatorVariant_setNumberOfPlanes(c_SeparatorVariant* a_self,
                                          const int* a_number_to_set);

void c_SeparatorVariant_setPlane(c_SeparatorVariant* a_self,
                                 const int* a_plane_index_to_set,
                                 const double* a_normal,
                                 const double* a_distance);

void c_SeparatorVariant_copy(
    c_SeparatorVariant* a_self,
    const c_SeparatorVariant* a_other_planar_separator);

int c_SeparatorVariant_getNumberOfPlanes(const c_SeparatorVariant* a_self);

void c_SeparatorVariant_getPlane(c_SeparatorVariant* a_self, const int* a_index,
                                 double* a_plane_listed);

bool c_SeparatorVariant_isFlipped(const c_SeparatorVariant* a_self);

void c_SeparatorVariant_printToScreen(const c_SeparatorVariant* a_self);

}  // end extern C

#endif  // IRL_C_INTERFACE_VARIANT_RECONSTRUCTION_C_SEPARATOR_VARIANT_H_
