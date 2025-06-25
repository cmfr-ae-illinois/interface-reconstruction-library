// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_DATA_STRUCTURES_C_OBJECT_ALLOCATION_SERVER_SEPARATOR_VARIANT_H_
#define IRL_C_INTERFACE_DATA_STRUCTURES_C_OBJECT_ALLOCATION_SERVER_SEPARATOR_VARIANT_H_

#include "irl/data_structures/object_allocation_server.h"
#include "irl/variant_reconstruction/separator_variant.h"

extern "C" {

struct c_ObjServer_SeparatorVariant {
  IRL::ObjectAllocationServer<IRL::SeparatorVariant>* obj_ptr = nullptr;
};

void c_ObjServer_SeparatorVariant_new(c_ObjServer_SeparatorVariant* a_self,
                                      const std::size_t* a_number_to_allocate);

void c_ObjServer_SeparatorVariant_delete(c_ObjServer_SeparatorVariant* a_self);
}

#endif  // IRL_C_INTERFACE_DATA_STRUCTURES_C_OBJECT_ALLOCATION_SERVER_SEPARATOR_VARIANT_H_
