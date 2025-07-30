// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/variant_reconstruction/c_localized_variant_link.h"

#include <cassert>

extern "C" {

void c_LocVariantLink_new(c_LocVariantLink* a_self,
                          const c_PlanarLoc* a_localizer,
                          const c_SeparatorVariant* a_separator) {
  assert(a_self->obj_ptr == nullptr);
  assert(a_localizer != nullptr);
  assert(a_localizer->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr = new IRL::LocalizedSeparatorVariantLink(
      a_localizer->obj_ptr, a_separator->obj_ptr);
}

void c_LocVariantLink_newFromObjectAllocationServer(
    c_LocVariantLink* a_self,
    c_ObjServer_LocVariantLink* a_object_allocation_server,
    const c_PlanarLoc* a_localizer, const c_SeparatorVariant* a_separator) {
  assert(a_self->obj_ptr == nullptr);
  assert(a_object_allocation_server != nullptr);
  assert(a_object_allocation_server->obj_ptr != nullptr);
  assert(a_localizer != nullptr);
  assert(a_localizer->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  a_self->is_owning = false;
  a_self->obj_ptr = a_object_allocation_server->obj_ptr->getNewObject();
  *a_self->obj_ptr = IRL::LocalizedSeparatorVariantLink(a_localizer->obj_ptr,
                                                        a_separator->obj_ptr);
}

void c_LocVariantLink_delete(c_LocVariantLink* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_LocVariantLink_setId(c_LocVariantLink* a_self, const int* a_id) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_id >= 0);
  a_self->obj_ptr->setId(static_cast<IRL::UnsignedIndex_t>(*a_id));
}

int c_LocVariantLink_getId(const c_LocVariantLink* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->getId());
}

void c_LocVariantLink_setEdgeConnectivity(
    c_LocVariantLink* a_self, const int* a_plane_index,
    const c_LocVariantLink* a_ptr_to_neighbor) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_ptr_to_neighbor != nullptr);
  assert(a_ptr_to_neighbor->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  a_self->obj_ptr->setEdgeConnectivity(
      static_cast<IRL::UnsignedIndex_t>(*a_plane_index),
      a_ptr_to_neighbor->obj_ptr);
}

void c_LocVariantLink_setEdgeConnectivityNull(c_LocVariantLink* a_self,
                                              const int* a_plane_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  a_self->obj_ptr->setEdgeConnectivity(
      static_cast<IRL::UnsignedIndex_t>(*a_plane_index), nullptr);
}

void c_LocVariantLink_printToScreen(const c_LocVariantLink* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << *a_self->obj_ptr << std::endl;
}
}
