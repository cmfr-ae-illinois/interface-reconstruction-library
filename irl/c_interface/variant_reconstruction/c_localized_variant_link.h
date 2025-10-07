// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_VARIANT_LINK_H_
#define IRL_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_VARIANT_LINK_H_

#include "irl/c_interface/data_structures/c_object_allocation_server_localized_variant_link.h"
#include "irl/c_interface/planar_reconstruction/c_localizers.h"
#include "irl/c_interface/variant_reconstruction/c_separator_variant.h"
#include "irl/data_structures/object_allocation_server.h"
#include "irl/variant_reconstruction/separator_variant.h"

extern "C" {

struct c_LocVariantLink {
  IRL::LocalizedSeparatorVariantLink* obj_ptr;
  bool is_owning = false;
};

void c_LocVariantLink_new(c_LocVariantLink* a_self,
                          const c_PlanarLoc* a_localizer,
                          const c_SeparatorVariant* a_separator);

void c_LocVariantLink_newFromObjectAllocationServer(
    c_LocVariantLink* a_self,
    c_ObjServer_LocVariantLink* a_object_allocation_server,
    const c_PlanarLoc* a_localizer, const c_SeparatorVariant* a_separator);

void c_LocVariantLink_delete(c_LocVariantLink* a_self);

void c_LocVariantLink_setId(c_LocVariantLink* a_self, const int* a_id);

int c_LocVariantLink_getId(const c_LocVariantLink* a_self);

void c_LocVariantLink_setEdgeConnectivity(
    c_LocVariantLink* a_self, const int* a_plane_index,
    const c_LocVariantLink* a_ptr_to_neighbor);

void c_LocVariantLink_setEdgeConnectivityNull(c_LocVariantLink* a_self,
                                              const int* a_plane_index);

void c_LocVariantLink_printToScreen(const c_LocVariantLink* a_self);
}

#endif  // IRL_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_VARIANT_LINK_H_
