// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_DATA_STRUCTURES_C_OBJECT_ALLOCATION_SERVER_BEZIER_TRIANGLES_H_
#define IRL_C_INTERFACE_DATA_STRUCTURES_C_OBJECT_ALLOCATION_SERVER_BEZIER_TRIANGLES_H_

#include "irl/data_structures/object_allocation_server.h"
#include "irl/surface_mesher/triangulated_surface.h"

extern "C" {

struct c_ObjServer_MixedPolygonBezierSurface {
  IRL::ObjectAllocationServer<IRL::MixedPolygonBezierSurface>* obj_ptr =
      nullptr;
};

void c_ObjServer_MixedPolygonBezierSurface_new(
    c_ObjServer_MixedPolygonBezierSurface* a_self,
    const std::size_t* a_number_to_allocate);

void c_ObjServer_MixedPolygonBezierSurface_delete(
    c_ObjServer_MixedPolygonBezierSurface* a_self);
}

#endif  // IRL_C_INTERFACE_DATA_STRUCTURES_C_OBJECT_ALLOCATION_SERVER_BEZIER_TRIANGLES_H_
