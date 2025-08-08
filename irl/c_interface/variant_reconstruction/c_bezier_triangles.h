// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_VARIANT_RECONSTRUCTION_C_BEZIER_TRIANGLES_H_
#define IRL_C_INTERFACE_VARIANT_RECONSTRUCTION_C_BEZIER_TRIANGLES_H_

#include "irl/c_interface/data_structures/c_object_allocation_server_bezier_triangles.h"
#include "irl/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "irl/c_interface/variant_reconstruction/c_separator_variant.h"
#include "irl/data_structures/object_allocation_server.h"
#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/paraboloid_reconstruction/paraboloid_parametrized_surface.h"
#include "irl/surface_mesher/triangulated_surface.h"

extern "C" {

struct c_MixedPolygonBezierSurface {
  IRL::MixedPolygonBezierSurface* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_MixedPolygonBezierSurface_new(c_MixedPolygonBezierSurface* a_self);

void c_MixedPolygonBezierSurface_newFromObjectAllocationServer(
    c_MixedPolygonBezierSurface* a_self,
    c_ObjServer_MixedPolygonBezierSurface* a_object_allocation_server);

void c_MixedPolygonBezierSurface_delete(c_MixedPolygonBezierSurface* a_self);

void c_MixedPolygonBezierSurface_clear(c_MixedPolygonBezierSurface* a_self);

int c_MixedPolygonBezierSurface_getNumberOfPoints(
    const c_MixedPolygonBezierSurface* a_self);

int c_MixedPolygonBezierSurface_getNumberOfPolygons(
    const c_MixedPolygonBezierSurface* a_self);

int c_MixedPolygonBezierSurface_getNumberOfTriangles(
    const c_MixedPolygonBezierSurface* a_self);

int c_MixedPolygonBezierSurface_getPolygonSize(
    const c_MixedPolygonBezierSurface* a_self, int* a_poly_size);

void c_MixedPolygonBezierSurface_getPolygonConnectivity(
    c_MixedPolygonBezierSurface* a_self, int* a_poly_index, int* a_poly_conn);

void c_MixedPolygonBezierSurface_getQuadraticTriangleConnectivity(
    c_MixedPolygonBezierSurface* a_self, int* a_tri_index, int* a_tri_conn);

void c_MixedPolygonBezierSurface_getPt(c_MixedPolygonBezierSurface* a_self,
                                       int* a_pt_index, double* a_pt);

void c_MixedPolygonBezierSurface_getSurface_RectCub_Variant(
    const c_RectCub* a_rectangular_cuboid,
    const c_SeparatorVariant* a_separator,
    c_MixedPolygonBezierSurface* a_surface);
}  // end extern C

#endif  // IRL_C_INTERFACE_VARIANT_RECONSTRUCTION_C_BEZIER_TRIANGLES_H_
