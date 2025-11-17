// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/variant_reconstruction/c_bezier_triangles.h"

#include <iostream>

extern "C" {

void c_MixedPolygonBezierSurface_new(c_MixedPolygonBezierSurface* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr = new IRL::MixedPolygonBezierSurface;
}

void c_MixedPolygonBezierSurface_newFromObjectAllocationServer(
    c_MixedPolygonBezierSurface* a_self,
    c_ObjServer_MixedPolygonBezierSurface* a_object_allocation_server) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  assert(a_object_allocation_server != nullptr);
  assert(a_object_allocation_server->obj_ptr != nullptr);
  a_self->is_owning = false;
  a_self->obj_ptr = a_object_allocation_server->obj_ptr->getNewObject();
}

void c_MixedPolygonBezierSurface_delete(c_MixedPolygonBezierSurface* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_MixedPolygonBezierSurface_clear(c_MixedPolygonBezierSurface* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->clearAll();
}

int c_MixedPolygonBezierSurface_getNumberOfPoints(
    const c_MixedPolygonBezierSurface* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->nPoints());
}

int c_MixedPolygonBezierSurface_getNumberOfPolygons(
    const c_MixedPolygonBezierSurface* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->nPolygons());
}

int c_MixedPolygonBezierSurface_getNumberOfTriangles(
    const c_MixedPolygonBezierSurface* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->nBezierTriangles());
}

int c_MixedPolygonBezierSurface_getPolygonSize(
    const c_MixedPolygonBezierSurface* a_self, int* a_poly_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_poly_index >= 0 && *a_poly_index < a_self->obj_ptr->nPolygons());
  return static_cast<int>(
      a_self->obj_ptr->getPolygonList().first[*a_poly_index]);
}

void c_MixedPolygonBezierSurface_getPolygonConnectivity(
    c_MixedPolygonBezierSurface* a_self, int* a_poly_index, int* a_poly_conn) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_poly_index >= 0 && *a_poly_index < a_self->obj_ptr->nPolygons());
  assert(a_poly_conn != nullptr);
  const IRL::UnsignedIndex_t nvertices =
      a_self->obj_ptr->getPolygonList().first[*a_poly_index];
  for (IRL::UnsignedIndex_t i = 0; i < nvertices; i++) {
    a_poly_conn[i] =
        a_self->obj_ptr->getPolygonList().second[*a_poly_index + i];
  }
}

void c_MixedPolygonBezierSurface_getQuadraticTriangleConnectivity(
    c_MixedPolygonBezierSurface* a_self, int* a_tri_index, int* a_tri_conn) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_tri_index >= 0 &&
         *a_tri_index < a_self->obj_ptr->nBezierTriangles());
  assert(a_tri_conn != nullptr);
  for (IRL::UnsignedIndex_t i = 0; i < 6; i++) {
    a_tri_conn[i] = a_self->obj_ptr->getBezierTriangleList()[*a_tri_index][i];
  }
}

void c_MixedPolygonBezierSurface_getPt(c_MixedPolygonBezierSurface* a_self,
                                       int* a_pt_index, double* a_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_pt_index >= 0 && *a_pt_index < a_self->obj_ptr->nPoints());
  assert(a_pt != nullptr);
  a_pt[0] = a_self->obj_ptr->getPointList()[*a_pt_index].first[0];
  a_pt[1] = a_self->obj_ptr->getPointList()[*a_pt_index].first[1];
  a_pt[2] = a_self->obj_ptr->getPointList()[*a_pt_index].first[2];
  a_pt[3] = a_self->obj_ptr->getPointList()[*a_pt_index].second;
}

void c_MixedPolygonBezierSurface_getSurface_RectCub_Variant(
    const c_RectCub* a_rectangular_cuboid,
    const c_SeparatorVariant* a_separator,
    c_MixedPolygonBezierSurface* a_surface) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_surface != nullptr);
  assert(a_surface->obj_ptr != nullptr);
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_separator->obj_ptr)) {
    for (IRL::UnsignedIndex_t n = 0; n < separator->getNumberOfPlanes(); n++) {
      const IRL::Polygon& polygon =
          IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
              *a_rectangular_cuboid->obj_ptr, *separator, (*separator)[n]);
      if (polygon.getNumberOfVertices() > 0) {
        a_surface->obj_ptr->addPolygon(polygon);
      }
    }
  } else if (IRL::Paraboloid* paraboloid =
                 std::get_if<IRL::Paraboloid>(a_separator->obj_ptr)) {
    using VolumeAndParaboloid =
        IRL::AddSurfaceOutput<IRL::Volume,
                              IRL::ParaboloidParametrizedSurfaceOutput>;
    VolumeAndParaboloid volume_and_surface =
        IRL::getVolumeMoments<VolumeAndParaboloid>(
            *a_rectangular_cuboid->obj_ptr, *paraboloid);
    a_surface->obj_ptr->addSurface(
        volume_and_surface.getSurface().getQuadraticBezierTriangleApprox());
  } else if (IRL::Cylinder* cylinder =
                 std::get_if<IRL::Cylinder>(a_separator->obj_ptr)) {
  using VolumeAndCylinder =
        IRL::AddSurfaceOutput<IRL::Volume,
                              IRL::CylinderParametrizedSurfaceOutput>;
    VolumeAndCylinder volume_and_surface =
        IRL::getVolumeMoments<VolumeAndCylinder>(
            *a_rectangular_cuboid->obj_ptr, *cylinder);
    a_surface->obj_ptr->addSurface(
        volume_and_surface.getSurface().getQuadraticBezierTriangleApprox());                
  }else {
    throw std::runtime_error("Unknown SeparatorVariant type in getSurface");
  }
}

}  // end extern C
