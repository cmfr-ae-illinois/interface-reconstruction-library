// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/generic_cutting/c_cut_polygon.h"

#include <cassert>

extern "C" {

void c_getPoly_RectCub_Poly_Sep(const c_RectCub* a_rectangular_cuboid,
                                const c_PlanarSep* a_separator,
                                const int* a_plane_index, c_Poly* a_polygon) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_polygon != nullptr);
  assert(a_polygon->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  (*a_polygon->obj_ptr) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
      *a_rectangular_cuboid->obj_ptr, *a_separator->obj_ptr,
      (*a_separator
            ->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index)]);
}

void c_getPoly_RectCub_Poly_Variant(const c_RectCub* a_rectangular_cuboid,
                                    const c_SeparatorVariant* a_separator,
                                    const int* a_plane_index,
                                    c_Poly* a_polygon) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_polygon != nullptr);
  assert(a_polygon->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_separator->obj_ptr)) {
    (*a_polygon->obj_ptr) =
        IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
            *a_rectangular_cuboid->obj_ptr, *separator,
            (*separator)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index)]);
  } else {
    throw std::runtime_error(
        "Polygon cannot be extrated from type different than PlanarSeparator");
  }
}

void c_getPoly_Tet_Poly(const c_Tet* a_tet, const c_PlanarSep* a_separator,
                        const int* a_plane_index, c_Poly* a_polygon) {
  assert(a_tet != nullptr);
  assert(a_tet->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_polygon != nullptr);
  assert(a_polygon->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  (*a_polygon->obj_ptr) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
      *a_tet->obj_ptr, *a_separator->obj_ptr,
      (*a_separator
            ->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index)]);
}

void c_getPoly_Hex_Poly(const c_Hex* a_hexahedron,
                        const c_PlanarSep* a_separator,
                        const int* a_plane_index, c_Poly* a_polygon) {
  assert(a_hexahedron != nullptr);
  assert(a_hexahedron->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_polygon != nullptr);
  assert(a_polygon->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  (*a_polygon->obj_ptr) = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
      *a_hexahedron->obj_ptr, *a_separator->obj_ptr,
      (*a_separator
            ->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index)]);
}

void c_getPoly_RectCub_DivPoly(const c_RectCub* a_rectangular_cuboid,
                               const c_PlanarSep* a_separator,
                               const int* a_plane_index,
                               c_DivPoly* a_divided_polygon) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_divided_polygon != nullptr);
  assert(a_divided_polygon->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  (*a_divided_polygon->obj_ptr) =
      IRL::getPlanePolygonFromReconstruction<IRL::DividedPolygon>(
          *a_rectangular_cuboid->obj_ptr, *a_separator->obj_ptr,
          (*a_separator
                ->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index)]);
  a_divided_polygon->obj_ptr->resetCentroid();
}

double c_getSA_RectCub_Sep(const c_RectCub* a_rectangular_cuboid,
                           const c_PlanarSep* a_separator) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  return IRL::getReconstructionSurfaceArea(*a_rectangular_cuboid->obj_ptr,
                                           *a_separator->obj_ptr);
}

double c_getSA_RectCub_Variant(const c_RectCub* a_rectangular_cuboid,
                               const c_SeparatorVariant* a_separator) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  if (IRL::PlanarSeparator* separator =
          std::get_if<IRL::PlanarSeparator>(a_separator->obj_ptr)) {
    return IRL::getReconstructionSurfaceArea(*a_rectangular_cuboid->obj_ptr,
                                             *separator);
  } else if (IRL::Paraboloid* paraboloid =
                 std::get_if<IRL::Paraboloid>(a_separator->obj_ptr)) {
    using VolumeAndSurface =
        IRL::AddSurfaceOutput<IRL::Volume,
                              IRL::ParaboloidParametrizedSurfaceOutput>;
    return IRL::getVolumeMoments<VolumeAndSurface>(
               *a_rectangular_cuboid->obj_ptr, *paraboloid)
        .getSurface()
        .getSurfaceArea();
  } else if (IRL::Cylinder* cylinder =
                 std::get_if<IRL::Cylinder>(a_separator->obj_ptr)) {
    using VolumeAndSurface =
        IRL::AddSurfaceOutput<IRL::Volume,
                              IRL::CylinderParametrizedSurfaceOutput>;
    return IRL::getVolumeMoments<VolumeAndSurface>(
               *a_rectangular_cuboid->obj_ptr, *cylinder)
        .getSurface()
        .getSurfaceArea();
  } else {
    throw std::runtime_error(
        "SeparatorVariant type not supported in c_getSA_RectCub_Variant");
  }
}

}  // end extern C
