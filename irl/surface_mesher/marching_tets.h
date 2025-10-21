// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SURFACE_MESH_MARCHING_TETS_H_
#define IRL_SURFACE_MESH_MARCHING_TETS_H_

#include <functional>
#include <vector>

#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/helpers/mymath.h"
#include "irl/surface_mesher/triangulated_surface.h"


namespace IRL {

/**
 * @class MarchingTets
 * @brief Surface extraction via the Marching Tetrahedra algorithm.
 *
 * Decomposes each cube in a structured 3D grid into six tetrahedra
 * and polygonizes each against an isosurface value (default = 0.0).
 * Produces a triangulated surface suitable for export (e.g., `.vtp` for ParaView).
 */
class MarchingTets {
 public:
  using ImplicitF = std::function<double(Pt)>;

  MarchingTets(void) = default;
  ~MarchingTets(void) = default;
  MarchingTets(const RectangularCuboid a_domain, const ImplicitF a_function);

  void setDomain(const RectangularCuboid a_domain);
  void setFunction(const ImplicitF a_function);
  TriangulatedSurfaceOutput triangulate(
      const UnsignedIndex_t a_subdivision = 1);

 private:
  ImplicitF function_m;
  RectangularCuboid domain_m;
};

}  // namespace IRL

#include "irl/surface_mesher/marching_tets.tpp"

#endif  // IRL_SURFACE_MESH_MARCHING_TETS_H_

