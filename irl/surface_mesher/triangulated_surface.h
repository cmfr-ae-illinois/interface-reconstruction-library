// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SURFACE_MESH_TRIANGULATED_SURFACE_H_
#define IRL_SURFACE_MESH_TRIANGULATED_SURFACE_H_

#include <fstream>
#include <functional>
#include <vector>

#include "irl/geometry/polygons/polygon.h"
#include "irl/geometry/polygons/tri.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

class TriangulatedSurfaceOutput {
 public:
  class PointStorage : public std::vector<Pt> {
   public:
    using pt_type = Pt;
  };
  using EdgeStorage = std::vector<std::pair<UnsignedIndex_t, UnsignedIndex_t>>;
  using TriangleStorage = std::vector<ProxyTri<PointStorage>>;

  /// \brief Default constructor.
  TriangulatedSurfaceOutput(void) = default;
  ~TriangulatedSurfaceOutput(void) = default;

  TriangulatedSurfaceOutput(const TriangulatedSurfaceOutput& a_rhs);
  TriangulatedSurfaceOutput(TriangulatedSurfaceOutput&& a_rhs);

  TriangulatedSurfaceOutput& operator=(const TriangulatedSurfaceOutput& a_rhs);
  TriangulatedSurfaceOutput& operator=(TriangulatedSurfaceOutput&& a_rhs);

  void addVertex(const Pt& a_vertex);
  void addBoundaryEdge(const UnsignedIndex_t a, const UnsignedIndex_t b);
  void addTriangle(const UnsignedIndex_t a, const UnsignedIndex_t b,
                   const UnsignedIndex_t c);

  PointStorage& getVertexList(void);
  const PointStorage& getVertexList(void) const;

  EdgeStorage& getBoundaryEdgeList(void);
  const EdgeStorage& getBoundaryEdgeList(void) const;

  TriangleStorage& getTriangleList(void);
  const TriangleStorage& getTriangleList(void) const;

  PointStorage::size_type nVertices(void) const;
  EdgeStorage::size_type nBoundaryEdges(void) const;
  TriangleStorage::size_type nTriangles(void) const;

  void clearVertices(void);
  void clearBoundaryEdges(void);
  void clearTriangles(void);
  void clearAll(void);
  void write(const std::string& filename);

  /// Refines triangles to be less than a_max_size. Only two of the dimensions
  /// are used, and the third (a_compute_dim) is comptued according to the
  /// provided function.
  void refineSize(
      const double a_max_size, const UnsignedIndex_t a_compute_dim,
      std::function<double(const double a_x, const double a_y)> a_func);

 private:
  // Includes both edge and hole vertices
  PointStorage vertices_m;
  EdgeStorage bdy_edges_m;
  TriangleStorage triangles_m;
};

class MixedPolygonBezierSurface {
 public:
  using PointStorage = std::vector<std::pair<Pt, double>>;
  using PolygonStorage =
      std::pair<std::vector<UnsignedIndex_t>, std::vector<UnsignedIndex_t>>;
  using BezierTriangle = std::vector<UnsignedIndex_t>;
  using BezierTriangleWeights = std::vector<double>;
  using BezierTriangleStorage = std::vector<BezierTriangle>;
  using BoundaryStorage =
      std::pair<std::vector<UnsignedIndex_t>, std::vector<UnsignedIndex_t>>;

  /// \brief Default constructor.
  MixedPolygonBezierSurface(void) = default;
  ~MixedPolygonBezierSurface(void) = default;

  MixedPolygonBezierSurface(const MixedPolygonBezierSurface& a_rhs);
  MixedPolygonBezierSurface(MixedPolygonBezierSurface&& a_rhs);

  MixedPolygonBezierSurface& operator=(const MixedPolygonBezierSurface& a_rhs);
  MixedPolygonBezierSurface& operator=(MixedPolygonBezierSurface&& a_rhs);

  void addSurface(const MixedPolygonBezierSurface& a_rhs);
  void addPoint(const Pt& a_pt, const double& weight = 1.0);
  void addPoints(const std::vector<Pt>& a_pts,
                 const std::vector<double>& a_weights = std::vector<double>(0));
  void addPolygon(const Polygon& a_polygon);
  void addPolygons(const std::vector<Polygon>& a_polygons);

  template <class ContainerPoints>
  void addBezierTriangle(const ContainerPoints& a_points);
  template <class ContainerPoints, class ContainerWeights>
  void addBezierTriangle(
      const ContainerPoints& a_points,
      const ContainerWeights& a_weights = ContainerWeights(1.0));

  template <class ContainerPoints>
  void addBezierTriangles(const std::vector<ContainerPoints>& a_points);
  template <class ContainerPoints, class ContainerWeights>
  void addBezierTriangles(const std::vector<ContainerPoints>& a_points,
                          const std::vector<ContainerWeights>& a_weights);

  template <class ContainerPoints>
  void addBoundary(const ContainerPoints& a_points);
  template <class ContainerPoints>
  void addBoundaries(const std::vector<ContainerPoints>& a_points);

  PointStorage& getPointList(void);
  const PointStorage& getPointList(void) const;

  PolygonStorage& getPolygonList(void);
  const PolygonStorage& getPolygonList(void) const;

  BezierTriangleStorage& getBezierTriangleList(void);
  const BezierTriangleStorage& getBezierTriangleList(void) const;

  BoundaryStorage& getBoundaryList(void);
  const BoundaryStorage& getBoundaryList(void) const;

  UnsignedIndex_t nPoints(void) const;
  UnsignedIndex_t nPolygons(void) const;
  UnsignedIndex_t nBezierTriangles(void) const;
  UnsignedIndex_t nBoundaries(void) const;

  void clearPoints(void);
  void clearPolygons(void);
  void clearBezierTriangles(void);
  void clearBoundaries(void);
  void clearAll(void);
  void write(const std::string& filename, const bool a_write_boundary = false);

 private:
  PointStorage points_m;
  PolygonStorage polygons_m;
  BezierTriangleStorage bezier_triangles_m;
  BoundaryStorage boundary_m;
};

inline std::ostream& operator<<(
    std::ostream& out, const TriangulatedSurfaceOutput& a_triangulated_surface);

}  // namespace IRL

#include "irl/surface_mesher/triangulated_surface.tpp"

#endif  // IRL_SURFACE_MESH_TRIANGULATED_SURFACE_H_
