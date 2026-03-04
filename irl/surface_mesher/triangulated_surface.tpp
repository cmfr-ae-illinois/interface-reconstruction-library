// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SURFACE_MESH_TRIANGULATED_SURFACE_TPP_
#define IRL_SURFACE_MESH_TRIANGULATED_SURFACE_TPP_

namespace IRL {

inline TriangulatedSurfaceOutput::TriangulatedSurfaceOutput(
    TriangulatedSurfaceOutput&& a_rhs)
    : vertices_m(std::move(a_rhs.vertices_m)),
      bdy_edges_m(std::move(a_rhs.bdy_edges_m)),
      triangles_m(std::move(a_rhs.triangles_m)) {}

inline TriangulatedSurfaceOutput::TriangulatedSurfaceOutput(
    const TriangulatedSurfaceOutput& a_rhs)
    : vertices_m(a_rhs.vertices_m),
      bdy_edges_m(a_rhs.bdy_edges_m),
      triangles_m(a_rhs.triangles_m) {}

inline TriangulatedSurfaceOutput& TriangulatedSurfaceOutput::operator=(
    TriangulatedSurfaceOutput&& a_rhs) {
  if (this != &a_rhs) {
    vertices_m = std::move(a_rhs.vertices_m);
    bdy_edges_m = std::move(a_rhs.bdy_edges_m);
    triangles_m = std::move(a_rhs.triangles_m);
  }
  return *this;
}

inline TriangulatedSurfaceOutput& TriangulatedSurfaceOutput::operator=(
    const TriangulatedSurfaceOutput& a_rhs) {
  if (this != &a_rhs) {
    vertices_m = a_rhs.vertices_m;
    bdy_edges_m = a_rhs.bdy_edges_m;
    triangles_m = a_rhs.triangles_m;
  }
  return *this;
}

inline void TriangulatedSurfaceOutput::addVertex(const Pt& a_vertex) {
  vertices_m.push_back(a_vertex);
}

inline void TriangulatedSurfaceOutput::addBoundaryEdge(
    const UnsignedIndex_t a, const UnsignedIndex_t b) {
  bdy_edges_m.push_back(std::pair<UnsignedIndex_t, UnsignedIndex_t>({a, b}));
}

inline void TriangulatedSurfaceOutput::addTriangle(const UnsignedIndex_t a,
                                                   const UnsignedIndex_t b,
                                                   const UnsignedIndex_t c) {
  triangles_m.push_back(
      ProxyTri<PointStorage>::fromNoExistencePlane(vertices_m, {a, b, c}));
}

inline TriangulatedSurfaceOutput::PointStorage&
TriangulatedSurfaceOutput::getVertexList(void) {
  return vertices_m;
}

inline const TriangulatedSurfaceOutput::PointStorage&
TriangulatedSurfaceOutput::getVertexList(void) const {
  return vertices_m;
}

inline TriangulatedSurfaceOutput::EdgeStorage&
TriangulatedSurfaceOutput::getBoundaryEdgeList(void) {
  return bdy_edges_m;
}

inline const TriangulatedSurfaceOutput::EdgeStorage&
TriangulatedSurfaceOutput::getBoundaryEdgeList(void) const {
  return bdy_edges_m;
}

inline TriangulatedSurfaceOutput::TriangleStorage&
TriangulatedSurfaceOutput::getTriangleList(void) {
  return triangles_m;
}

inline const TriangulatedSurfaceOutput::TriangleStorage&
TriangulatedSurfaceOutput::getTriangleList(void) const {
  return triangles_m;
}

inline TriangulatedSurfaceOutput::PointStorage::size_type
TriangulatedSurfaceOutput::nVertices(void) const {
  return vertices_m.size();
}

inline TriangulatedSurfaceOutput::EdgeStorage::size_type
TriangulatedSurfaceOutput::nBoundaryEdges(void) const {
  return bdy_edges_m.size();
}

inline TriangulatedSurfaceOutput::TriangleStorage::size_type
TriangulatedSurfaceOutput::nTriangles(void) const {
  return triangles_m.size();
}

inline void TriangulatedSurfaceOutput::clearVertices(void) {
  vertices_m.clear();
}
inline void TriangulatedSurfaceOutput::clearBoundaryEdges(void) {
  bdy_edges_m.clear();
}
inline void TriangulatedSurfaceOutput::clearTriangles(void) {
  triangles_m.clear();
}
inline void TriangulatedSurfaceOutput::clearAll(void) {
  vertices_m.clear();
  bdy_edges_m.clear();
  triangles_m.clear();
}

inline void TriangulatedSurfaceOutput::refineSize(
    const double a_max_size, const UnsignedIndex_t a_compute_dim,
    std::function<double(const double a_x, const double a_y)> a_func) {
  const auto original_number_of_tris = this->nTriangles();
  UnsignedIndex_t new_pos = static_cast<UnsignedIndex_t>(this->nVertices());

  std::array<UnsignedIndex_t, 2> dims =
      a_compute_dim == 0   ? std::array<UnsignedIndex_t, 2>{{1, 2}}
      : a_compute_dim == 1 ? std::array<UnsignedIndex_t, 2>{{0, 2}}
                           : std::array<UnsignedIndex_t, 2>{{0, 1}};

  for (std::size_t i = 0; i < this->nTriangles(); ++i) {
    while (triangles_m[i].calculateAbsoluteVolume() > a_max_size) {
      //    if (triangle.calculateAbsoluteVolume() > a_max_size) {
      const auto triangle = triangles_m[i];
      const auto& indices = triangle.getIndexMapping();
      vertices_m.resize(new_pos + 3);

      vertices_m[new_pos] = 0.5 * (triangle[1] - triangle[0]) + triangle[0];
      vertices_m[new_pos + 1] = 0.5 * (triangle[2] - triangle[0]) + triangle[0];
      vertices_m[new_pos + 2] = 0.5 * (triangle[2] - triangle[1]) + triangle[1];

      for (UnsignedIndex_t n = 0; n < 3; ++n) {
        auto& vertex = vertices_m[new_pos + n];
        vertex[a_compute_dim] = a_func(vertex[dims[0]], vertex[dims[1]]);
      }

      triangles_m[i] = ProxyTri<PointStorage>::fromNoExistencePlane(
          vertices_m, {indices[0], new_pos, new_pos + 1});
      triangles_m.push_back(ProxyTri<PointStorage>::fromNoExistencePlane(
          vertices_m, {indices[1], new_pos + 2, new_pos}));
      triangles_m.push_back(ProxyTri<PointStorage>::fromNoExistencePlane(
          vertices_m, {new_pos, new_pos + 2, new_pos + 1}));
      triangles_m.push_back(ProxyTri<PointStorage>::fromNoExistencePlane(
          vertices_m, {indices[2], new_pos + 1, new_pos + 2}));

      new_pos += 3;
    }
  }
}

inline void TriangulatedSurfaceOutput::write(const std::string& filename) {
  // binary file
  std::string header_info = filename;
  char head[80];
  std::strncpy(head, header_info.c_str(), sizeof(header_info) - 1);
  char attribute[2] = "0";
  char dummy[4] = "0";
  const unsigned long nTriLong = triangles_m.size();
  std::ofstream myfile;
  myfile.open(filename + ".stl", std::ios::out | std::ios::binary);
  myfile.write(head, sizeof(head));
  myfile.write(reinterpret_cast<const char*>(&nTriLong), 4);

  // write down every triangle
  for (std::size_t i = 0; i < triangles_m.size(); ++i) {
    // normal vector coordinates
    myfile.write(dummy, 4);
    myfile.write(dummy, 4);
    myfile.write(dummy, 4);

    const auto& triangle = triangles_m[i];
    auto points = std::array<float, 9>(
        {static_cast<float>(triangle[0][0]), static_cast<float>(triangle[0][1]),
         static_cast<float>(triangle[0][2]), static_cast<float>(triangle[1][0]),
         static_cast<float>(triangle[1][1]), static_cast<float>(triangle[1][2]),
         static_cast<float>(triangle[2][0]), static_cast<float>(triangle[2][1]),
         static_cast<float>(triangle[2][2])});
    // Write all coordinates
    myfile.write(reinterpret_cast<char*>(points.data()), sizeof(float) * 9);
    myfile.write(attribute, 2);
  }
  myfile.close();
}

inline MixedPolygonBezierSurface::MixedPolygonBezierSurface(
    MixedPolygonBezierSurface&& a_rhs)
    : points_m(std::move(a_rhs.points_m)),
      polygons_m(std::move(a_rhs.polygons_m)),
      bezier_triangles_m(std::move(a_rhs.bezier_triangles_m)),
      boundary_m(std::move(a_rhs.boundary_m)) {}

inline MixedPolygonBezierSurface::MixedPolygonBezierSurface(
    const MixedPolygonBezierSurface& a_rhs)
    : points_m(a_rhs.points_m),
      polygons_m(a_rhs.polygons_m),
      bezier_triangles_m(a_rhs.bezier_triangles_m),
      boundary_m(a_rhs.boundary_m) {}

inline MixedPolygonBezierSurface& MixedPolygonBezierSurface::operator=(
    MixedPolygonBezierSurface&& a_rhs) {
  if (this != &a_rhs) {
    points_m = std::move(a_rhs.points_m);
    polygons_m = std::move(a_rhs.polygons_m);
    bezier_triangles_m = std::move(a_rhs.bezier_triangles_m);
    boundary_m = std::move(a_rhs.boundary_m);
  }
  return *this;
}

inline MixedPolygonBezierSurface& MixedPolygonBezierSurface::operator=(
    const MixedPolygonBezierSurface& a_rhs) {
  if (this != &a_rhs) {
    points_m = a_rhs.points_m;
    polygons_m = a_rhs.polygons_m;
    bezier_triangles_m = a_rhs.bezier_triangles_m;
    boundary_m = a_rhs.boundary_m;
  }
  return *this;
}

inline void MixedPolygonBezierSurface::addSurface(
    const MixedPolygonBezierSurface& a_rhs) {
  if (this != &a_rhs) {
    const UnsignedIndex_t npoints = points_m.size();
    const UnsignedIndex_t npolygons = polygons_m.second.size();
    const UnsignedIndex_t nboundaries = boundary_m.second.size();
    const UnsignedIndex_t nbeziertriangles = bezier_triangles_m.size();
    points_m.insert(points_m.end(), a_rhs.points_m.begin(),
                    a_rhs.points_m.end());
    polygons_m.first.insert(polygons_m.first.end(),
                            a_rhs.polygons_m.first.begin(),
                            a_rhs.polygons_m.first.end());
    polygons_m.second.insert(polygons_m.second.end(),
                             a_rhs.polygons_m.second.begin(),
                             a_rhs.polygons_m.second.end());
    bezier_triangles_m.insert(bezier_triangles_m.end(),
                              a_rhs.bezier_triangles_m.begin(),
                              a_rhs.bezier_triangles_m.end());
    boundary_m.first.insert(boundary_m.first.end(),
                            a_rhs.boundary_m.first.begin(),
                            a_rhs.boundary_m.first.end());
    boundary_m.second.insert(boundary_m.second.end(),
                             a_rhs.boundary_m.second.begin(),
                             a_rhs.boundary_m.second.end());
    for (UnsignedIndex_t i = npolygons; i < polygons_m.second.size(); i++) {
      polygons_m.second[i] += npoints;
    }
    for (UnsignedIndex_t i = nboundaries; i < boundary_m.second.size(); i++) {
      boundary_m.second[i] += npoints;
    }
    for (UnsignedIndex_t i = nbeziertriangles; i < bezier_triangles_m.size();
         i++) {
      for (UnsignedIndex_t j = 0; j < 6; j++) {
        bezier_triangles_m[i][j] += npoints;
      }
    }
  }
}

inline void MixedPolygonBezierSurface::addPoint(const Pt& a_pt,
                                                const double& a_weight) {
  points_m.push_back(std::pair<Pt, double>(a_pt, a_weight));
}

inline void MixedPolygonBezierSurface::addPoints(
    const std::vector<Pt>& a_pts, const std::vector<double>& a_weights) {
  if (a_weights.size() == 0 || a_weights.size() != a_pts.size()) {
    for (UnsignedIndex_t i = 0; i < a_pts.size(); i++) {
      points_m.push_back(std::pair<Pt, double>(a_pts[i], 1.0));
    }
  } else {
    for (UnsignedIndex_t i = 0; i < a_pts.size(); i++) {
      points_m.push_back(std::pair<Pt, double>(a_pts[i], a_weights[i]));
    }
  }
}

inline void MixedPolygonBezierSurface::addPolygon(const Polygon& a_polygon) {
  if (a_polygon.getNumberOfVertices() > 0) {
    polygons_m.first.push_back(a_polygon.getNumberOfVertices());
    for (const auto& pt : a_polygon) {
      polygons_m.second.push_back(points_m.size());
      points_m.push_back(std::pair<Pt, double>(pt, 1.0));
    }
  }
}

inline void MixedPolygonBezierSurface::addPolygons(
    const std::vector<Polygon>& a_polygons) {
  for (UnsignedIndex_t i = 0; i < a_polygons.size(); i++) {
    if (a_polygons[i].getNumberOfVertices() > 0) {
      polygons_m.first.push_back(a_polygons[i].getNumberOfVertices());
      for (const auto& pt : a_polygons[i]) {
        polygons_m.second.push_back(points_m.size());
        points_m.push_back(std::pair<Pt, double>(pt, 1.0));
      }
    }
  }
}

inline void MixedPolygonBezierSurface::addCellScalar(
    const std::string& name, const std::vector<double>& values) {
  cell_scalar_data_m.emplace_back(name, values);
}

template <class ContainerPoints>
inline void MixedPolygonBezierSurface::addBezierTriangle(
    const ContainerPoints& a_points) {
  bezier_triangles_m.resize(bezier_triangles_m.size() + 1);
  const UnsignedIndex_t size_container = a_points.size();
  bezier_triangles_m.back().resize(size_container);
  for (UnsignedIndex_t i = 0; i < size_container; i++) {
    bezier_triangles_m.back()[i] = a_points[i];
  }
}

template <class ContainerPoints, class ContainerWeights>
inline void MixedPolygonBezierSurface::addBezierTriangle(
    const ContainerPoints& a_points, const ContainerWeights& a_weights) {
  bezier_triangles_m.resize(bezier_triangles_m.size() + 1);
  const UnsignedIndex_t size_container = a_points.size();
  bezier_triangles_m.back().resize(size_container);
  for (UnsignedIndex_t i = 0; i < size_container; i++) {
    bezier_triangles_m.back()[i] = a_points[i];
    points_m[a_points[i]].second = a_weights[i];
  }
}

template <class ContainerPoints>
inline void MixedPolygonBezierSurface::addBezierTriangles(
    const std::vector<ContainerPoints>& a_points) {
  const UnsignedIndex_t triangle_size = bezier_triangles_m.size();
  if (a_points.size() > 0) {
    const UnsignedIndex_t size_container = a_points.back().size();
    bezier_triangles_m.resize(triangle_size + a_points.size());
    for (UnsignedIndex_t i = 0; i < a_points.size(); i++) {
      bezier_triangles_m[triangle_size + i].resize(size_container);
      for (UnsignedIndex_t j = 0; j < size_container; j++) {
        bezier_triangles_m[triangle_size + i][j] = a_points[i][j];
      }
    }
  }
}

template <class ContainerPoints, class ContainerWeights>
inline void MixedPolygonBezierSurface::addBezierTriangles(
    const std::vector<ContainerPoints>& a_points,
    const std::vector<ContainerWeights>& a_weights) {
  const UnsignedIndex_t triangle_size = bezier_triangles_m.size();
  if (a_points.size() > 0) {
    const UnsignedIndex_t size_container = a_points.back().size();
    bezier_triangles_m.resize(triangle_size + a_points.size());
    for (UnsignedIndex_t i = 0; i < a_points.size(); i++) {
      bezier_triangles_m[triangle_size + i].resize(size_container);
      for (UnsignedIndex_t j = 0; j < size_container; j++) {
        bezier_triangles_m[triangle_size + i][j] = a_points[i][j];
        points_m[a_points[i][j]].second = a_weights[i][j];
      }
    }
  }
}

template <class ContainerPoints>
inline void MixedPolygonBezierSurface::addBoundary(
    const ContainerPoints& a_points) {
  const UnsignedIndex_t size_container = a_points.size();
  if (size_container > 0) {
    boundary_m.first.push_back(size_container);
    for (UnsignedIndex_t j = 0; j < size_container; j++) {
      boundary_m.second.push_back(a_points[j]);
    }
  }
}

template <class ContainerPoints>
inline void MixedPolygonBezierSurface::addBoundaries(
    const std::vector<ContainerPoints>& a_points) {
  if (a_points.size() > 0) {
    for (UnsignedIndex_t i = 0; i < a_points.size(); i++) {
      const UnsignedIndex_t size_container = a_points[i].size();
      if (size_container > 0) {
        boundary_m.first.push_back(size_container);
        for (UnsignedIndex_t j = 0; j < size_container; j++) {
          boundary_m.second.push_back(a_points[i][j]);
        }
      }
    }
  }
}

inline MixedPolygonBezierSurface::PointStorage&
MixedPolygonBezierSurface::getPointList(void) {
  return points_m;
}

inline const MixedPolygonBezierSurface::PointStorage&
MixedPolygonBezierSurface::getPointList(void) const {
  return points_m;
}

inline MixedPolygonBezierSurface::PolygonStorage&
MixedPolygonBezierSurface::getPolygonList(void) {
  return polygons_m;
}

inline const MixedPolygonBezierSurface::PolygonStorage&
MixedPolygonBezierSurface::getPolygonList(void) const {
  return polygons_m;
}

inline MixedPolygonBezierSurface::BezierTriangleStorage&
MixedPolygonBezierSurface::getBezierTriangleList(void) {
  return bezier_triangles_m;
}

inline const MixedPolygonBezierSurface::BezierTriangleStorage&
MixedPolygonBezierSurface::getBezierTriangleList(void) const {
  return bezier_triangles_m;
}

inline MixedPolygonBezierSurface::BoundaryStorage&
MixedPolygonBezierSurface::getBoundaryList(void) {
  return boundary_m;
}

inline const MixedPolygonBezierSurface::BoundaryStorage&
MixedPolygonBezierSurface::getBoundaryList(void) const {
  return boundary_m;
}

inline UnsignedIndex_t MixedPolygonBezierSurface::nPoints(void) const {
  return points_m.size();
}

inline UnsignedIndex_t MixedPolygonBezierSurface::nPolygons(void) const {
  return polygons_m.first.size();
}

inline UnsignedIndex_t MixedPolygonBezierSurface::nBezierTriangles(void) const {
  return bezier_triangles_m.size();
}

inline UnsignedIndex_t MixedPolygonBezierSurface::nBoundaries(void) const {
  return boundary_m.first.size();
}

inline void MixedPolygonBezierSurface::clearPoints(void) { points_m.clear(); }

inline void MixedPolygonBezierSurface::clearPolygons(void) {
  polygons_m.first.clear();
  polygons_m.second.clear();
}

inline void MixedPolygonBezierSurface::clearBezierTriangles(void) {
  bezier_triangles_m.clear();
}

inline void MixedPolygonBezierSurface::clearBoundaries(void) {
  boundary_m.first.clear();
  boundary_m.second.clear();
}

inline void MixedPolygonBezierSurface::clearAll(void) {
  points_m.clear();
  polygons_m.first.clear();
  polygons_m.second.clear();
  bezier_triangles_m.clear();
  boundary_m.first.clear();
  boundary_m.second.clear();
}

inline void MixedPolygonBezierSurface::clearCellData(void) {
  cell_scalar_data_m.clear();
}

inline void MixedPolygonBezierSurface::write(const std::string& filename,
                                             const bool a_write_boundary) {
  FILE* file;
  file = fopen(std::string(filename + ".vtu").c_str(), "w");

  const UnsignedIndex_t number_of_points = points_m.size();
  const UnsignedIndex_t number_of_polygons = polygons_m.first.size();
  const UnsignedIndex_t number_of_triangles = bezier_triangles_m.size();
  const UnsignedIndex_t number_of_boundaries = boundary_m.first.size();

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file,
          "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
          "byte_order=\"LittleEndian\" header_type=\"UInt32\" "
          "compressor=\"vtkZLibDataCompressor\">\n");
  fprintf(file, "<UnstructuredGrid>\n");
  fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          number_of_points, number_of_triangles + number_of_polygons);
  fprintf(file,
          "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n");
  for (UnsignedIndex_t i = 0; i < number_of_points; ++i) {
    fprintf(file, "%15.8E %15.8E %15.8E ", std::get<0>(points_m[i])[0],
            std::get<0>(points_m[i])[1], std::get<0>(points_m[i])[2]);
  }
  fprintf(file, "\n</DataArray>\n</Points>\n");
  fprintf(file,
          "<PointData RationalWeights=\"RationalWeights\">\n<DataArray "
          "type=\"Float64\" Name=\"RationalWeights\" "
          "format=\"ascii\">\n");
  for (UnsignedIndex_t i = 0; i < number_of_points; ++i) {
    fprintf(file, "%15.8E ", std::get<1>(points_m[i]));
  }
  fprintf(file, "\n</DataArray>\n</PointData>\n");

  fprintf(file, "<Cells>\n");
  fprintf(
      file,
      "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (UnsignedIndex_t i = 0; i < number_of_triangles; ++i) {
    for (UnsignedIndex_t j = 0; j < bezier_triangles_m[i].size(); j++) {
      fprintf(file, "%d ", bezier_triangles_m[i][j]);
    }
  }
  for (UnsignedIndex_t i = 0; i < polygons_m.second.size(); ++i) {
    fprintf(file, "%d ", polygons_m.second[i]);
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file,
          "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  UnsignedIndex_t count = 0;
  for (UnsignedIndex_t i = 0; i < number_of_triangles; ++i) {
    count += bezier_triangles_m[i].size();
    fprintf(file, "%d ", count);
  }
  for (UnsignedIndex_t i = 0; i < number_of_polygons; ++i) {
    count += polygons_m.first[i];
    fprintf(file, "%d ", count);
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (UnsignedIndex_t i = 0; i < number_of_triangles; ++i) {
    fprintf(file, "76 ");
  }
  for (UnsignedIndex_t i = 0; i < number_of_polygons; ++i) {
    fprintf(file, "7 ");
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file, "</Cells>\n");

  // scalars per cell
  const UnsignedIndex_t ncells = number_of_triangles + number_of_polygons;

  if (!cell_scalar_data_m.empty()) {
    fprintf(file, "<CellData Scalars=\"%s\">\n",
            cell_scalar_data_m.front().first.c_str());

    for (const auto& kv : cell_scalar_data_m) {
      const std::string& name = kv.first;
      const std::vector<double>& vals = kv.second;

      fprintf(file,
              "<DataArray type=\"Float64\" Name=\"%s\" "
              "NumberOfComponents=\"1\" format=\"ascii\">\n",
              name.c_str());

      for (UnsignedIndex_t i = 0; i < ncells; ++i) {
        fprintf(file, "%15.8E ", vals[i]);
      }
      fprintf(file, "\n</DataArray>\n");
    }

    fprintf(file, "</CellData>\n");
  }

  fprintf(file, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
  fclose(file);

  if (number_of_boundaries > 0 && a_write_boundary) {
    file = fopen(std::string("boundary_" + filename + ".vtu").c_str(), "w");
    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file,
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\" header_type=\"UInt32\" "
            "compressor=\"vtkZLibDataCompressor\">\n");
    fprintf(file, "<UnstructuredGrid>\n");
    fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
            number_of_points, number_of_polygons + number_of_boundaries);
    fprintf(
        file,
        "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n");
    for (UnsignedIndex_t i = 0; i < number_of_points; ++i) {
      fprintf(file, "%15.8E %15.8E %15.8E ", std::get<0>(points_m[i])[0],
              std::get<0>(points_m[i])[1], std::get<0>(points_m[i])[2]);
    }
    fprintf(file, "\n</DataArray>\n</Points>\n");
    fprintf(file,
            "<PointData RationalWeights=\"RationalWeights\">\n<DataArray "
            "type=\"Float64\" Name=\"RationalWeights\" "
            "format=\"ascii\">\n");
    for (UnsignedIndex_t i = 0; i < number_of_points; ++i) {
      fprintf(file, "%15.8E ", std::get<1>(points_m[i]));
    }
    fprintf(file, "\n</DataArray>\n</PointData>\n");

    fprintf(file, "<Cells>\n");
    fprintf(
        file,
        "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    UnsignedIndex_t count = 0;
    for (UnsignedIndex_t i = 0; i < polygons_m.first.size(); ++i) {
      for (UnsignedIndex_t j = 0; j < polygons_m.first[i]; ++j) {
        fprintf(file, "%d ", polygons_m.second[count++]);
      }
      fprintf(file, "%d ", polygons_m.second[count - polygons_m.first[i]]);
    }
    for (UnsignedIndex_t i = 0; i < boundary_m.second.size(); ++i) {
      fprintf(file, "%d ", boundary_m.second[i]);
    }
    fprintf(file, "\n</DataArray>\n");

    fprintf(file,
            "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    count = 0;
    for (UnsignedIndex_t i = 0; i < number_of_polygons; ++i) {
      count += polygons_m.first[i] + 1;
      fprintf(file, "%d ", count);
    }
    for (UnsignedIndex_t i = 0; i < number_of_boundaries; ++i) {
      count += boundary_m.first[i];
      fprintf(file, "%d ", count);
    }
    fprintf(file, "\n</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (UnsignedIndex_t i = 0; i < number_of_polygons; ++i) {
      fprintf(file, "4 ");
    }
    for (UnsignedIndex_t i = 0; i < number_of_boundaries; ++i) {
      fprintf(file, "75 ");
    }
    fprintf(file, "\n</DataArray>\n");
    fprintf(file, "</Cells>\n");
    fprintf(file, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
    fclose(file);
  }
}

inline std::ostream& operator<<(
    std::ostream& out,
    const TriangulatedSurfaceOutput& a_triangulated_surface) {
  out << "triangulated surface has " << a_triangulated_surface.nTriangles()
      << " triangles\n";
  out << "a total of " << a_triangulated_surface.nVertices() << " vertices and "
      << a_triangulated_surface.nBoundaryEdges() << " boundaries\n";
  return out;
}
}  // namespace IRL

#endif  // IRL_SURFACE_MESH_TRIANGULATED_SURFACE_TPP_
