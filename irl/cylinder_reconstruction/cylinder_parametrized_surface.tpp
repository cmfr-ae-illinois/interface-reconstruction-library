// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_CYLINDER_RECONSTRUCTION_CYLINDER_PARAMETRIZED_SURFACE_TPP_
#define IRL_CYLINDER_RECONSTRUCTION_CYLINDER_PARAMETRIZED_SURFACE_TPP_

#include <fstream>
#include <iomanip>

// #define IRL_USE_EARCUT
// #define IRL_USE_TRIANGLE
#define IRL_USE_CGAL
// #define IRL_USE_GEOGRAM

#include "external/NumericalIntegration/NumericalIntegration.h"
#include "irl/cylinder_reconstruction/cylinder_parametrized_surface.h"

// #define VALDEBUG2
#include "irl/helpers/mymath.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

namespace IRL {

template <class VertexList>
void projectOnSurface(VertexList& vertices, const AlignedCylinder cylinder,
                      const UnsignedIndex_t fixed_vertices) {
  if (vertices.size() > fixed_vertices) {
    using ScalarType = AlignedCylinder::value_type;
    for (UnsignedIndex_t i = fixed_vertices; i < vertices.size(); i++) {
      const double y = vertices[i][1];
      vertices[i][2] =
          sqrt(maximum(cylinder.r() - cylinder.b() * y * y, ScalarType(0)));
    }
  }
}

template <class VertexList, class EdgeList, class TriList>
void reMeshPolygon(VertexList& vertices, EdgeList& edges, TriList& triangles,
                   const double length_scale, const AlignedCylinder cylinder) {
  const double low = 4.0 / 5.0 * length_scale;
  const double high = 4.0 / 3.0 * length_scale;
  const int fixed_vertices = vertices.size();
  const int original_tris = triangles.size() / 3;

  // Construct edges and connectivity
  std::vector<std::array<int, 3>> tri_neigh;
  tri_neigh.resize(original_tris, std::array<int, 3>({-1, -1, -1}));
  std::vector<std::array<int, 3>> tri_edges;
  tri_edges.resize(original_tris, std::array<int, 3>({-1, -1, -1}));
  std::vector<std::vector<int>> vert_tri;
  vert_tri.resize(fixed_vertices);
  std::vector<std::vector<int>> vert_neigh;
  vert_neigh.resize(fixed_vertices);

  for (int i = 0; i < fixed_vertices; i++) {
    vert_tri[i].resize(0);
    vert_neigh[i].resize(0);
  }

  // List of triangles linked to a vertex
  for (int i = 0; i < original_tris; ++i) {
    for (int d = 0; d < 3; ++d) {
      vert_tri[triangles[3 * i + d]].push_back(i);
      tri_neigh[i][d] = -1;
      tri_edges[i][d] = -1;
    }
  }

  // Triangle neighbours
  for (int i = 0; i < original_tris; ++i) {
    // Loop over the 3 edges
    for (int d = 0; d < 3; ++d) {
      // tri_neigh[i][d] = -1;
      const int v0 = triangles[3 * i + d];
      const int v1 = triangles[3 * i + (d + 1) % 3];
      // Loop over triangles attached to start point
      bool found = false;
      for (int j = 0; j < vert_tri[v0].size(); ++j) {
        const int neigh = vert_tri[v0][j];
        if (neigh != i) {
          for (int k = 0; k < 3; ++k) {
            if (triangles[3 * neigh + k] == v1) {
              tri_neigh[i][d] = neigh;
              found = true;
              assert(triangles[3 * neigh + (k + 1) % 3] == v0 ||
                     triangles[3 * neigh + (k + 2) % 3] == v0);
              break;
            }
          }
        }
        if (found) {
          break;
        }
      }
      // if (!found) {
      //   std::cout << "Triangle " << i << " has boundary" << std::endl;
      // }
      assert(tri_neigh[i][d] != i);
    }
  }

  // Build edges
  for (int i = 0; i < original_tris; ++i) {
    // Loop over the 3 edges
    for (int d = 0; d < 3; ++d) {
      const int neigh = tri_neigh[i][d];
      if (i > neigh) {
        tri_edges[i][d] = edges.size();
        if (neigh >= 0) {
          bool found_edge = false;
          for (int k = 0; k < 3; k++) {
            if (tri_neigh[neigh][k] == i) {
              tri_edges[neigh][k] = edges.size();
              found_edge = true;
              break;
            }
          }
          // assert(found_edge);
        }
        edges.push_back({static_cast<int>(triangles[3 * i + d]),
                         static_cast<int>(triangles[3 * i + (d + 1) % 3]),
                         static_cast<int>(i), tri_neigh[i][d]});
      }
    }
  }

  bool edges_are_valid = true;
  for (int i = 0; i < original_tris; ++i) {
    for (int d = 0; d < 3; ++d) {
      if (tri_edges[i][d] < 0) {
        edges_are_valid = false;
      }
      // assert(tri_edges[i][d] >= 0);
    }
  }

  if (edges_are_valid) {
    // Build valence
    const int original_edges = edges.size();
    for (int i = 0; i < original_edges; ++i) {
      vert_neigh[edges[i][0]].push_back(edges[i][1]);
      vert_neigh[edges[i][1]].push_back(edges[i][0]);
    }

    bool no_duplicates = true;
    for (int i = 0; i < fixed_vertices; i++) {
      if (!noDuplicates(vert_neigh[i])) {
        no_duplicates = false;
        // std::cout << "Duplicate neighbours" << std::endl;
        break;
      }
      if (!noDuplicates(vert_tri[i])) {
        no_duplicates = false;
        // std::cout << "Duplicate triangles" << std::endl;
        break;
      }
    }

    if (no_duplicates) {
      for (int i = 0; i < 20; ++i) {
        splitLongEdges(vertices, triangles, edges, tri_neigh, tri_edges,
                       vert_tri, vert_neigh, high);
        // TODO: fix bug in edge collapse
        // collapseShortEdges(vertices, triangles, edges, tri_neigh, tri_edges,
        //                    vert_tri, vert_neigh, low, high, fixed_vertices);
        equalizeValence(vertices, triangles, edges, tri_neigh, tri_edges,
                        vert_tri, vert_neigh, fixed_vertices);
        relaxVertices(vertices, edges, vert_neigh, fixed_vertices, 30, 0.5);
        projectOnSurface(vertices, cylinder, fixed_vertices);
      }
    }
  }
}

inline CylinderParametrizedSurfaceOutput::CylinderParametrizedSurfaceOutput()
    : ParametrizedSurfaceOutput{},
      indexes_of_flip{},
      cylinder_m{},
      scale_m{} {}

inline CylinderParametrizedSurfaceOutput::CylinderParametrizedSurfaceOutput(
    const Cylinder& a_cylinder)
    : cylinder_m{{a_cylinder}},
      ParametrizedSurfaceOutput{},
      indexes_of_flip{},
      scale_m{} {}

inline CylinderParametrizedSurfaceOutput::CylinderParametrizedSurfaceOutput(
    CylinderParametrizedSurfaceOutput&& a_rhs)
    : ParametrizedSurfaceOutput(a_rhs), cylinder_m(a_rhs.cylinder_m),
      indexes_of_flip(a_rhs.indexes_of_flip),
      scale_m(a_rhs.scale_m) {}

inline CylinderParametrizedSurfaceOutput::CylinderParametrizedSurfaceOutput(
    const CylinderParametrizedSurfaceOutput& a_rhs)
    : ParametrizedSurfaceOutput(a_rhs), cylinder_m(a_rhs.cylinder_m),
      indexes_of_flip(a_rhs.indexes_of_flip),
      scale_m(a_rhs.scale_m) {}

inline CylinderParametrizedSurfaceOutput&
CylinderParametrizedSurfaceOutput::operator=(
    CylinderParametrizedSurfaceOutput&& a_rhs) {
  if (this != &a_rhs) {
    cylinder_m = a_rhs.cylinder_m;
    indexes_of_flip = a_rhs.indexes_of_flip;
    scale_m = a_rhs.scale_m;
    arc_list_m = std::move(a_rhs.arc_list_m);
    knows_surface_area_m = a_rhs.knows_surface_area_m;
    surface_area_m = a_rhs.surface_area_m;
    knows_avg_normal_m = a_rhs.knows_avg_normal_m;
    avg_normal_m = a_rhs.avg_normal_m;
    knows_int_mean_curv_m = a_rhs.knows_int_mean_curv_m;
    int_mean_curv_m = a_rhs.int_mean_curv_m;
    knows_int_gaussian_curv_m = a_rhs.knows_int_gaussian_curv_m;
    int_gaussian_curv_m = a_rhs.int_gaussian_curv_m;
    length_scale_m = a_rhs.length_scale_m;
    pt_from_bezier_split_m = std::move(a_rhs.pt_from_bezier_split_m);
  }
  return *this;
}

inline CylinderParametrizedSurfaceOutput&
CylinderParametrizedSurfaceOutput::operator=(
    const CylinderParametrizedSurfaceOutput& a_rhs) {
  if (this != &a_rhs) {
    cylinder_m = a_rhs.cylinder_m;
    indexes_of_flip = a_rhs.indexes_of_flip;
    scale_m = a_rhs.scale_m;
    arc_list_m = a_rhs.arc_list_m;
    knows_surface_area_m = a_rhs.knows_surface_area_m;
    surface_area_m = a_rhs.surface_area_m;
    knows_avg_normal_m = a_rhs.knows_avg_normal_m;
    avg_normal_m = a_rhs.avg_normal_m;
    knows_int_mean_curv_m = a_rhs.knows_int_mean_curv_m;
    int_mean_curv_m = a_rhs.int_mean_curv_m;
    knows_int_gaussian_curv_m = a_rhs.knows_int_gaussian_curv_m;
    int_gaussian_curv_m = a_rhs.int_gaussian_curv_m;
    length_scale_m = a_rhs.length_scale_m;
  }
  return *this;
}

inline void CylinderParametrizedSurfaceOutput::setCylinder(
    const Cylinder& a_cylinder) {
  indexes_of_flip.push_back(arc_list_m.size());
  cylinder_m.push_back(a_cylinder);
}

inline void CylinderParametrizedSurfaceOutput::setScale(
    double a_scale) {
  scale_m = a_scale;
}

inline void CylinderParametrizedSurfaceOutput::resetCylinder(void) {
  indexes_of_flip.resize(0);
  cylinder_m.resize(0);
}

// inline CylinderParametrizedSurfaceOutput::~CylinderParametrizedSurfaceOutput(
//     void) {
//   for (auto elem : pt_from_bezier_split_m) {
//     delete elem;
//   }
// }

inline double CylinderParametrizedSurfaceOutput::getSurfaceArea(void) {
  return 0.0;
};
inline double CylinderParametrizedSurfaceOutput::getMeanCurvatureIntegral(
    void) {
  return 0.0;
};
inline double CylinderParametrizedSurfaceOutput::getGaussianCurvatureIntegral(
    void) {
  return 0.0;
};
inline Normal CylinderParametrizedSurfaceOutput::getAverageNormal(void) {
  return Normal(0.0, 0.0, 0.0);
};
inline Normal CylinderParametrizedSurfaceOutput::getAverageNormalNonAligned(
    void) {
  return Normal(0.0, 0.0, 0.0);
};
inline Normal CylinderParametrizedSurfaceOutput::getNormalAligned(
    const Pt a_pt) {
  return Normal(0.0, 0.0, 0.0);
};
inline Normal CylinderParametrizedSurfaceOutput::getNormalNonAligned(
    const Pt a_pt) {
  return Normal(0.0, 0.0, 0.0);
};
inline double CylinderParametrizedSurfaceOutput::getMeanCurvatureAligned(
    const Pt a_pt) {
  return 0.0;
};
inline double CylinderParametrizedSurfaceOutput::getMeanCurvatureNonAligned(
    const Pt a_pt) {
  return 0.0;
};
inline double CylinderParametrizedSurfaceOutput::getGaussianCurvatureAligned(
    const Pt a_pt) {
  return 0.0;
};
inline double CylinderParametrizedSurfaceOutput::getGaussianCurvatureNonAligned(
    const Pt a_pt) {
  return 0.0;
};

inline TriangulatedSurfaceOutput CylinderParametrizedSurfaceOutput::triangulate(
    const double a_length_scale, const UnsignedIndex_t a_nsplit) const {
  TriangulatedSurfaceOutput returned_surface;
  this->triangulate_fromPtr(a_length_scale, a_nsplit, &returned_surface);
  return returned_surface;
}

inline void CylinderParametrizedSurfaceOutput::triangulate_fromPtr(
    const double a_length_scale, const UnsignedIndex_t a_nsplit,
    TriangulatedSurfaceOutput* returned_surface) const {
  const UnsignedIndex_t nArcs = this->size();
  double length_scale, length_scale_ref = length_scale_m;

  double inv_scale_sqr = double(1) / (scale_m * scale_m);

  auto an_indexes_of_flip = indexes_of_flip;
  int nb_rotation = an_indexes_of_flip.size();
  an_indexes_of_flip.push_back(nArcs);

  if (a_length_scale > 0.0) {
    length_scale_ref = a_length_scale;
  } else if (length_scale_ref <= 0.0) {
    auto surf = (*this);
    const double avg_length = std::sqrt(std::abs(surf.getSurfaceArea())) / 3.0;
    const double curv = std::fabs(surf.getAverageMeanCurvature());
    length_scale_ref = std::min(0.1 / curv, avg_length);
  }

#ifdef VALDEBUG2
  std::cout << "triangulating a cylinder surface\nThere is " << nArcs
            << " arcs\nLet's find the close curves\n";
  std::cout << "indexes of flip are : ";
  for (int indice : an_indexes_of_flip) {
    std::cout << indice << " ";
  }
  std::cout << std::endl;
#endif

  std::vector<std::vector<RationalBezierArc>> list_of_closed_curves;
  std::vector<bool> visited(nArcs, false);
  std::vector<int> indexes_of_close({0});

  // First, we need to order the arcs so as to form closed
  // curves
  bool valid_curves = true;
  for (int i = 0; i < nb_rotation; i++)
  {
    if (an_indexes_of_flip[i+1] == an_indexes_of_flip[i]) {
      indexes_of_close.push_back(list_of_closed_curves.size());
      continue;
    }
    for (std::size_t t = an_indexes_of_flip[i]; t < an_indexes_of_flip[i+1]; ++t) {
      if (visited[t]) {
        continue;
      }
      visited[t] = true;

      // Start with next available arc
      
      #ifdef VALDEBUG2
          std::cout << "starting a curve with an arc going from "
                    << arc_list_m[t].start_point() << " to "
                    << arc_list_m[t].end_point() << std::endl;
      #endif
      if (arc_list_m[t].arc_length() < DBL_EPSILON ) {
        #ifdef VALDEBUG2
            std::cout << "the arc is too small, it is not added" << std::endl;
        #endif
        list_of_closed_curves.push_back(
            std::vector<RationalBezierArc>());
      } else if (arc_list_m[t].weight() > 1.0e15) {
        const Pt p0 = arc_list_m[t].start_point();
        const Pt p1 = arc_list_m[t].control_point();
        const Pt p2 = arc_list_m[t].end_point();
        list_of_closed_curves.push_back(std::vector<RationalBezierArc>(
            {RationalBezierArc(p0, 0.5 * (p0 + p1), p1, 0.0),
            RationalBezierArc(p1, 0.5 * (p1 + p2), p2, 0.0)}));
      } else {
        list_of_closed_curves.push_back(
            std::vector<RationalBezierArc>({arc_list_m[t]}));
      }
      const std::uintptr_t start_id = arc_list_m[t].start_point_id();
      std::uintptr_t end_id = arc_list_m[t].end_point_id();
      Pt pe = arc_list_m[t].end_point();
      int counter = 0;
      while (end_id != start_id) {
        bool next_edge_found = false;
        for (std::size_t e = t + 1; e < an_indexes_of_flip[i+1]; ++e) {
          if (!visited[e]) {
            if (arc_list_m[e].start_point_id() == end_id) {
              visited[e] = true;
              next_edge_found = true;
              #ifdef VALDEBUG2
                  std::cout << "next arc is going from "
                            << arc_list_m[e].start_point() << " to "
                            << arc_list_m[e].end_point() << std::endl;
              #endif
              if (arc_list_m[e].arc_length() < DBL_EPSILON ) {
                #ifdef VALDEBUG2
                    std::cout << "the arc is too small, skipping" << std::endl;
                #endif
              } else if (arc_list_m[e].weight() > 1.0e15) {
                const Pt p0 = arc_list_m[e].start_point();
                const Pt p1 = arc_list_m[e].control_point();
                const Pt p2 = arc_list_m[e].end_point();
                list_of_closed_curves.back().push_back(
                    RationalBezierArc(p0, 0.5 * (p0 + p1), p1, 0.0));
                list_of_closed_curves.back().push_back(
                    RationalBezierArc(p1, 0.5 * (p1 + p2), p2, 0.0));
              } else {
                list_of_closed_curves.back().push_back(arc_list_m[e]);
              }
              end_id = arc_list_m[e].end_point_id();
              pe = arc_list_m[e].end_point();
              break;
            }
          }
        }
        if (!next_edge_found) {
          #ifdef VALDEBUG2
          std::cout << "could not find the next arc by comparing id. "
                        "trying to find by compaing distance" << std::endl;
          #endif
          double min_dist = DBL_MAX;
          std::size_t min_index = -1;
          for (std::size_t e = t + 1; e < an_indexes_of_flip[i+1]; ++e) {
            if (!visited[e]) {
              double dist = squaredMagnitude(pe - arc_list_m[e].start_point());
              #ifdef VALDEBUG2
              std::cout << "points : " << arc_list_m[e].start_point() << ", dist : " << dist << std::endl;
              #endif
              if (dist < min_dist) {
                min_dist = dist;
                min_index = e;
              }
            }
          }
          if (min_dist <= DBL_EPSILON) {
            next_edge_found = true;
            visited[min_index] = true;
            #ifdef VALDEBUG2
                std::cout << "next arc is going from "
                          << arc_list_m[min_index].start_point() << " to "
                          << arc_list_m[min_index].end_point() << std::endl;
            #endif
            if (arc_list_m[min_index].arc_length() < DBL_EPSILON ) {
              #ifdef VALDEBUG2
                  std::cout << "the arc is too small, skipping" << std::endl;
              #endif
            } else if (arc_list_m[min_index].weight() > 1.0e15) {
              const Pt p0 = arc_list_m[min_index].start_point();
              const Pt p1 = arc_list_m[min_index].control_point();
              const Pt p2 = arc_list_m[min_index].end_point();
              list_of_closed_curves.back().push_back(
                  RationalBezierArc(p0, 0.5 * (p0 + p1), p1, 0.0));
              list_of_closed_curves.back().push_back(
                  RationalBezierArc(p1, 0.5 * (p1 + p2), p2, 0.0));
            } else {
              list_of_closed_curves.back().push_back(arc_list_m[min_index]);
            }
            end_id = arc_list_m[min_index].end_point_id();
            pe = arc_list_m[min_index].end_point();
          } else {
            #ifdef VALDEBUG2
            std::cout << "could not find the next arc by compaing distance" << std::endl;
            #endif
          }
        }
        if (!next_edge_found) {
          #ifdef VALDEBUG2
          std::cout << "could not find the next arc by any means\n" 
                      "current ending point is " << pe << "\n"
                      "start point to join is " << arc_list_m[t].start_point() << "\n"
                      "dist : " << squaredMagnitude(pe - arc_list_m[t].start_point()) << std::endl;
          #endif
          // quick hack to end if the start and end point are at the same position but doenst have the same id
          if (squaredMagnitude(pe - arc_list_m[t].start_point()) <= DBL_EPSILON) {
            end_id = start_id;
            #ifdef VALDEBUG2
            std::cout << "because the points it is closed but with defferent ids, force end" << std::endl;
            #endif
          } else {
            #ifdef VALDEBUG2
            std::cout << "invalid curve" << std::endl;
            #endif
            valid_curves = false;
            break;
          }
        }
      }
  #ifdef VALDEBUG2
      std::cout << "end of that close curve\n";
  #endif
    }
    indexes_of_close.push_back(list_of_closed_curves.size());
  }
#ifdef VALDEBUG2
  std::cout << "in the end, there is " << list_of_closed_curves.size()
            << " close curves\n";
#endif

  returned_surface->clearAll();

  if (valid_curves) {
#ifdef IRL_USE_EARCUT
    // The number type to use for tessellation
    using Coord = double;
    auto aligned_cylinder = cylinder_m[0].getAlignedCylinder();
    // The index type. Defaults to uint32_t, but you can also
    // pass uint16_t if you know that your data won't have
    // more than 65536 vertices.

    length_scale = DBL_MAX;
    const UnsignedIndex_t nCurves =
        static_cast<UnsignedIndex_t>(list_of_closed_curves.size());

    // Create array
    using Point = std::array<Coord, 2>;
    std::vector<std::vector<Point>> polygon;
    polygon.resize(nCurves);

#ifdef VALDEBUG2
    std::cout << "let's define out polygons, there are " << nCurves
              << " of them\n";
#endif

    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
#ifdef VALDEBUG2
      std::cout << "polygon number " << i << " has " << nLocalArcs << " arcs\n";
#endif
      // Loop over arcs of curve
      UnsignedIndex_t added_points = 0;
      double signed_area = 0.0;
      for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
        // Compute approximate arc length
        const RationalBezierArc& arc = list_of_closed_curves[i][j];
        const double arc_length = arc.arc_length();
        // Split arc
        UnsignedIndex_t nSplit = a_nsplit <= 0 ? 1 : a_nsplit;
        if (length_scale_ref > 0.0) {
          nSplit = static_cast<UnsignedIndex_t>(arc_length / length_scale_ref);
          nSplit = nSplit < a_nsplit ? a_nsplit : nSplit;
        }
        const double step = 1.0 / static_cast<double>(nSplit);
        length_scale = std::min(length_scale, step * arc_length);
        if (length_scale_ref > 0.0) length_scale = length_scale_ref;
        Pt previous_pt = arc.point(0.0);
        for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
          const double t = static_cast<double>(k) * step;
          const auto pt = arc.point(t);
#ifdef VALDEBUG2
          std::cout << "adding point " << pt << "\n";
#endif
          polygon[i].push_back({pt[0], pt[1]});
          signed_area +=
              0.5 * (previous_pt[0] * pt[1] - pt[0] * previous_pt[1]);
          previous_pt = pt;
        }
      }

#ifdef VALDEBUG2
      std::cout << "polygon number " << i << " has " << polygon[i].size()
                << " points\n";
#endif
    }

    if (a_length_scale > 0.0) {
      length_scale = a_length_scale;
    }
    // Run tessellation
    // Returns array of indices that refer to the vertices of
    // the input polygon. e.g: the index 6 would refer to {25,
    // 75} in this example. Three subsequent indices form a
    // triangle. Output triangles are clockwise.
    std::vector<int> indices = mapbox::earcut<int>(polygon);
#ifdef VALDEBUG2
    std::cout << "this is indices after earcut : \n";
    for (int i = 0; i < indices.size(); i++) {
      std::cout << indices[i];
      if (i != indices.size() - 1) {
        std::cout << ", ";
      } else {
        std::cout << "\n";
      }
    }
#endif

    auto& vlist = returned_surface->getVertexList();
    auto& tlist = returned_surface->getTriangleList();
    std::vector<std::array<int, 4>> elist;

    UnsignedIndex_t count = 0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      count += polygon[i].size();
    }
    vlist.resize(count);
#ifdef VALDEBUG2
    std::cout << "let's map all the points (" << count << ")\n";
#endif

    count = 0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      for (UnsignedIndex_t j = 0; j < polygon[i].size(); ++j) {
        double x = polygon[i][j][0];
        double y = polygon[i][j][1];
        double z = sqrt(maximum(
            aligned_cylinder.r() * inv_scale_sqr - aligned_cylinder.b() * y * y, double(0)));
#ifdef VALDEBUG2
        std::cout << "point " << Pt(x, y, z) << " added to vlist\n";
#endif
        vlist[count++] = Pt(x, y, z);
      }
    }

    reMeshPolygon(vlist, elist, indices, length_scale, aligned_cylinder);

    count = 0;
    for (UnsignedIndex_t i = 0; i < indices.size() / 3; ++i) {
      if (indices[3 * i + 0] >= 0 && indices[3 * i + 1] >= 0 &&
          indices[3 * i + 2] >= 0) {
        count++;
      }
    }
#ifdef VALDEBUG2
    std::cout << "after remesh, vlist has " << vlist.size()
              << " elements, and count has a value of : " << count << "\n";
#endif
    // assert(count == indices.size() / 3);
    tlist.resize(count, TriangulatedSurfaceOutput::TriangleStorage::value_type::
                            fromNoExistencePlane(vlist, {0, 0, 0}));
    count = 0;
    for (UnsignedIndex_t i = 0; i < indices.size() / 3; ++i) {
      if (indices[3 * i + 0] >= 0 && indices[3 * i + 1] >= 0 &&
          indices[3 * i + 2] >= 0) {
        assert(indices[3 * i + 0] != indices[3 * i + 1]);
        assert(indices[3 * i + 1] != indices[3 * i + 2]);
        assert(indices[3 * i + 2] != indices[3 * i + 0]);
        tlist[count++] = TriangulatedSurfaceOutput::TriangleStorage::
            value_type::fromNoExistencePlane(
                vlist, {static_cast<UnsignedIndex_t>(indices[3 * i]),
                        static_cast<UnsignedIndex_t>(indices[3 * i + 1]),
                        static_cast<UnsignedIndex_t>(indices[3 * i + 2])});
      }
    }

    for (UnsignedIndex_t i = 0; i < elist.size(); ++i) {
      if ((elist[i][0] >= 0 || elist[i][1] >= 0) &&
          (elist[i][2] < 0 || elist[i][3] < 0)) {
        returned_surface->addBoundaryEdge(elist[i][0], elist[i][1]);
      }
    }

    // Translate and rotate triangulated surface vertices
    const auto& datum = cylinder_m[0].getDatum();
    const auto& ref_frame = cylinder_m[0].getReferenceFrame();
    for (auto& vertex : vlist) {
      const Pt base_pt = vertex;
      vertex = Pt(0.0, 0.0, 0.0);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          vertex[n] += ref_frame[d][n] * base_pt[d];
        }
      }
      vertex += datum;
    }

#elif defined IRL_USE_CGAL
    // Didn't do that
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kexact;
    typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Exact_predicates_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
    typedef CDT::Vertex_handle Vertex_handle;
    typedef CDT::Point Point;
    typedef CGAL::Arr_segment_traits_2<Kexact> Traits_2;
    typedef Traits_2::Curve_2 SegmentExact;
    typedef Kexact::Point_2 PointExact;

    auto& vlist = returned_surface->getVertexList();
    auto& tlist = returned_surface->getTriangleList();

    int vertix_offset = 0;
    int triangle_offset = 0;

    for (int r = 0; r < nb_rotation; r++)
    {
      if (indexes_of_close[r] == indexes_of_close[r+1]) {
        continue;
      }
      CDT cdt;
      const auto& aligned_cylinder = cylinder_m[r].getAlignedCylinder();

      #ifdef VALDEBUG2
          std::cout << "CGAL triangulation for cylinder " << aligned_cylinder << std::endl;
      #endif

      // std::ofstream myfile;
      // myfile.open("triangulation_log.txt");
      // myfile << "Starting triangulating surface.\n";
      // myfile << std::setprecision(16) << std::scientific
      //        << "Paraboloid: " << aligned_paraboloid << "\n";

      // Create boundaries
      std::vector<Point> points;
      std::list<Point> list_of_seeds;
      int nCurves = indexes_of_close[r+1]-indexes_of_close[r];
      UnsignedIndex_t start_points = 0;
      double total_signed_area = 0.0;
      double xmin = DBL_MAX, xmax = -DBL_MAX;
      double ymin = DBL_MAX, ymax = -DBL_MAX;
      UnsignedIndex_t vertex_count = 0;
      bool previous_valid = false;
      for (UnsignedIndex_t i = indexes_of_close[r]; i < indexes_of_close[r+1]; ++i) {
        #ifdef VALDEBUG2
            std::cout << "triangulating close loop # " << i << std::endl;
        #endif
        points.resize(0);
        const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
        // Loop over arcs of curve
        UnsignedIndex_t added_points = 0;
        double signed_area = 0.0;
        for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
          // Compute approximate arc length
          const RationalBezierArc& arc = list_of_closed_curves[i][j];
          #ifdef VALDEBUG2
              std::cout << "sub arc # " << j << std::endl;
              std::cout << "arc : " << arc << std::endl;
          #endif
          const double arc_length = arc.arc_length();
          // myfile << std::setprecision(16) << std::scientific << "Curve "
          // <<
          // i
          //        << " has arc: " << arc << "\n";
          // Split arc
          UnsignedIndex_t nSplit = a_nsplit <= 0 ? 1 : a_nsplit;
          if (length_scale_ref > 0.0) {
            nSplit = static_cast<UnsignedIndex_t>(arc_length / length_scale_ref);
            nSplit = nSplit < a_nsplit ? a_nsplit : nSplit;
          }
          const double step = 1.0 / static_cast<double>(nSplit);
          length_scale = std::min(length_scale, step * arc_length);
          if (length_scale_ref > 0.0) length_scale = length_scale_ref;
          Pt previous_pt = arc.point(0.0);
          for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
            const double t = static_cast<double>(k) * step;
            const auto pt = arc.point(t);
            // myfile << std::setprecision(16) << std::scientific << "Adding
            // vertex "
            //        << vertex_count++ << " at " << pt[0] << ", " << pt[1]
            //        <<
            //        ".\n";
            #ifdef VALDEBUG2
              std::cout << "ajout du point " << pt << std::endl;
            #endif
            points.push_back(Point(pt[0], pt[1]));
            signed_area +=
                0.5 * (previous_pt[0] * pt[1] - pt[0] * previous_pt[1]);
            previous_pt = pt;
          }
        }

        /* Remove duplicates */
        UnsignedIndex_t id0 = 0;
        do {
          UnsignedIndex_t id1 = (id0 + 1) % points.size();
          if ((points[id1].x() - points[id0].x()) *
                      (points[id1].x() - points[id0].x()) +
                  (points[id1].y() - points[id0].y()) *
                      (points[id1].y() - points[id0].y()) <
              1.0e12 * DBL_EPSILON * DBL_EPSILON) {
            // myfile << std::setprecision(16) << std::scientific
            //        << "Removing duplicate " << id1 << " at " <<
            //        points[id1].x()
            //        << ", " << points[id1].y() << " too close to " << id0
            //        << " at "
            //        << points[id0].x() << ", " << points[id0].y() <<
            //        ".\n";
            points.erase(points.begin() + id1);
            continue;
          } else {
            id0 = id1;
          }
        } while (id0 != 0);

        /* Create constraints */
        if (points.size() >= 3 &&
            std::fabs(signed_area) >
                std::max(1.0e-4 * length_scale * length_scale, 1.0e-14)) {
          // Construct the input segments.
          std::vector<SegmentExact> segments;
          segments.resize(points.size());
          segments[0] = SegmentExact(PointExact(points[points.size() - 1].x(),
                                                points[points.size() - 1].y()),
                                    PointExact(points[0].x(), points[0].y()));
          for (UnsignedIndex_t j = 0; j < points.size() - 1; ++j) {
            segments[j + 1] =
                SegmentExact(PointExact(points[j].x(), points[j].y()),
                            PointExact(points[j + 1].x(), points[j + 1].y()));
          }

          if (!CGAL::do_curves_intersect(segments.begin(), segments.end())) {
            if (nCurves > 1 && signed_area < 0.0) {
              // Add hole
              const auto p1x = CGAL::to_double(points[0].x());
              const auto p1y = CGAL::to_double(points[0].y());
              const auto p2x = CGAL::to_double(points[1].x());
              const auto p2y = CGAL::to_double(points[1].y());
              std::array<double, 2> hole_location{
                  {0.5 * (p1x + p2x), 0.5 * (p1y + p2y)}};
              Normal shift_dir = Normal(p2y - p1y, p1x - p2x, 0.0);
              shift_dir.normalize();
              // myfile << std::setprecision(16) << std::scientific <<
              // "Adding hole "
              //        << hole_location[0] + (1.0e3 * DBL_EPSILON) *
              //        shift_dir[0]
              //        << ", "
              //        << hole_location[1] + (1.0e3 * DBL_EPSILON) *
              //        shift_dir[1]
              //        << ".\n";
              list_of_seeds.push_back(
                  Point(hole_location[0] + (1.0e3 * DBL_EPSILON) * shift_dir[0],
                        hole_location[1] + (1.0e3 * DBL_EPSILON) * shift_dir[1]));
            }

            // Create segments
            // myfile << "Adding constraint " << points.size() - 1 << " -- "
            // << 0
            //        << ".\n";
            cdt.insert_constraint(points[points.size() - 1], points[0]);

            for (UnsignedIndex_t j = 0; j < points.size() - 1; ++j) {
              // myfile << "Adding constraint " << j << " -- " << j + 1 <<
              // ".\n";
              cdt.insert_constraint(points[j], points[j + 1]);
            }
          }
          start_points += added_points;
          total_signed_area += 0.5 * signed_area;
        }
      }

      // myfile << "Surface has area " << total_signed_area << "\n";
      // myfile << "Mesh has " << cdt.number_of_vertices() << "
      // vertices.\n"; myfile << "Refining with length-scale " <<
      // length_scale << ".\n"; sleep(1.0e-4);
      CGAL::refine_Delaunay_mesh_2(cdt,
                                  CGAL::parameters::seeds(list_of_seeds)
                                      .criteria(Criteria(0.15, length_scale)));
      // , CGAL::parameters::seeds_are_in_domain(false));
      // myfile << "Mesh has " << cdt.number_of_vertices() << "
      // vertices.\n"; myfile << "Mesh has " << cdt.number_of_faces() << "
      // faces.\n"; CGAL::lloyd_optimize_mesh_2(cdt,
      //                             CGAL::parameters::max_iteration_number
      //                             = 20);
      UnsignedIndex_t count = 0;
      CDT::Finite_faces_iterator face;
      // myfile << "Counting faces.\n";
      for (face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
          face++) {
        if (face->is_in_domain()) {
          count++;
        }
      }

      #ifdef VALDEBUG2
          std::cout << "adding " << count << " triangles" << std::endl;
      #endif

      vlist.resize(vertix_offset + 3 * count);
      tlist.resize(triangle_offset + count, TriangulatedSurfaceOutput::TriangleStorage::value_type::
                              fromNoExistencePlane(vlist, {0, 0, 0}));
      count = 0;
      // myfile << "Adding faces and vertices.\n";
      for (face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
          face++) {
        if (face->is_in_domain()) {
          // myfile << "Adding face " << count << ".\n";
          tlist[triangle_offset + count] = TriangulatedSurfaceOutput::TriangleStorage::value_type::
              fromNoExistencePlane(vlist,
                                  {vertix_offset + 3 * count, vertix_offset + 3 * count + 1, vertix_offset + 3 * count + 2});
          for (UnsignedIndex_t d = 0; d < 3; d++) {
            const double x = CGAL::to_double(face->vertex(d)->point().x());
            const double y = CGAL::to_double(face->vertex(d)->point().y());
            const double z = std::sqrt(std::max(
                0.0, aligned_cylinder.r() * inv_scale_sqr - aligned_cylinder.b() * y * y));
            vlist[vertix_offset + 3 * count + d] = Pt(x, y, z);
            auto neigh = face->neighbor(d);
            if (!neigh->is_in_domain()) {
              returned_surface->addBoundaryEdge(3 * count + (d + 1) % 3,
                                                3 * count + (d + 2) % 3);
            }
          }
          count++;
        }
      }

      // Translate and rotate triangulated surface vertices
      const auto& datum = cylinder_m[r].getDatum();
      const auto& ref_frame = cylinder_m[r].getReferenceFrame();
      for (int i = vertix_offset; i < vlist.size(); i++) {
        auto& vertex = vlist[i];
        const Pt base_pt = vertex;
        vertex = Pt(0.0, 0.0, 0.0);
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          for (UnsignedIndex_t n = 0; n < 3; ++n) {
            vertex[n] += ref_frame[d][n] * base_pt[d];
          }
        }
        vertex += datum;
      }
      vertix_offset = vlist.size();
      triangle_offset = tlist.size();
    }

// myfile << "Finished triangulating surface.\n";
// myfile.close();
#elif defined IRL_USE_TRIANGLE
// Didn't do that
// // Second, we approximate the arc length of the arc, so as
// // to know how many times it needs to be split
// std::vector<REAL> input_points;
// std::vector<REAL> input_holes;
// std::vector<int> input_segments;
// const UnsignedIndex_t nCurves =
//     static_cast<UnsignedIndex_t>(list_of_closed_curves.size());
// // Loop over curves
// UnsignedIndex_t start_points = 0;
// double total_signed_area = 0.0;
// for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
//   const UnsignedIndex_t nLocalArcs =
//   list_of_closed_curves[i].size();
//   // Loop over arcs of curve
//   UnsignedIndex_t added_points = 0;
//   double signed_area = 0.0;
//   for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
//     // Compute approximate arc length
//     const RationalBezierArc& arc = list_of_closed_curves[i][j];
//     // const auto& sp = arc.start_point();
//     // const auto& ep = arc.start_point();
//     // signed_area += (sp[0] * ep[1] - ep[0] * sp[1]);
//     const double arc_length = arc.arc_length();

//     // Split arc
//     UnsignedIndex_t nSplit = a_nsplit <= 0 ? 1 : a_nsplit;
//     if (length_scale_ref > 0.0) {
//       nSplit = static_cast<UnsignedIndex_t>(arc_length /
//       length_scale_ref); nSplit = nSplit < a_nsplit ? a_nsplit :
//       nSplit;
//     }
//     const double step = 1.0 / static_cast<double>(nSplit);
//     length_scale = std::min(length_scale, step * arc_length);
//     if (length_scale_ref > 0.0) length_scale = length_scale_ref;
//     // added_points += nSplit;
//     // const auto start_ind = input_points.size();
//     // input_points.resize(start_ind + 2 * nSplit);
//     // auto loc = input_points.begin() + start_ind;
//     Pt previous_pt = arc.point(0.0);
//     for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
//       const double t = static_cast<double>(k) * step;
//       const auto pt = arc.point(t);
//       if (squaredMagnitude(pt - previous_pt) >
//           1.0e8 * DBL_EPSILON * DBL_EPSILON) {
//         input_points.push_back(pt[0]);
//         input_points.push_back(pt[1]);
//         previous_pt = pt;
//         added_points++;
//       }
//     }
//   }

//   if (added_points >= 3) {
//     signed_area += (input_points[start_points + 2 * added_points -
//     2] *
//                         input_points[start_points + 1] -
//                     input_points[start_points + 0] *
//                         input_points[start_points + 2 *
//                         added_points - 1]);
//     for (UnsignedIndex_t j = 0; j < added_points - 1; ++j) {
//       signed_area += (input_points[start_points + 2 * j + 0] *
//                           input_points[start_points + 2 * j + 3] -
//                       input_points[start_points + 2 * j + 2] *
//                           input_points[start_points + 2 * j + 1]);
//     }

//     if (nCurves > 1 && signed_area < 0.0) {
//       // Add hole
//       const auto p1x = input_points[start_points];
//       const auto p1y = input_points[start_points + 1];
//       const auto p2x = input_points[start_points + 2];
//       const auto p2y = input_points[start_points + 3];
//       std::array<double, 2> hole_location{
//           {0.5 * (p1x + p2x), 0.5 * (p1y + p2y)}};
//       Normal shift_dir = Normal(p2y - p1y, p1x - p2x, 0.0);
//       shift_dir.normalize();
//       const auto start_ind = input_holes.size();
//       input_holes.resize(start_ind + 2);
//       input_holes[start_ind] = 0.0;
//       // hole_location[0] - (500.0 * DBL_EPSILON) *
//       // shift_dir[0];
//       input_holes[start_ind + 1] = 0.0;
//       // hole_location[1] - (500.0 * DBL_EPSILON) *
//       // shift_dir[1];
//     }

//     // Create segments
//     const int seg_size = input_segments.size();
//     input_segments.resize(seg_size + 2 * (added_points));
//     auto seg_loc = input_segments.begin() + seg_size;
//     *(seg_loc++) = start_points + added_points - 1;
//     *(seg_loc++) = start_points;
//     for (UnsignedIndex_t j = start_points;
//          j < start_points + added_points - 1; ++j) {
//       *(seg_loc++) = j;
//       *(seg_loc++) = j + 1;
//     }
//     start_points += added_points;
//     total_signed_area += 0.5 * signed_area;
//   }
// }

// // Below section is for Triangle library
// if (input_points.size() > 0) {
//   // std::cout << " Total area = " << total_signed_area <<
//   // " compared to "
//   //           << 2.0 * length_scale * length_scale <<
//   //           std::endl;
//   if (std::fabs(total_signed_area) > length_scale * length_scale) {
//     // Calling triangulation library
//     struct triangulateio in = {0}, out = {0};
//     in.numberofpoints = input_points.size() / 2;
//     in.pointlist = input_points.data();

//     std::vector<int> pointmarkerlist(in.numberofpoints, 1);
//     in.pointmarkerlist = pointmarkerlist.data();

//     in.numberofsegments = input_segments.size() / 2;
//     in.segmentlist = input_segments.data();
//     std::vector<int> segmentmarkerlist(in.numberofsegments, 1);
//     in.segmentmarkerlist = segmentmarkerlist.data();

//     in.numberofholes = input_holes.size() / 2;
//     if (in.numberofholes > 0) {
//       in.holelist = input_holes.data();
//     }

//     char flags[50];
//     sprintf(flags, "pza%.15feiQ", 0.5 * length_scale *
//     length_scale);

//     // std::cout << "Calling triangle with flags " <<
//     // flags << " and with "
//     //           << in.numberofpoints << " points and " <<
//     //           in.numberofsegments
//     //           << " segments and " << in.numberofholes
//     //           << " holes and max area = "
//     //           << 0.5 * length_scale * length_scale <<
//     //           std::endl;

//     // for (UnsignedIndex_t i = 0; i < in.numberofpoints;
//     // ++i) {
//     //   const double x = in.pointlist[2 * i + 0];
//     //   const double y = in.pointlist[2 * i + 1];
//     //   std::cout << "Point " << i << " = (" << x << ", "
//     //   << y << ")"
//     //             << std::endl;
//     // }
//     // for (UnsignedIndex_t i = 0; i <
//     // in.numberofsegments; ++i) {
//     //   const int j = in.segmentlist[2 * i + 0];
//     //   const int k = in.segmentlist[2 * i + 1];
//     //   std::cout << "Segment " << i << " = (" << j << ",
//     //   " << k << ")"
//     //             << std::endl;
//     // }

//     try {
//       triangulate_from_lib(flags, &in, &out, (struct
//       triangulateio*)NULL);
//       // std::cout << "Triangle finished" << std::endl;

//     } catch (std::runtime_error& e) {
//       // std::cerr << e.what() << std::endl;
//       // free(in.pointlist);
//       free(in.pointattributelist);
//       // free(in.pointmarkerlist);
//       free(in.trianglelist);
//       free(in.triangleattributelist);
//       free(in.trianglearealist);
//       free(in.neighborlist);
//       // free(in.segmentlist);
//       // free(in.segmentmarkerlist);
//       // free(in.holelist);
//       free(in.regionlist);
//       free(in.edgelist);
//       free(in.edgemarkerlist);
//       free(in.normlist);
//       free(out.pointlist);
//       free(out.pointattributelist);
//       free(out.pointmarkerlist);
//       free(out.trianglelist);
//       free(out.triangleattributelist);
//       free(out.trianglearealist);
//       free(out.neighborlist);
//       free(out.segmentlist);
//       free(out.segmentmarkerlist);
//       free(out.regionlist);
//       free(out.edgelist);
//       free(out.edgemarkerlist);
//       free(out.normlist);
//     }

//     auto& vlist = returned_surface->getVertexList();
//     vlist.resize(out.numberofpoints);
//     for (UnsignedIndex_t i = 0; i < out.numberofpoints; ++i) {
//       const double x = out.pointlist[2 * i + 0];
//       const double y = out.pointlist[2 * i + 1];
//       const double z =
//           -aligned_paraboloid.a() * x * x - aligned_paraboloid.b()
//           * y * y;
//       vlist[i] = Pt(x, y, z);
//     }

//     // Translate and rotate triangulated surface vertices
//     const auto& datum = paraboloid_m.getDatum();
//     const auto& ref_frame = paraboloid_m.getReferenceFrame();
//     for (auto& vertex : vlist) {
//       const Pt base_pt = vertex;
//       vertex = Pt(0.0, 0.0, 0.0);
//       for (UnsignedIndex_t d = 0; d < 3; ++d) {
//         for (UnsignedIndex_t n = 0; n < 3; ++n) {
//           vertex[n] += ref_frame[d][n] * base_pt[d];
//         }
//       }
//       vertex += datum;
//     }

//     for (UnsignedIndex_t i = 0; i < out.numberofedges; ++i) {
//       if (out.edgemarkerlist[i] == 1) {
//         returned_surface->addBoundaryEdge(out.edgelist[2 * i],
//                                           out.edgelist[2 * i + 1]);
//       }
//     }

//     auto& tlist = returned_surface->getTriangleList();
//     tlist.resize(out.numberoftriangles,
//                  TriangulatedSurfaceOutput::TriangleStorage::value_type::
//                      fromNoExistencePlane(vlist, {0, 0, 0}));
//     for (UnsignedIndex_t i = 0; i < out.numberoftriangles; ++i) {
//       tlist[i] =
//       TriangulatedSurfaceOutput::TriangleStorage::value_type::
//           fromNoExistencePlane(
//               vlist,
//               {static_cast<UnsignedIndex_t>(out.trianglelist[3 * i
//               + 0]),
//                static_cast<UnsignedIndex_t>(out.trianglelist[3 * i
//                + 1]),
//                static_cast<UnsignedIndex_t>(out.trianglelist[3 * i
//                + 2])});
//     }

//     /* free all allocated arrays, including those
//      * allocated by Triangle.
//      */
//     free(out.pointlist);
//     free(out.pointattributelist);
//     free(out.pointmarkerlist);
//     free(out.trianglelist);
//     free(out.triangleattributelist);
//     free(out.trianglearealist);
//     free(out.neighborlist);
//     free(out.segmentlist);
//     free(out.segmentmarkerlist);
//     free(out.regionlist);
//     free(out.edgelist);
//     free(out.edgemarkerlist);
//     free(out.normlist);
//     // free(in.pointlist);
//     free(in.pointattributelist);
//     // free(in.pointmarkerlist);
//     free(in.trianglelist);
//     free(in.triangleattributelist);
//     free(in.trianglearealist);
//     free(in.neighborlist);
//     // free(in.segmentlist);
//     // free(in.segmentmarkerlist);
//     // free(in.holelist);
//     free(in.regionlist);
//     free(in.edgelist);
//     free(in.edgemarkerlist);
//     free(in.normlist);
//   } else {  // Triangulate by hand
//     auto& vlist = returned_surface->getVertexList();
//     vlist.resize(input_points.size() / 2);
//     for (UnsignedIndex_t i = 0; i < input_points.size() / 2; ++i) {
//       const double x = input_points[2 * i + 0];
//       const double y = input_points[2 * i + 1];
//       const double z =
//           -aligned_paraboloid.a() * x * x - aligned_paraboloid.b()
//           * y * y;
//       vlist[i] = Pt(x, y, z);
//     }

//     // Translate and rotate triangulated surface vertices
//     const auto& datum = paraboloid_m.getDatum();
//     const auto& ref_frame = paraboloid_m.getReferenceFrame();
//     for (auto& vertex : vlist) {
//       const Pt base_pt = vertex;
//       vertex = Pt(0.0, 0.0, 0.0);
//       for (UnsignedIndex_t d = 0; d < 3; ++d) {
//         for (UnsignedIndex_t n = 0; n < 3; ++n) {
//           vertex[n] += ref_frame[d][n] * base_pt[d];
//         }
//       }
//       vertex += datum;
//     }

//     returned_surface->addBoundaryEdge(input_points.size() / 2 - 1,
//     0); for (UnsignedIndex_t i = 0; i < input_points.size() / 2 -
//     1; ++i) {
//       returned_surface->addBoundaryEdge(i, i + 1);
//     }

//     auto& tlist = returned_surface->getTriangleList();
//     tlist.resize(input_points.size() / 2 - 2,
//                  TriangulatedSurfaceOutput::TriangleStorage::value_type::
//                      fromNoExistencePlane(vlist, {0, 0, 0}));
//     for (UnsignedIndex_t i = 0; i < input_points.size() / 2 - 2;
//     ++i) {
//       tlist[i] =
//       TriangulatedSurfaceOutput::TriangleStorage::value_type::
//           fromNoExistencePlane(vlist, {0, i + 1, i + 2});
//     }
//   }
// }
#elif defined IRL_USE_GEOGRAM
    // Didn't do that
    // GEO::initialize();
    // auto cdt = GEO::CDT2d();
    // cdt.clear();
    // cdt.create_enclosing_quad(GEO::vec2(-0.5, -1.5), GEO::vec2(1.5,
    // -0.5),
    //                           GEO::vec2(1.5, 1.5),
    //                           GEO::vec2(-0.5, 1.5));

    // // Loop over curves
    // const UnsignedIndex_t nCurves =
    //     static_cast<UnsignedIndex_t>(list_of_closed_curves.size());
    // UnsignedIndex_t start_points = 0;
    // std::vector<double> points;
    // for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
    //   UnsignedIndex_t added_points = 0;
    //   const UnsignedIndex_t nLocalArcs =
    //   list_of_closed_curves[i].size();
    //   // Loop over arcs of curve
    //   for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
    //     // Compute approximate arc length
    //     const RationalBezierArc& arc = list_of_closed_curves[i][j];
    //     const double arc_length = arc.arc_length();

    //     // Split arc
    //     UnsignedIndex_t nSplit = a_nsplit <= 0 ? 1 : a_nsplit;
    //     if (length_scale_ref > 0.0) {
    //       nSplit = static_cast<UnsignedIndex_t>(arc_length /
    //       length_scale_ref); nSplit = nSplit < a_nsplit ? a_nsplit :
    //       nSplit;
    //     }
    //     const double step = 1.0 / static_cast<double>(nSplit);
    //     length_scale = std::min(length_scale, step * arc_length);
    //     if (length_scale_ref > 0.0) length_scale = length_scale_ref;
    //     for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
    //       const double t = static_cast<double>(k) * step;
    //       const auto pt = arc.point(t);
    //       points.push_back(pt[0]);
    //       points.push_back(pt[1]);
    //       added_points++;
    //     }
    //   }
    //   std::cout << "Inserting " << added_points << " points " <<
    //   std::endl; cdt.insert(static_cast<GEO::index_t>(added_points),
    //   points.data());

    //   if (added_points >= 3) {
    //     std::cout << "Inserting constraint (" << start_points +
    //     added_points
    //     - 1
    //               << ", " << start_points << ")" << std::endl;
    //     cdt.insert_constraint(
    //         static_cast<GEO::index_t>(4 + start_points + added_points -
    //         1), static_cast<GEO::index_t>(4 + start_points));
    //     for (UnsignedIndex_t j = start_points;
    //          j < start_points + added_points - 1; ++j) {
    //       std::cout << "Inserting constraint (" << j << ", " << j + 1
    //       << ")"
    //                 << std::endl;
    //       cdt.insert_constraint(static_cast<GEO::index_t>(4 + j),
    //                             static_cast<GEO::index_t>(4 + j + 1));
    //     }
    //     start_points += added_points;
    //   }
    // }

    // cdt.remove_external_triangles();
    // cdt.refine(18.0 * M_PI / 180.0, 0.00025);
    // // cdt.remove_external_triangles();
    // // cdt.remove_holes();

    // auto& vlist = returned_surface->getVertexList();
    // vlist.resize(cdt.nv());
    // for (UnsignedIndex_t i = 0; i < cdt.nv(); ++i) {
    //   const double x = cdt.point(i).x;
    //   const double y = cdt.point(i).y;
    //   const double z =
    //       -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y
    //       * y;
    //   vlist[i] = Pt(x, y, z);
    // }

    // // Translate and rotate triangulated surface vertices
    // const auto& datum = paraboloid_m.getDatum();
    // const auto& ref_frame = paraboloid_m.getReferenceFrame();
    // for (auto& vertex : vlist) {
    //   const Pt base_pt = vertex;
    //   vertex = Pt(0, 0, 0);
    //   for (UnsignedIndex_t d = 0; d < 3; ++d) {
    //     for (UnsignedIndex_t n = 0; n < 3; ++n) {
    //       vertex[n] += ref_frame[d][n] * base_pt[d];
    //     }
    //   }
    //   vertex += datum;
    // }

    // auto& tlist = returned_surface->getTriangleList();
    // tlist.resize(cdt.nT(),
    //              TriangulatedSurfaceOutput::TriangleStorage::value_type::
    //                  fromNoExistencePlane(vlist, {0, 0, 0}));
    // for (UnsignedIndex_t i = 0; i < cdt.nT(); ++i) {
    //   tlist[i] =
    //   TriangulatedSurfaceOutput::TriangleStorage::value_type::
    //       fromNoExistencePlane(vlist,
    //                            {static_cast<UnsignedIndex_t>(cdt.Tv(i,
    //                            0)),
    //                             static_cast<UnsignedIndex_t>(cdt.Tv(i,
    //                             1)),
    //                             static_cast<UnsignedIndex_t>(cdt.Tv(i,
    //                             2))});
    // }

#endif
  }
}

inline std::ostream& operator<<(
    std::ostream& out,
    const CylinderParametrizedSurfaceOutput& a_parametrized_surface) {
  for (UnsignedIndex_t i = 0; i < a_parametrized_surface.size(); ++i) {
    out << a_parametrized_surface[i];
    if (i < a_parametrized_surface.size() - 1) out << std::endl;
  }
  return out;
}

}  // namespace IRL

#endif  // IRL_CYLINDER_RECONSTRUCTION_CYLINDER_PARAMETRIZED_SURFACE_TPP_
