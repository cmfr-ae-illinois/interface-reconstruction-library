// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_PARAMETRIZED_SURFACE_TPP_
#define IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_PARAMETRIZED_SURFACE_TPP_

#include <Eigen/Dense>
#include <fstream>
#include <iomanip>

#include "external/NumericalIntegration/NumericalIntegration.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.h"
#include "irl/generic_cutting/quadratic_intersection/quadratic_intersection.h"
#include "irl/paraboloid_reconstruction/paraboloid_parametrized_surface.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

namespace IRL {

inline Normal computeNormalizedTangentAtPoint(
    const AlignedParaboloid& a_paraboloid, const Normal& a_plane_normal,
    const Pt& a_pt) {
  Normal surface_normal = getParaboloidSurfaceNormal(a_paraboloid, a_pt);
  surface_normal.approximatelyNormalize();
  Normal tangent_at_pt = crossProduct(a_plane_normal, surface_normal);
  if (squaredMagnitude(tangent_at_pt) < DBL_EPSILON * DBL_EPSILON) {
    return Normal(0.0, 0.0, 0.0);
  }
  const double normal_correction = tangent_at_pt * a_plane_normal;
  tangent_at_pt = tangent_at_pt - normal_correction * a_plane_normal;
  tangent_at_pt.normalize();
  return tangent_at_pt;
}

inline RationalBezierArc computeVerticalRationalBezierArc(
    const AlignedParaboloid& a_paraboloid, const Pt& pt_0, const Pt& pt_1) {
  const double DISTANCE_EPSILON = 1.0e2 * DBL_EPSILON;
  const double ANGLE_EPSILON = 1.0e6 * DBL_EPSILON;

  // Calculate edge vector and its normalized version
  const Normal edge_vector = pt_1 - pt_0;
  Normal edge_vector_normalized = edge_vector;
  edge_vector_normalized.normalize();

  const Normal plane_normal =
      crossProduct(edge_vector_normalized, Normal(0.0, 0.0, 1.0));
  Normal tangent_0 =
      computeNormalizedTangentAtPoint(a_paraboloid, plane_normal, pt_0);
  Normal tangent_1 =
      computeNormalizedTangentAtPoint(a_paraboloid, plane_normal, pt_1);

  // Compute dot product between normalized edge and end-point tangents
  double tgt0_dot_edge = tangent_0 * edge_vector_normalized;
  double tgt1_dot_edge = tangent_1 * edge_vector_normalized;
  if (tgt0_dot_edge < 0.0) {
    tgt0_dot_edge = -tgt0_dot_edge;
    tangent_0 = -tangent_0;
  }
  if (tgt1_dot_edge > 0.0) {
    tgt1_dot_edge = -tgt1_dot_edge;
    tangent_1 = -tangent_1;
  }

  if (magnitude(tangent_0) < 0.9 || magnitude(tangent_0) < 0.9 ||
      (magnitude(edge_vector) < DISTANCE_EPSILON &&
       fabs(1.0 - tgt0_dot_edge) < ANGLE_EPSILON &&
       fabs(1.0 + tgt1_dot_edge) < ANGLE_EPSILON)) {
    return RationalBezierArc(pt_0, 0.5 * (pt_0 + pt_1), pt_1, 0.0);
  }
  return RationalBezierArc(pt_0, tangent_0, pt_1, tangent_1, plane_normal,
                           a_paraboloid);
}

template <class VertexList>
void projectOnSurface(VertexList& vertices, const AlignedParaboloid paraboloid,
                      const UnsignedIndex_t fixed_vertices) {
  if (vertices.size() > fixed_vertices) {
    for (UnsignedIndex_t i = fixed_vertices; i < vertices.size(); i++) {
      const double x = vertices[i][0];
      const double y = vertices[i][1];
      vertices[i][2] = -paraboloid.a() * x * x - paraboloid.b() * y * y;
    }
  }
}

template <class VertexList, class EdgeList, class TriList>
void reMeshPolygon(VertexList& vertices, EdgeList& edges, TriList& triangles,
                   const double length_scale,
                   const AlignedParaboloid paraboloid) {
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
        projectOnSurface(vertices, paraboloid, fixed_vertices);
      }
    }
  }
}

inline ParaboloidParametrizedSurfaceOutput::
    ParaboloidParametrizedSurfaceOutput()
    : ParametrizedSurfaceOutput{} {}

inline ParaboloidParametrizedSurfaceOutput::ParaboloidParametrizedSurfaceOutput(
    const Paraboloid& a_paraboloid)
    : paraboloid_m{a_paraboloid}, ParametrizedSurfaceOutput{} {}

inline ParaboloidParametrizedSurfaceOutput::ParaboloidParametrizedSurfaceOutput(
    ParaboloidParametrizedSurfaceOutput&& a_rhs)
    : ParametrizedSurfaceOutput(a_rhs), paraboloid_m(a_rhs.paraboloid_m) {}

inline ParaboloidParametrizedSurfaceOutput::ParaboloidParametrizedSurfaceOutput(
    const ParaboloidParametrizedSurfaceOutput& a_rhs)
    : ParametrizedSurfaceOutput(a_rhs), paraboloid_m(a_rhs.paraboloid_m) {}

inline ParaboloidParametrizedSurfaceOutput&
ParaboloidParametrizedSurfaceOutput::operator=(
    ParaboloidParametrizedSurfaceOutput&& a_rhs) {
  if (this != &a_rhs) {
    paraboloid_m = a_rhs.paraboloid_m;
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

inline ParaboloidParametrizedSurfaceOutput&
ParaboloidParametrizedSurfaceOutput::operator=(
    const ParaboloidParametrizedSurfaceOutput& a_rhs) {
  if (this != &a_rhs) {
    paraboloid_m = a_rhs.paraboloid_m;
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

inline void ParaboloidParametrizedSurfaceOutput::setParaboloid(
    const Paraboloid& a_paraboloid) {
  paraboloid_m = a_paraboloid;
}

inline const Paraboloid& ParaboloidParametrizedSurfaceOutput::getParaboloid(
    void) const {
  return paraboloid_m;
}

class ArcContributionToParaboloidSurfaceArea_Functor {
 public:
  ArcContributionToParaboloidSurfaceArea_Functor(
      const RationalBezierArc& a_arc, const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive0 =
            (2. * pt0[0] *
                 std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                           4. * (b * b) * (pt0[1] * pt0[1])) -
             (1. + 4. * (b * b) * (pt0[1] * pt0[1])) *
                 std::log(-2. * std::fabs(a) * pt0[0] +
                          std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                                    4. * (b * b) * (pt0[1] * pt0[1]))) /
                 std::fabs(a)) /
            4.;
        const double primitive1 =
            (2. * pt1[0] *
                 std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                           4. * (b * b) * (pt1[1] * pt1[1])) -
             (1. + 4. * (b * b) * (pt1[1] * pt1[1])) *
                 std::log(-2. * std::fabs(a) * pt1[0] +
                          std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                                    4. * (b * b) * (pt1[1] * pt1[1]))) /
                 std::fabs(a)) /
            4.;
        return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
      } else {
        const double primitive0 =
            -(2. * pt0[1] *
                  std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                            4. * (b * b) * (pt0[1] * pt0[1])) -
              (1. + 4. * (a * a) * (pt0[0] * pt0[0])) *
                  std::log(-2. * std::fabs(b) * pt0[1] +
                           std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                                     4. * (b * b) * (pt0[1] * pt0[1]))) /
                  std::fabs(b)) /
            (4.);
        const double primitive1 =
            -(2. * pt1[1] *
                  std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                            4. * (b * b) * (pt1[1] * pt1[1])) -
              (1. + 4. * (a * a) * (pt1[0] * pt1[0])) *
                  std::log(-2. * std::fabs(b) * pt1[1] +
                           std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                                     4. * (b * b) * (pt1[1] * pt1[1]))) /
                  std::fabs(b)) /
            (4.);
        return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
      }
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return pt[0] * der[1];
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive =
            (2. * pt[0] *
                 std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                           4. * (b * b) * (pt[1] * pt[1])) -
             (1. + 4. * (b * b) * (pt[1] * pt[1])) *
                 std::log(-2. * std::fabs(a) * pt[0] +
                          std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                                    4. * (b * b) * (pt[1] * pt[1]))) /
                 std::fabs(a)) /
            4.;
        return primitive * der[1];
      } else {
        const double primitive =
            -(2. * pt[1] *
                  std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                            4. * (b * b) * (pt[1] * pt[1])) -
              (1. + 4. * (a * a) * (pt[0] * pt[0])) *
                  std::log(-2. * std::fabs(b) * pt[1] +
                           std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                                     4. * (b * b) * (pt[1] * pt[1]))) /
                  std::fabs(b)) /
            (4.);
        return primitive * der[0];
      }
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToParaboloidNormalX_Functor {
 public:
  ArcContributionToParaboloidNormalX_Functor(
      const RationalBezierArc& a_arc, const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive0 = a * pt0[0] * pt0[0];
      const double primitive1 = a * pt1[0] * pt1[0];
      return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive = a * pt[0] * pt[0];
      return primitive * der[1];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToParaboloidNormalY_Functor {
 public:
  ArcContributionToParaboloidNormalY_Functor(
      const RationalBezierArc& a_arc, const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive0 = -b * pt0[1] * pt0[1];
      const double primitive1 = -b * pt1[1] * pt1[1];
      return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive = -b * pt[1] * pt[1];
      return primitive * der[0];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToParaboloidNormalZ_Functor {
 public:
  ArcContributionToParaboloidNormalZ_Functor(
      const RationalBezierArc& a_arc, const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double primitive0 = pt0[0];
      const double primitive1 = pt1[0];
      return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double primitive = pt[0];
      return primitive * der[1];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToParaboloidMeanCurvature_Functor {
 public:
  ArcContributionToParaboloidMeanCurvature_Functor(
      const RationalBezierArc& a_arc, const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive0 =
            2. * b * pt0[0] +
            ((a + 4. * a * (b * b) * (pt0[1] * pt0[1]) -
              4. * (b * b * b) * (pt0[1] * pt0[1])) *
             std::atan((2. * a * pt0[0]) /
                       std::sqrt(1. + 4. * (b * b) * (pt0[1] * pt0[1])))) /
                (a * std::sqrt(1. + 4. * (b * b) * (pt0[1] * pt0[1])));
        const double primitive1 =
            2. * b * pt1[0] +
            ((a + 4. * a * (b * b) * (pt1[1] * pt1[1]) -
              4. * (b * b * b) * (pt1[1] * pt1[1])) *
             std::atan((2. * a * pt1[0]) /
                       std::sqrt(1. + 4. * (b * b) * (pt1[1] * pt1[1])))) /
                (a * std::sqrt(1. + 4. * (b * b) * (pt1[1] * pt1[1])));
        return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
      } else {
        const double primitive0 =
            -2. * a * pt0[1] -
            ((b - 4. * (a * a * a) * (pt0[0] * pt0[0]) +
              4. * (a * a) * b * (pt0[0] * pt0[0])) *
             std::atan((2. * b * pt0[1]) /
                       std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0])))) /
                (b * std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0])));
        const double primitive1 =
            -2. * a * pt1[1] -
            ((b - 4. * (a * a * a) * (pt1[0] * pt1[0]) +
              4. * (a * a) * b * (pt1[0] * pt1[0])) *
             std::atan((2. * b * pt1[1]) /
                       std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0])))) /
                (b * std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0])));
        return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
      }
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive =
            2. * b * pt[0] +
            ((a + 4. * a * (b * b) * (pt[1] * pt[1]) -
              4. * (b * b * b) * (pt[1] * pt[1])) *
             std::atan((2. * a * pt[0]) /
                       std::sqrt(1. + 4. * (b * b) * (pt[1] * pt[1])))) /
                (a * std::sqrt(1. + 4. * (b * b) * (pt[1] * pt[1])));
        return primitive * der[1];
      } else {
        const double primitive =
            -2. * a * pt[1] -
            ((b - 4. * (a * a * a) * (pt[0] * pt[0]) +
              4. * (a * a) * b * (pt[0] * pt[0])) *
             std::atan((2. * b * pt[1]) /
                       std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0])))) /
                (b * std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0])));
        return primitive * der[0];
      }
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};  // namespace IRL

class ArcContributionToParaboloidGaussianCurvature_Functor {
 public:
  ArcContributionToParaboloidGaussianCurvature_Functor(
      const RationalBezierArc& a_arc, const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else {
        const double primitive0 =
            4.0 * a * b * pt0[0] /
            ((1.0 + 4.0 * b * b * pt0[1] * pt0[1]) *
             std::sqrt(1.0 + 4.0 * a * a * pt0[0] * pt0[0] +
                       4.0 * b * b * pt0[1] * pt0[1]));
        const double primitive1 =
            4.0 * a * b * pt1[0] /
            ((1.0 + 4.0 * b * b * pt1[1] * pt1[1]) *
             std::sqrt(1.0 + 4.0 * a * a * pt1[0] * pt1[0] +
                       4.0 * b * b * pt1[1] * pt1[1]));
        return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
      }
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else {
        const double primitive = 4.0 * a * b * pt[0] /
                                 ((1.0 + 4.0 * b * b * pt[1] * pt[1]) *
                                  std::sqrt(1.0 + 4.0 * a * a * pt[0] * pt[0] +
                                            4.0 * b * b * pt[1] * pt[1]));
        return primitive * der[1];
      }
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};  // namespace IRL

inline double ParaboloidParametrizedSurfaceOutput::getSurfaceArea(void) {
  if (!knows_surface_area_m) {
    const UnsignedIndex_t nArcs = this->size();
    surface_area_m = 0.0;
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToParaboloidSurfaceArea_Functor functor(
          arc_list_m[t], aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      surface_area_m += integrator.quadratureAdaptive(functor, 0.0, 1.0, epsabs,
                                                      epsrel, quadrature_rule);
    }
    knows_surface_area_m = true;
  }
  return surface_area_m;
}

inline Normal ParaboloidParametrizedSurfaceOutput::getAverageNormal(void) {
  if (!knows_avg_normal_m) {
    const UnsignedIndex_t nArcs = this->size();
    avg_normal_m = Normal();
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToParaboloidNormalX_Functor functorx(arc_list_m[t],
                                                          aligned_paraboloid);
      ArcContributionToParaboloidNormalY_Functor functory(arc_list_m[t],
                                                          aligned_paraboloid);
      ArcContributionToParaboloidNormalZ_Functor functorz(arc_list_m[t],
                                                          aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      avg_normal_m[0] += integrator.quadratureAdaptive(
          functorx, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
      avg_normal_m[1] += integrator.quadratureAdaptive(
          functory, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
      avg_normal_m[2] += integrator.quadratureAdaptive(
          functorz, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
    }
    avg_normal_m.normalize();
    knows_avg_normal_m = true;
  }
  return avg_normal_m;
}

inline Normal ParaboloidParametrizedSurfaceOutput::getAverageNormalNonAligned(
    void) {
  auto aligned_normal = this->getAverageNormal();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  auto normal = Normal();
  for (std::size_t d = 0; d < 3; ++d) {
    for (std::size_t n = 0; n < 3; ++n) {
      normal[n] += ref_frame[d][n] * aligned_normal[d];
    }
  }
  return normal;
}

inline double ParaboloidParametrizedSurfaceOutput::getMeanCurvatureIntegral(
    void) {
  if (!knows_int_mean_curv_m) {
    const UnsignedIndex_t nArcs = this->size();
    int_mean_curv_m = 0.0;
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToParaboloidMeanCurvature_Functor functor(
          arc_list_m[t], aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      int_mean_curv_m += integrator.quadratureAdaptive(
          functor, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
    }
    knows_int_mean_curv_m = true;
  }
  return int_mean_curv_m;
}

inline double ParaboloidParametrizedSurfaceOutput::getGaussianCurvatureIntegral(
    void) {
  if (!knows_int_gaussian_curv_m) {
    const UnsignedIndex_t nArcs = this->size();
    int_gaussian_curv_m = 0.0;
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToParaboloidGaussianCurvature_Functor functor(
          arc_list_m[t], aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      int_gaussian_curv_m += integrator.quadratureAdaptive(
          functor, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
    }
    knows_int_gaussian_curv_m = true;
  }
  return int_gaussian_curv_m;
}

inline Normal ParaboloidParametrizedSurfaceOutput::getNormalAligned(
    const Pt a_pt) {
  auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
  auto aligned_normal = getParaboloidSurfaceNormal(aligned_paraboloid, a_pt);
  aligned_normal.normalize();
  return aligned_normal;
}

inline Normal ParaboloidParametrizedSurfaceOutput::getNormalNonAligned(
    const Pt a_pt) {
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());
  const Pt original_pt = a_pt - datum;
  auto aligned_pt = a_pt;
  for (std::size_t n = 0; n < 3; ++n) {
    aligned_pt[n] = ref_frame[n] * original_pt;
  }
  auto aligned_normal = this->getNormalAligned(aligned_pt);
  auto normal = Normal();
  for (std::size_t d = 0; d < 3; ++d) {
    for (std::size_t n = 0; n < 3; ++n) {
      normal[n] += ref_frame[d][n] * aligned_normal[d];
    }
  }
  return normal;
}

inline double ParaboloidParametrizedSurfaceOutput::getMeanCurvatureAligned(
    const Pt a_pt) {
  auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
  return (2. * (aligned_paraboloid.a() + aligned_paraboloid.b() +
                4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
                    aligned_paraboloid.b() * (a_pt[0] * a_pt[0]) +
                4. * aligned_paraboloid.a() *
                    (aligned_paraboloid.b() * aligned_paraboloid.b()) *
                    (a_pt[1] * a_pt[1]))) /
         std::pow(1. +
                      4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
                          (a_pt[0] * a_pt[0]) +
                      4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
                          (a_pt[1] * a_pt[1]),
                  1.5);
}

inline double ParaboloidParametrizedSurfaceOutput::getMeanCurvatureNonAligned(
    const Pt a_pt) {
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());
  const Pt original_pt = a_pt - datum;
  auto aligned_pt = a_pt;
  for (std::size_t n = 0; n < 3; ++n) {
    aligned_pt[n] = ref_frame[n] * original_pt;
  }
  return this->getMeanCurvatureAligned(aligned_pt);
}

inline double ParaboloidParametrizedSurfaceOutput::getGaussianCurvatureAligned(
    const Pt a_pt) {
  auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
  return 4. * aligned_paraboloid.a() * aligned_paraboloid.b() /
         ((1. +
           4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
               (a_pt[0] * a_pt[0]) +
           4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
               (a_pt[1] * a_pt[1])) *
          (1. +
           4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
               (a_pt[0] * a_pt[0]) +
           4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
               (a_pt[1] * a_pt[1])));
}

inline double
ParaboloidParametrizedSurfaceOutput::getGaussianCurvatureNonAligned(
    const Pt a_pt) {
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());
  const Pt original_pt = a_pt - datum;
  auto aligned_pt = a_pt;
  for (std::size_t n = 0; n < 3; ++n) {
    aligned_pt[n] = ref_frame[n] * original_pt;
  }
  return this->getGaussianCurvatureAligned(aligned_pt);
}

inline double ParaboloidParametrizedSurfaceOutput::getIntegrator(
    const F a_F, const bool useAdaptive,
    const Eigen::Integrator<double, 2>::QuadratureRule quadratureRule,
    const int npts) {
  double result = 0.0;

  // paraboloid params
  const auto& datum = this->getParaboloid().getDatum();
  const auto& frame = this->getParaboloid().getReferenceFrame();
  const auto& a = this->getParaboloid().getAlignedParaboloid().a();
  const auto& b = this->getParaboloid().getAlignedParaboloid().b();

  // reference point for pseduo-triangles
  Pt x_ref(0., 0., 0.);
  for (const auto& arc : arc_list_m) {
    x_ref += arc.start_point();
  }
  x_ref /= static_cast<double>(arc_list_m.size());

  // intgerating f(x) over clipped paraboloid
  for (const auto& arc : arc_list_m) {
    //  points for bezier quadrilaterals
    Pt A = arc.start_point();
    Pt B = arc.point(0.5);
    Pt C = arc.end_point();
    Pt E = x_ref;
    Pt D = 0.5 * (E + C);
    Pt F = 0.5 * (A + E);
    Pt G = 1. / 3. * (A + C + E);

    // arcs for bezier quadrilaterals
    auto [arc_1, arc_2] = arc.split();
    RationalBezierArc arc_3(C, 0.5 * (C + D), D, 1.0);
    RationalBezierArc arc_4(D, 0.5 * (D + E), E, 1.0);
    RationalBezierArc arc_5(E, 0.5 * (E + F), F, 1.0);
    RationalBezierArc arc_6(F, 0.5 * (F + A), A, 1.0);
    RationalBezierArc arc_7(B, 0.5 * (B + G), G, 1.0);
    RationalBezierArc arc_8(D, 0.5 * (D + G), G, 1.0);
    RationalBezierArc arc_9(F, 0.5 * (F + G), G, 1.0);

    // bezier quadrilaterals
    std::vector<RationalBezierArc> quad_1 = {arc_1, arc_7, arc_9, -arc_6};
    std::vector<RationalBezierArc> quad_2 = {arc_2, arc_3, -arc_8, arc_7};
    std::vector<RationalBezierArc> quad_3 = {arc_4, arc_5, -arc_9, arc_8};
    std::vector<std::vector<RationalBezierArc>> quad_list = {quad_1, quad_2,
                                                             quad_3};

    for (const auto& quad : quad_list) {
      const double a_bound = 0.0;
      const double b_bound = 1.0;
      auto functor = [=](const double u, const double v) {
        CoonsPatch coons(quad[0], quad[1], quad[2], quad[3]);
        Pt coons_val = coons.evaluate(u, v);
        double det_J = coons.detJacobian(u, v);
        double dS_factor =
            std::sqrt(1. + 4. * a * a * coons_val[0] * coons_val[0] +
                      4. * b * b * coons_val[1] * coons_val[1]);
        return (a_F(coons_val) * dS_factor * det_J);
      };
      double integral = 0.0;
      if (useAdaptive) {
        const std::size_t max_nsubdivisions = 256;
        const double epsabs = 10.0 * DBL_EPSILON;
        const double epsrel = 10.0 * DBL_EPSILON;
        for (std::size_t i = 1; i <= max_nsubdivisions; i *= 4) {
          Eigen::Integrator<double, 2> integrator(i);
          integral = integrator.quadratureAdaptive(
              functor, a_bound, b_bound, epsabs, epsrel, quadratureRule);
        }
      } else {
        GaussLegendreIntegrator<double, 2> integrator(npts);
        integral =
            integrator.integrate(functor, a_bound, a_bound, b_bound, b_bound);
      }
      result += integral;
    }
  }
  return result;
}

template <std::size_t ORDER>
inline GeneralSurfaceMoments3D<ORDER>
ParaboloidParametrizedSurfaceOutput::getSurfaceMoments() {
  static_assert(ORDER >= 0 && ORDER <= 2,
                "ONLY ORDER = 0, 1, or 2 supported for paraboloids");
  GeneralSurfaceMoments3D<ORDER> moments;

  const auto& a = this->getParaboloid().getAlignedParaboloid().a();
  const auto& b = this->getParaboloid().getAlignedParaboloid().b();
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();

  auto z = [a, b](const Pt& p) { return -a * p[0] * p[0] - b * p[1] * p[1]; };

  const double M0 = this->getIntegrator([](const Pt&) { return 1.0; });
  moments[0] = M0;

  if constexpr (ORDER == 0) return moments;

  const double M1x = this->getIntegrator([](const Pt& p) { return p[0]; });
  const double M1y = this->getIntegrator([](const Pt& p) { return p[1]; });
  const double M1z = this->getIntegrator([&](const Pt& p) { return z(p); });
  moments[1] = M1x;
  moments[2] = M1y;
  moments[3] = M1z;

  if constexpr (ORDER == 1)
    return (moments.moveAndRotateMoments(datum, ref_frame));

  const double Mxx =
      this->getIntegrator([](const Pt& p) { return p[0] * p[0]; });
  const double Mxy =
      this->getIntegrator([](const Pt& p) { return p[0] * p[1]; });
  const double Mxz =
      this->getIntegrator([&](const Pt& p) { return p[0] * z(p); });
  const double Myy =
      this->getIntegrator([](const Pt& p) { return p[1] * p[1]; });
  const double Myz =
      this->getIntegrator([&](const Pt& p) { return p[1] * z(p); });
  const double Mzz =
      this->getIntegrator([&](const Pt& p) { return z(p) * z(p); });
  moments[4] = Mxx;
  moments[5] = Mxy;
  moments[6] = Mxz;
  moments[7] = Myy;
  moments[8] = Myz;
  moments[9] = Mzz;
  moments.moveAndRotateMoments(datum, ref_frame);

  return moments;
}

inline MixedPolygonBezierSurface
ParaboloidParametrizedSurfaceOutput::getQuadraticBezierTriangleApprox(void) {
  return std::move(this->getBezierTriangleApprox(2));
}

inline MixedPolygonBezierSurface
ParaboloidParametrizedSurfaceOutput::getCubicBezierTriangleApprox(void) {
  return std::move(this->getBezierTriangleApprox(3));
}

inline MixedPolygonBezierSurface
ParaboloidParametrizedSurfaceOutput::getBezierTriangleApprox(
    const UnsignedIndex_t a_order) {
  // Initialize bezier surface
  MixedPolygonBezierSurface bezier_surface;

  // First, let's generate list of closed curves
  const UnsignedIndex_t nArcs = arc_list_m.size();
  std::vector<std::vector<UnsignedIndex_t>> list_of_closed_curves(0);
  std::vector<bool> visited(nArcs, false);
  bool valid_curves = true;
  for (UnsignedIndex_t t = 0; t < nArcs; ++t) {
    if (visited[t]) {
      continue;
    }
    visited[t] = true;
    // Start with next available arc
    list_of_closed_curves.push_back(std::vector<UnsignedIndex_t>({t}));
    const std::uintptr_t start_id = arc_list_m[t].start_point_id();
    std::uintptr_t end_id = arc_list_m[t].end_point_id();
    UnsignedIndex_t counter = 0;
    while (end_id != start_id) {
      for (UnsignedIndex_t e = t + 1; e < nArcs; ++e) {
        if (arc_list_m[e].start_point_id() == end_id) {
          visited[e] = true;
          list_of_closed_curves.back().push_back(e);
          end_id = arc_list_m[e].end_point_id();
          break;
        }
      }
      if (++counter > nArcs) {
        valid_curves = false;
        break;
      }
    }
  }

  // Only if we have found closed curves, then produce bezier triangle
  // approximation
  if (valid_curves) {
    const UnsignedIndex_t nCurves = list_of_closed_curves.size();
    std::vector<std::vector<std::array<double, 2>>> polygon(nCurves);

    UnsignedIndex_t npoints = 0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      const int nLocalArcs = list_of_closed_curves[i].size();
      for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
        const UnsignedIndex_t arc_id = list_of_closed_curves[i][j];
        const Pt& pt = arc_list_m[arc_id].start_point();
        const UnsignedIndex_t next_id = npoints + (j + 1) % nLocalArcs;
        polygon[i].push_back({pt[0], pt[1]});
      }
      npoints += nLocalArcs;
    }

    std::vector<Pt> points(npoints);
    std::vector<double> weights(npoints);
    std::vector<std::tuple<UnsignedIndex_t, UnsignedIndex_t, UnsignedIndex_t>>
        info(npoints);
    UnsignedIndex_t count = 0;
    npoints = 0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      const int nLocalArcs = list_of_closed_curves[i].size();
      for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
        const UnsignedIndex_t arc_id = list_of_closed_curves[i][j];
        points[count] = arc_list_m[arc_id].start_point();
        weights[count] = 1.0;
        info[count++] = std::make_tuple(npoints + (j + 1) % nLocalArcs, i, j);
      }
      npoints += nLocalArcs;
    }

    // Compute earcut triangulation of region constrained by closed curves
    std::vector<UnsignedIndex_t> indices =
        mapbox::earcut<UnsignedIndex_t>(polygon);

    // Convert flat triangles into quadratic rational Bezier triangle
    const UnsignedIndex_t ntriangles = indices.size() / 3;
    if (a_order == 2) {
      std::vector<std::array<UnsignedIndex_t, 6>> bezier_triangles(ntriangles);
      std::vector<std::array<UnsignedIndex_t, 3>> boundaries(npoints);
      const auto& aligned_p = paraboloid_m.getAlignedParaboloid();
      for (UnsignedIndex_t i = 0; i < indices.size() / 3; ++i) {
        for (int j = 0; j < 3; j++) {
          bezier_triangles[i][j] = indices[3 * i + j];
        }
        for (int j = 0; j < 3; j++) {
          const int v0 = indices[3 * i + j];
          const int v1 = indices[3 * i + (j + 1) % 3];
          if (v1 == std::get<0>(info[v0])) {
            const int i_id = std::get<1>(info[v0]);
            const int j_id = std::get<2>(info[v0]);
            const UnsignedIndex_t arc_id = list_of_closed_curves[i_id][j_id];
            bezier_triangles[i][3 + j] = points.size();
            boundaries[v0][0] = v0;
            boundaries[v0][1] = v1;
            boundaries[v0][2] = points.size();
            points.push_back(arc_list_m[arc_id].control_point());
            weights.push_back(arc_list_m[arc_id].weight());
          } else {
            const Pt& pt_0 = points[v0];
            const Pt& pt_1 = points[v1];
            const auto arc =
                computeVerticalRationalBezierArc(aligned_p, pt_0, pt_1);
            bezier_triangles[i][3 + j] = points.size();
            points.push_back(arc.control_point());
            weights.push_back(arc.weight());
          }
        }
      }

      const auto& datum = paraboloid_m.getDatum();
      const auto& frame = paraboloid_m.getReferenceFrame();
      for (int i = 0; i < points.size(); i++) {
        const Pt base_pt = points[i];
        points[i] = Pt(0.0, 0.0, 0.0);
        for (int d = 0; d < 3; ++d) {
          for (int n = 0; n < 3; ++n) {
            points[i][n] += frame[d][n] * base_pt[d];
          }
        }
        points[i] += datum;
      }

      bezier_surface.addPoints(points, weights);
      bezier_surface.addBezierTriangles(bezier_triangles);
      bezier_surface.addBoundaries(boundaries);
    } else if (a_order == 3) {
      std::vector<std::array<UnsignedIndex_t, 10>> bezier_triangles(ntriangles);
      const auto& aligned_p = paraboloid_m.getAlignedParaboloid();
      for (UnsignedIndex_t i = 0; i < indices.size() / 3; ++i) {
        auto V = Pt(0.0, 0.0, 0.0);
        auto E = Pt(0.0, 0.0, 0.0);
        for (int j = 0; j < 3; j++) {
          bezier_triangles[i][j] = indices[3 * i + j];
          const Pt& pt = points[indices[3 * i + j]];
          V += pt;
        }
        V *= 1.0 / 3.0;
        for (int j = 0; j < 3; j++) {
          const int v0 = indices[3 * i + j];
          const int v1 = indices[3 * i + (j + 1) % 3];
          const Pt& pt_0 = points[v0];
          const Pt& pt_2 = points[v1];
          if (v1 == std::get<0>(info[v0])) {
            const int i_id = std::get<1>(info[v0]);
            const int j_id = std::get<2>(info[v0]);
            const UnsignedIndex_t arc_id = list_of_closed_curves[i_id][j_id];
            const Pt& pt_1 = arc_list_m[arc_id].control_point();
            const double w = arc_list_m[arc_id].weight();
            const double w1 = (1.0 + 2.0 * w) / 3.0;
            const double w1_inv = 1.0 / (3.0 * w1);
            const Pt pt_1_n = (pt_0 + w * 2.0 * pt_1) * w1_inv;
            const Pt pt_2_n = (pt_2 + w * 2.0 * pt_1) * w1_inv;
            bezier_triangles[i][3 + 2 * j] = points.size();
            points.push_back(pt_1_n);
            weights.push_back(w1);
            bezier_triangles[i][3 + 2 * j + 1] = points.size();
            points.push_back(pt_2_n);
            weights.push_back(w1);
            E += pt_1_n + pt_2_n;
          } else {
            const auto arc =
                computeVerticalRationalBezierArc(aligned_p, pt_0, pt_2);
            const Pt& pt_1 = arc.control_point();
            const double w = arc.weight();
            const double w1 = (1.0 + 2.0 * w) / 3.0;
            const double w1_inv = 1.0 / (3.0 * w1);
            const Pt pt_1_n = (pt_0 + w * 2.0 * pt_1) * w1_inv;
            const Pt pt_2_n = (pt_2 + w * 2.0 * pt_1) * w1_inv;
            bezier_triangles[i][3 + 2 * j] = points.size();
            points.push_back(pt_1_n);
            weights.push_back(w1);
            bezier_triangles[i][3 + 2 * j + 1] = points.size();
            points.push_back(pt_2_n);
            weights.push_back(w1);
            E += pt_1_n + pt_2_n;
          }
        }
        E *= 1.0 / 6.0;
        Pt pt_ctrl = E + 0.5 * (E - V);
        bezier_triangles[i][9] = points.size();
        points.push_back(pt_ctrl);
        weights.push_back(1.0);
      }

      const auto& datum = paraboloid_m.getDatum();
      const auto& frame = paraboloid_m.getReferenceFrame();
      for (int i = 0; i < points.size(); i++) {
        const Pt base_pt = points[i];
        points[i] = Pt(0.0, 0.0, 0.0);
        for (int d = 0; d < 3; ++d) {
          for (int n = 0; n < 3; ++n) {
            points[i][n] += frame[d][n] * base_pt[d];
          }
        }
        points[i] += datum;
      }

      bezier_surface.addPoints(points, weights);
      bezier_surface.addBezierTriangles(bezier_triangles);
    }
  }

  return std::move(bezier_surface);
}

inline TriangulatedSurfaceOutput
ParaboloidParametrizedSurfaceOutput::triangulate(
    const double a_length_scale, const UnsignedIndex_t a_nsplit) const {
  TriangulatedSurfaceOutput returned_surface;
  this->triangulate_fromPtr(a_length_scale, a_nsplit, &returned_surface);
  return returned_surface;
}

inline void ParaboloidParametrizedSurfaceOutput::triangulate_fromPtr(
    const double a_length_scale, const UnsignedIndex_t a_nsplit,
    TriangulatedSurfaceOutput* returned_surface) const {
  const UnsignedIndex_t nArcs = this->size();
  double length_scale, length_scale_ref = length_scale_m;
  if (a_length_scale > 0.0) {
    length_scale_ref = a_length_scale;
  } else if (length_scale_ref <= 0.0) {
    auto surf = (*this);
    const double avg_length = std::sqrt(std::abs(surf.getSurfaceArea())) / 3.0;
    const double curv = std::fabs(surf.getAverageMeanCurvature());
    length_scale_ref = std::min(0.1 / curv, avg_length);
  }
  const auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();

  std::vector<std::vector<RationalBezierArc>> list_of_closed_curves;
  std::vector<bool> visited(nArcs, false);

  // First, we need to order the arcs so as to form closed
  // curves
  double min_arc_length = DBL_MAX;
  bool valid_curves = true;
  for (std::size_t t = 0; t < nArcs; ++t) {
    if (visited[t]) {
      continue;
    }
    visited[t] = true;
    // Start with next available arc
    if (arc_list_m[t].weight() > 1.0e15) {
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
    int counter = 0;
    while (end_id != start_id) {
      for (std::size_t e = t + 1; e < nArcs; ++e) {
        if (arc_list_m[e].start_point_id() == end_id) {
          visited[e] = true;
          if (arc_list_m[e].weight() > 1.0e15) {
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
          break;
        }
      }
      if (++counter > nArcs) {
        valid_curves = false;
        break;
      }
    }
  }

  returned_surface->clearAll();

  if (valid_curves) {
#ifdef IRL_USE_EARCUT
    // The number type to use for tessellation
    using Coord = double;
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

    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
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
          polygon[i].push_back({pt[0], pt[1]});
          signed_area +=
              0.5 * (previous_pt[0] * pt[1] - pt[0] * previous_pt[1]);
          previous_pt = pt;
        }
      }
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

    auto& vlist = returned_surface->getVertexList();
    auto& tlist = returned_surface->getTriangleList();
    std::vector<std::array<int, 4>> elist;

    UnsignedIndex_t count = 0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      count += polygon[i].size();
    }
    vlist.resize(count);

    count = 0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      for (UnsignedIndex_t j = 0; j < polygon[i].size(); ++j) {
        double x = polygon[i][j][0];
        double y = polygon[i][j][1];
        double z =
            -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
        vlist[count++] = Pt(x, y, z);
      }
    }

    reMeshPolygon(vlist, elist, indices, length_scale, aligned_paraboloid);

    count = 0;
    for (UnsignedIndex_t i = 0; i < indices.size() / 3; ++i) {
      if (indices[3 * i + 0] >= 0 && indices[3 * i + 1] >= 0 &&
          indices[3 * i + 2] >= 0) {
        count++;
      }
    }
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
    const auto& datum = paraboloid_m.getDatum();
    const auto& ref_frame = paraboloid_m.getReferenceFrame();
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

    CDT cdt;

    // std::ofstream myfile;
    // myfile.open("triangulation_log.txt");
    // myfile << "Starting triangulating surface.\n";
    // myfile << std::setprecision(16) << std::scientific
    //        << "Paraboloid: " << aligned_paraboloid << "\n";

    // Create boundaries
    std::vector<Point> points;
    std::list<Point> list_of_seeds;
    const UnsignedIndex_t nCurves =
        static_cast<UnsignedIndex_t>(list_of_closed_curves.size());
    UnsignedIndex_t start_points = 0;
    double total_signed_area = 0.0;
    double xmin = DBL_MAX, xmax = -DBL_MAX;
    double ymin = DBL_MAX, ymax = -DBL_MAX;
    UnsignedIndex_t vertex_count = 0;
    bool previous_valid = false;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      points.resize(0);
      const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
      // Loop over arcs of curve
      UnsignedIndex_t added_points = 0;
      double signed_area = 0.0;
      for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
        // Compute approximate arc length
        const RationalBezierArc& arc = list_of_closed_curves[i][j];
        const double arc_length = arc.arc_length();
        // myfile << std::setprecision(16) << std::scientific << "Curve " << i
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
          //        << vertex_count++ << " at " << pt[0] << ", " << pt[1] <<
          //        ".\n";
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
          //        << "Removing duplicate " << id1 << " at " << points[id1].x()
          //        << ", " << points[id1].y() << " too close to " << id0 << "
          //        at "
          //        << points[id0].x() << ", " << points[id0].y() << ".\n";
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
            // myfile << std::setprecision(16) << std::scientific << "Adding
            // hole "
            //        << hole_location[0] + (1.0e3 * DBL_EPSILON) * shift_dir[0]
            //        << ", "
            //        << hole_location[1] + (1.0e3 * DBL_EPSILON) * shift_dir[1]
            //        << ".\n";
            list_of_seeds.push_back(
                Point(hole_location[0] + (1.0e3 * DBL_EPSILON) * shift_dir[0],
                      hole_location[1] + (1.0e3 * DBL_EPSILON) * shift_dir[1]));
          }

          // Create segments
          // myfile << "Adding constraint " << points.size() - 1 << " -- " << 0
          //        << ".\n";
          cdt.insert_constraint(points[points.size() - 1], points[0]);

          for (UnsignedIndex_t j = 0; j < points.size() - 1; ++j) {
            // myfile << "Adding constraint " << j << " -- " << j + 1 << ".\n";
            cdt.insert_constraint(points[j], points[j + 1]);
          }
        }
        start_points += added_points;
        total_signed_area += 0.5 * signed_area;
      }
    }

    // myfile << "Surface has area " << total_signed_area << "\n";
    // myfile << "Mesh has " << cdt.number_of_vertices() << " vertices.\n";
    // myfile << "Refining with length-scale " << length_scale << ".\n";
    // sleep(1.0e-4);
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 CGAL::parameters::seeds(list_of_seeds)
                                     .criteria(Criteria(0.15, length_scale)));
    // , CGAL::parameters::seeds_are_in_domain(false));
    // myfile << "Mesh has " << cdt.number_of_vertices() << " vertices.\n";
    // myfile << "Mesh has " << cdt.number_of_faces() << " faces.\n";
    // CGAL::lloyd_optimize_mesh_2(cdt,
    //                             CGAL::parameters::max_iteration_number
    //                             = 20);
    auto& vlist = returned_surface->getVertexList();
    auto& tlist = returned_surface->getTriangleList();
    UnsignedIndex_t count = 0;
    CDT::Finite_faces_iterator face;
    // myfile << "Counting faces.\n";
    for (face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
         face++) {
      if (face->is_in_domain()) {
        count++;
      }
    }
    vlist.resize(3 * count);
    tlist.resize(count, TriangulatedSurfaceOutput::TriangleStorage::value_type::
                            fromNoExistencePlane(vlist, {0, 0, 0}));
    count = 0;
    // myfile << "Adding faces and vertices.\n";
    for (face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
         face++) {
      if (face->is_in_domain()) {
        // myfile << "Adding face " << count << ".\n";
        tlist[count] = TriangulatedSurfaceOutput::TriangleStorage::value_type::
            fromNoExistencePlane(vlist,
                                 {3 * count, 3 * count + 1, 3 * count + 2});
        for (UnsignedIndex_t d = 0; d < 3; d++) {
          const double x = CGAL::to_double(face->vertex(d)->point().x());
          const double y = CGAL::to_double(face->vertex(d)->point().y());
          const double z =
              -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
          vlist[3 * count + d] = Pt(x, y, z);
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
    // myfile << "Moving vertices in canonical frame.\n";
    const auto& datum = paraboloid_m.getDatum();
    const auto& ref_frame = paraboloid_m.getReferenceFrame();
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

    // myfile << "Finished triangulating surface.\n";
    // myfile.close();
#elif defined IRL_USE_TRIANGLE
    // Second, we approximate the arc length of the arc, so as
    // to know how many times it needs to be split
    std::vector<REAL> input_points;
    std::vector<REAL> input_holes;
    std::vector<int> input_segments;
    const UnsignedIndex_t nCurves =
        static_cast<UnsignedIndex_t>(list_of_closed_curves.size());
    // Loop over curves
    UnsignedIndex_t start_points = 0;
    double total_signed_area = 0.0;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
      // Loop over arcs of curve
      UnsignedIndex_t added_points = 0;
      double signed_area = 0.0;
      for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
        // Compute approximate arc length
        const RationalBezierArc& arc = list_of_closed_curves[i][j];
        // const auto& sp = arc.start_point();
        // const auto& ep = arc.start_point();
        // signed_area += (sp[0] * ep[1] - ep[0] * sp[1]);
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
        // added_points += nSplit;
        // const auto start_ind = input_points.size();
        // input_points.resize(start_ind + 2 * nSplit);
        // auto loc = input_points.begin() + start_ind;
        Pt previous_pt = arc.point(0.0);
        for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
          const double t = static_cast<double>(k) * step;
          const auto pt = arc.point(t);
          if (squaredMagnitude(pt - previous_pt) >
              1.0e8 * DBL_EPSILON * DBL_EPSILON) {
            input_points.push_back(pt[0]);
            input_points.push_back(pt[1]);
            previous_pt = pt;
            added_points++;
          }
        }
      }

      if (added_points >= 3) {
        signed_area += (input_points[start_points + 2 * added_points - 2] *
                            input_points[start_points + 1] -
                        input_points[start_points + 0] *
                            input_points[start_points + 2 * added_points - 1]);
        for (UnsignedIndex_t j = 0; j < added_points - 1; ++j) {
          signed_area += (input_points[start_points + 2 * j + 0] *
                              input_points[start_points + 2 * j + 3] -
                          input_points[start_points + 2 * j + 2] *
                              input_points[start_points + 2 * j + 1]);
        }

        if (nCurves > 1 && signed_area < 0.0) {
          // Add hole
          const auto p1x = input_points[start_points];
          const auto p1y = input_points[start_points + 1];
          const auto p2x = input_points[start_points + 2];
          const auto p2y = input_points[start_points + 3];
          std::array<double, 2> hole_location{
              {0.5 * (p1x + p2x), 0.5 * (p1y + p2y)}};
          Normal shift_dir = Normal(p2y - p1y, p1x - p2x, 0.0);
          shift_dir.normalize();
          const auto start_ind = input_holes.size();
          input_holes.resize(start_ind + 2);
          input_holes[start_ind] = 0.0;
          // hole_location[0] - (500.0 * DBL_EPSILON) *
          // shift_dir[0];
          input_holes[start_ind + 1] = 0.0;
          // hole_location[1] - (500.0 * DBL_EPSILON) *
          // shift_dir[1];
        }

        // Create segments
        const int seg_size = input_segments.size();
        input_segments.resize(seg_size + 2 * (added_points));
        auto seg_loc = input_segments.begin() + seg_size;
        *(seg_loc++) = start_points + added_points - 1;
        *(seg_loc++) = start_points;
        for (UnsignedIndex_t j = start_points;
             j < start_points + added_points - 1; ++j) {
          *(seg_loc++) = j;
          *(seg_loc++) = j + 1;
        }
        start_points += added_points;
        total_signed_area += 0.5 * signed_area;
      }
    }

    // Below section is for Triangle library
    if (input_points.size() > 0) {
      // std::cout << " Total area = " << total_signed_area <<
      // " compared to "
      //           << 2.0 * length_scale * length_scale <<
      //           std::endl;
      if (std::fabs(total_signed_area) > length_scale * length_scale) {
        // Calling triangulation library
        struct triangulateio in = {0}, out = {0};
        in.numberofpoints = input_points.size() / 2;
        in.pointlist = input_points.data();

        std::vector<int> pointmarkerlist(in.numberofpoints, 1);
        in.pointmarkerlist = pointmarkerlist.data();

        in.numberofsegments = input_segments.size() / 2;
        in.segmentlist = input_segments.data();
        std::vector<int> segmentmarkerlist(in.numberofsegments, 1);
        in.segmentmarkerlist = segmentmarkerlist.data();

        in.numberofholes = input_holes.size() / 2;
        if (in.numberofholes > 0) {
          in.holelist = input_holes.data();
        }

        char flags[50];
        sprintf(flags, "pza%.15feiQ", 0.5 * length_scale * length_scale);

        // std::cout << "Calling triangle with flags " <<
        // flags << " and with "
        //           << in.numberofpoints << " points and " <<
        //           in.numberofsegments
        //           << " segments and " << in.numberofholes
        //           << " holes and max area = "
        //           << 0.5 * length_scale * length_scale <<
        //           std::endl;

        // for (UnsignedIndex_t i = 0; i < in.numberofpoints;
        // ++i) {
        //   const double x = in.pointlist[2 * i + 0];
        //   const double y = in.pointlist[2 * i + 1];
        //   std::cout << "Point " << i << " = (" << x << ", "
        //   << y << ")"
        //             << std::endl;
        // }
        // for (UnsignedIndex_t i = 0; i <
        // in.numberofsegments; ++i) {
        //   const int j = in.segmentlist[2 * i + 0];
        //   const int k = in.segmentlist[2 * i + 1];
        //   std::cout << "Segment " << i << " = (" << j << ",
        //   " << k << ")"
        //             << std::endl;
        // }

        try {
          triangulate_from_lib(flags, &in, &out, (struct triangulateio*)NULL);
          // std::cout << "Triangle finished" << std::endl;

        } catch (std::runtime_error& e) {
          // std::cerr << e.what() << std::endl;
          // free(in.pointlist);
          free(in.pointattributelist);
          // free(in.pointmarkerlist);
          free(in.trianglelist);
          free(in.triangleattributelist);
          free(in.trianglearealist);
          free(in.neighborlist);
          // free(in.segmentlist);
          // free(in.segmentmarkerlist);
          // free(in.holelist);
          free(in.regionlist);
          free(in.edgelist);
          free(in.edgemarkerlist);
          free(in.normlist);
          free(out.pointlist);
          free(out.pointattributelist);
          free(out.pointmarkerlist);
          free(out.trianglelist);
          free(out.triangleattributelist);
          free(out.trianglearealist);
          free(out.neighborlist);
          free(out.segmentlist);
          free(out.segmentmarkerlist);
          free(out.regionlist);
          free(out.edgelist);
          free(out.edgemarkerlist);
          free(out.normlist);
        }

        auto& vlist = returned_surface->getVertexList();
        vlist.resize(out.numberofpoints);
        for (UnsignedIndex_t i = 0; i < out.numberofpoints; ++i) {
          const double x = out.pointlist[2 * i + 0];
          const double y = out.pointlist[2 * i + 1];
          const double z =
              -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
          vlist[i] = Pt(x, y, z);
        }

        // Translate and rotate triangulated surface vertices
        const auto& datum = paraboloid_m.getDatum();
        const auto& ref_frame = paraboloid_m.getReferenceFrame();
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

        for (UnsignedIndex_t i = 0; i < out.numberofedges; ++i) {
          if (out.edgemarkerlist[i] == 1) {
            returned_surface->addBoundaryEdge(out.edgelist[2 * i],
                                              out.edgelist[2 * i + 1]);
          }
        }

        auto& tlist = returned_surface->getTriangleList();
        tlist.resize(out.numberoftriangles,
                     TriangulatedSurfaceOutput::TriangleStorage::value_type::
                         fromNoExistencePlane(vlist, {0, 0, 0}));
        for (UnsignedIndex_t i = 0; i < out.numberoftriangles; ++i) {
          tlist[i] = TriangulatedSurfaceOutput::TriangleStorage::value_type::
              fromNoExistencePlane(
                  vlist,
                  {static_cast<UnsignedIndex_t>(out.trianglelist[3 * i + 0]),
                   static_cast<UnsignedIndex_t>(out.trianglelist[3 * i + 1]),
                   static_cast<UnsignedIndex_t>(out.trianglelist[3 * i + 2])});
        }

        /* free all allocated arrays, including those
         * allocated by Triangle.
         */
        free(out.pointlist);
        free(out.pointattributelist);
        free(out.pointmarkerlist);
        free(out.trianglelist);
        free(out.triangleattributelist);
        free(out.trianglearealist);
        free(out.neighborlist);
        free(out.segmentlist);
        free(out.segmentmarkerlist);
        free(out.regionlist);
        free(out.edgelist);
        free(out.edgemarkerlist);
        free(out.normlist);
        // free(in.pointlist);
        free(in.pointattributelist);
        // free(in.pointmarkerlist);
        free(in.trianglelist);
        free(in.triangleattributelist);
        free(in.trianglearealist);
        free(in.neighborlist);
        // free(in.segmentlist);
        // free(in.segmentmarkerlist);
        // free(in.holelist);
        free(in.regionlist);
        free(in.edgelist);
        free(in.edgemarkerlist);
        free(in.normlist);
      } else {  // Triangulate by hand
        auto& vlist = returned_surface->getVertexList();
        vlist.resize(input_points.size() / 2);
        for (UnsignedIndex_t i = 0; i < input_points.size() / 2; ++i) {
          const double x = input_points[2 * i + 0];
          const double y = input_points[2 * i + 1];
          const double z =
              -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
          vlist[i] = Pt(x, y, z);
        }

        // Translate and rotate triangulated surface vertices
        const auto& datum = paraboloid_m.getDatum();
        const auto& ref_frame = paraboloid_m.getReferenceFrame();
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

        returned_surface->addBoundaryEdge(input_points.size() / 2 - 1, 0);
        for (UnsignedIndex_t i = 0; i < input_points.size() / 2 - 1; ++i) {
          returned_surface->addBoundaryEdge(i, i + 1);
        }

        auto& tlist = returned_surface->getTriangleList();
        tlist.resize(input_points.size() / 2 - 2,
                     TriangulatedSurfaceOutput::TriangleStorage::value_type::
                         fromNoExistencePlane(vlist, {0, 0, 0}));
        for (UnsignedIndex_t i = 0; i < input_points.size() / 2 - 2; ++i) {
          tlist[i] = TriangulatedSurfaceOutput::TriangleStorage::value_type::
              fromNoExistencePlane(vlist, {0, i + 1, i + 2});
        }
      }
    }
#elif defined IRL_USE_GEOGRAM
    GEO::initialize();
    auto cdt = GEO::CDT2d();
    cdt.clear();
    cdt.create_enclosing_quad(GEO::vec2(-0.5, -1.5), GEO::vec2(1.5, -0.5),
                              GEO::vec2(1.5, 1.5), GEO::vec2(-0.5, 1.5));

    // Loop over curves
    const UnsignedIndex_t nCurves =
        static_cast<UnsignedIndex_t>(list_of_closed_curves.size());
    UnsignedIndex_t start_points = 0;
    std::vector<double> points;
    for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
      UnsignedIndex_t added_points = 0;
      const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
      // Loop over arcs of curve
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
        for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
          const double t = static_cast<double>(k) * step;
          const auto pt = arc.point(t);
          points.push_back(pt[0]);
          points.push_back(pt[1]);
          added_points++;
        }
      }
      std::cout << "Inserting " << added_points << " points " << std::endl;
      cdt.insert(static_cast<GEO::index_t>(added_points), points.data());

      if (added_points >= 3) {
        std::cout << "Inserting constraint (" << start_points + added_points - 1
                  << ", " << start_points << ")" << std::endl;
        cdt.insert_constraint(
            static_cast<GEO::index_t>(4 + start_points + added_points - 1),
            static_cast<GEO::index_t>(4 + start_points));
        for (UnsignedIndex_t j = start_points;
             j < start_points + added_points - 1; ++j) {
          std::cout << "Inserting constraint (" << j << ", " << j + 1 << ")"
                    << std::endl;
          cdt.insert_constraint(static_cast<GEO::index_t>(4 + j),
                                static_cast<GEO::index_t>(4 + j + 1));
        }
        start_points += added_points;
      }
    }

    cdt.remove_external_triangles();
    cdt.refine(18.0 * M_PI / 180.0, 0.00025);
    // cdt.remove_external_triangles();
    // cdt.remove_holes();

    auto& vlist = returned_surface->getVertexList();
    vlist.resize(cdt.nv());
    for (UnsignedIndex_t i = 0; i < cdt.nv(); ++i) {
      const double x = cdt.point(i).x;
      const double y = cdt.point(i).y;
      const double z =
          -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
      vlist[i] = Pt(x, y, z);
    }

    // Translate and rotate triangulated surface vertices
    const auto& datum = paraboloid_m.getDatum();
    const auto& ref_frame = paraboloid_m.getReferenceFrame();
    for (auto& vertex : vlist) {
      const Pt base_pt = vertex;
      vertex = Pt(0, 0, 0);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          vertex[n] += ref_frame[d][n] * base_pt[d];
        }
      }
      vertex += datum;
    }

    auto& tlist = returned_surface->getTriangleList();
    tlist.resize(cdt.nT(),
                 TriangulatedSurfaceOutput::TriangleStorage::value_type::
                     fromNoExistencePlane(vlist, {0, 0, 0}));
    for (UnsignedIndex_t i = 0; i < cdt.nT(); ++i) {
      tlist[i] = TriangulatedSurfaceOutput::TriangleStorage::value_type::
          fromNoExistencePlane(vlist,
                               {static_cast<UnsignedIndex_t>(cdt.Tv(i, 0)),
                                static_cast<UnsignedIndex_t>(cdt.Tv(i, 1)),
                                static_cast<UnsignedIndex_t>(cdt.Tv(i, 2))});
    }

#endif
  }
}

inline std::ostream& operator<<(
    std::ostream& out,
    const ParaboloidParametrizedSurfaceOutput& a_parametrized_surface) {
  const auto& aligned_paraboloid =
      a_parametrized_surface.getParaboloid().getAlignedParaboloid();
  out.precision(16);
  out << std::scientific << aligned_paraboloid.a() << " "
      << aligned_paraboloid.b() << std::endl;
  for (UnsignedIndex_t i = 0; i < a_parametrized_surface.size(); ++i) {
    out << a_parametrized_surface[i];
    if (i < a_parametrized_surface.size() - 1) out << std::endl;
  }
  return out;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_PARAMETRIZED_SURFACE_TPP_
