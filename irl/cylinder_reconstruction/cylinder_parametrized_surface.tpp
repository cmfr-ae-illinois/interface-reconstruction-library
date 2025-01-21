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

#define IRL_USE_EARCUT
// #define IRL_USE_TRIANGLE
// #define IRL_USE_CGAL
// #define IRL_USE_GEOGRAM

#include "irl/helpers/mymath.h"
#include "external/NumericalIntegration/NumericalIntegration.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"
#include "irl/cylinder_reconstruction/cylinder_parametrized_surface.h"

//#define VALDEBUG2

namespace IRL {

template <class VertexList>
void projectOnSurface(VertexList& vertices, const AlignedCylinder cylinder,
                      const UnsignedIndex_t fixed_vertices) {
  if (vertices.size() > fixed_vertices) {
    using ScalarType = AlignedCylinder::value_type;
    for (UnsignedIndex_t i = fixed_vertices; i < vertices.size(); i++) {
      const double y = vertices[i][1];
      vertices[i][2] = sqrt(maximum(cylinder.r() - cylinder.b() * y * y, ScalarType(0)));
    }
  }
}

template <class VertexList, class EdgeList, class TriList>
void reMeshPolygon(VertexList& vertices, EdgeList& edges, TriList& triangles,
                   const double length_scale,
                   const AlignedCylinder cylinder) {
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
    : ParametrizedSurfaceOutput{} {}

inline CylinderParametrizedSurfaceOutput::CylinderParametrizedSurfaceOutput(
    const Cylinder& a_cylinder)
    : cylinder_m{a_cylinder},
      ParametrizedSurfaceOutput{} {}

inline CylinderParametrizedSurfaceOutput::CylinderParametrizedSurfaceOutput(
    CylinderParametrizedSurfaceOutput&& a_rhs)
    : ParametrizedSurfaceOutput(a_rhs),
    cylinder_m(a_rhs.cylinder_m) {}

inline CylinderParametrizedSurfaceOutput::CylinderParametrizedSurfaceOutput(
    const CylinderParametrizedSurfaceOutput& a_rhs)
    : ParametrizedSurfaceOutput(a_rhs),
    cylinder_m(a_rhs.cylinder_m) {}

inline CylinderParametrizedSurfaceOutput& CylinderParametrizedSurfaceOutput::operator=(
    CylinderParametrizedSurfaceOutput&& a_rhs) {
  if (this != &a_rhs) {
    cylinder_m = a_rhs.cylinder_m;
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

inline CylinderParametrizedSurfaceOutput& CylinderParametrizedSurfaceOutput::operator=(
    const CylinderParametrizedSurfaceOutput& a_rhs) {
  if (this != &a_rhs) {
    cylinder_m = a_rhs.cylinder_m;
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
  cylinder_m = a_cylinder;
}

inline const Cylinder& CylinderParametrizedSurfaceOutput::getCylinder(void) const {
  return cylinder_m;
}

inline CylinderParametrizedSurfaceOutput::~CylinderParametrizedSurfaceOutput(void) {
  for (auto elem : pt_from_bezier_split_m) {
    delete elem;
  }
}

// class ArcContributionToCylinderSurfaceArea_Functor {
//  public:
//   ArcContributionToCylinderSurfaceArea_Functor(const RationalBezierArc& a_arc,
//                                        const AlignedCylinder& a_cylinder)
//       : arc_m(a_arc), cylinder_m(a_cylinder) {}

//   double operator()(double a_t) const {
//     const auto weight = arc_m.weight();
//     if (weight > 1.0e15) {
//       const auto pt0 = arc_m.point(0.5 * a_t);
//       const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
//       const auto der0 = arc_m.derivative(0.5 * a_t);
//       const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
//       const double b = cylinder_m.b();
//       const double r = cylinder_m.r();

//       const double primitive0 = 
//         pt0[0] * std::sqrt(r + (b - 1) * b * (pt0[1] * pt0[1])) /
//           sqrt(r - b * (pt0[1] * pt0[1]));

//       const double primitive1 = 
//         pt1[0] * std::sqrt(r + (b - 1) * b * (pt1[1] * pt1[1])) /
//           sqrt(r - b * (pt1[1] * pt1[1]));

//       return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);

//       //if (std::fabs(a) < 10.0 * DBL_EPSILON &&
//       //     std::fabs(b) < 10.0 * DBL_EPSILON) {
//       //   return 0.0;
//       // } else if (std::fabs(a) > std::fabs(b)) {
//       //   const double primitive0 =
//       //       (2. * pt0[0] *
//       //            std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
//       //                      4. * (b * b) * (pt0[1] * pt0[1])) -
//       //        (1. + 4. * (b * b) * (pt0[1] * pt0[1])) *
//       //            std::log(-2. * std::fabs(a) * pt0[0] +
//       //                     std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
//       //                               4. * (b * b) * (pt0[1] * pt0[1]))) /
//       //            std::fabs(a)) /
//       //       4.;
//       //   const double primitive1 =
//       //       (2. * pt1[0] *
//       //            std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
//       //                      4. * (b * b) * (pt1[1] * pt1[1])) -
//       //        (1. + 4. * (b * b) * (pt1[1] * pt1[1])) *
//       //            std::log(-2. * std::fabs(a) * pt1[0] +
//       //                     std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
//       //                               4. * (b * b) * (pt1[1] * pt1[1]))) /
//       //            std::fabs(a)) /
//       //       4.;
//       //   return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
//       // } else {
//       //   const double primitive0 =
//       //       -(2. * pt0[1] *
//       //             std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
//       //                       4. * (b * b) * (pt0[1] * pt0[1])) -
//       //         (1. + 4. * (a * a) * (pt0[0] * pt0[0])) *
//       //             std::log(-2. * std::fabs(b) * pt0[1] +
//       //                      std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
//       //                                4. * (b * b) * (pt0[1] * pt0[1]))) /
//       //             std::fabs(b)) /
//       //       (4.);
//       //   const double primitive1 =
//       //       -(2. * pt1[1] *
//       //             std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
//       //                       4. * (b * b) * (pt1[1] * pt1[1])) -
//       //         (1. + 4. * (a * a) * (pt1[0] * pt1[0])) *
//       //             std::log(-2. * std::fabs(b) * pt1[1] +
//       //                      std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
//       //                                4. * (b * b) * (pt1[1] * pt1[1]))) /
//       //             std::fabs(b)) /
//       //       (4.);
//       //   return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
//       // }
//     } else {
//       const auto pt = arc_m.point(a_t);
//       const auto der = arc_m.derivative(a_t);
//       const double b = cylinder_m.b();
//       const double r = cylinder_m.r();
//       const double primitive =
//           pt[0] * std::sqrt(r + (b - 1) * b * (pt[1] * pt[1])) /
//         sqrt(r - b * (pt[1] * pt[1]));
//       if (std::isnan(primitive)) {
//         std::cout << "Pr = " << pt << std::endl;
//         std::cout << "Der = " << der << std::endl;
//         std::cout << "b = " << b << std::endl;
//         std::cout << "r = " << r << std::endl;
//         std::cout << "Arc: weight = " << arc_m.weight() << std::endl;
//         std::cout << "Arc: start = " << arc_m.start_point() << std::endl;
//         std::cout << "Arc: ctrl  = " << arc_m.control_point() << std::endl;
//         std::cout << "Arc: end   = " << arc_m.end_point() << std::endl;
//         std::cout << "Primitive is NaN" << std::endl;
//         exit(1);
//       }
//       if (std::isnan(der[0])) {
//         std::cout << "der[0] is NaN" << std::endl;
//         exit(1);
//       }

//       return primitive * der[0];
//     }
//   }

//  private:
//   const RationalBezierArc& arc_m;
//   const AlignedCylinder& cylinder_m;
// };

class ArcContributionToCylinderNormalX_Functor {
 public:
  ArcContributionToCylinderNormalX_Functor(const RationalBezierArc& a_arc,
                                   const AlignedCylinder& a_cylinder)
      : arc_m(a_arc), cylinder_m(a_cylinder) {}

  double operator()(double a_t) const {
    // there is no X component to the Normal vecteur
    return 0.;
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedCylinder& cylinder_m;
};

class ArcContributionToCylinderNormalY_Functor {
 public:
  ArcContributionToCylinderNormalY_Functor(const RationalBezierArc& a_arc,
                                   const AlignedCylinder& a_cylinder)
      : arc_m(a_arc), cylinder_m(a_cylinder) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double b = cylinder_m.b();
      const double r = cylinder_m.r();
      const double primitive0 = b * pt0[1] * pt0[0] / 
        std::sqrt(r - b * (pt0[1] * pt0[1]));
      const double primitive1 = b * pt1[1] * pt1[0] / 
        std::sqrt(r - b * (pt1[1] * pt1[1]));
      return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double b = cylinder_m.b();
      const double r = cylinder_m.r();
      const double primitive = b * pt[1] * pt[0] / 
        std::sqrt(r - b * (pt[1] * pt[1]));
      return primitive * der[0];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedCylinder& cylinder_m;
};

class ArcContributionToCylinderNormalZ_Functor {
 public:
  ArcContributionToCylinderNormalZ_Functor(const RationalBezierArc& a_arc,
                                   const AlignedCylinder& a_cylinder)
      : arc_m(a_arc), cylinder_m(a_cylinder) {}

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
  const AlignedCylinder& cylinder_m;
};

// class ArcContributionToMeanCurvature_Functor {
//  public:
//   ArcContributionToMeanCurvature_Functor(const RationalBezierArc& a_arc,
//                                          const AlignedParaboloid& a_paraboloid)
//       : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

//   double operator()(double a_t) const {
//     const auto weight = arc_m.weight();
//     if (weight > 1.0e15) {
//       const auto pt0 = arc_m.point(0.5 * a_t);
//       const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
//       const auto der0 = arc_m.derivative(0.5 * a_t);
//       const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
//       const double a = paraboloid_m.a();
//       const double b = paraboloid_m.b();
//       if (std::fabs(a) < 10.0 * DBL_EPSILON &&
//           std::fabs(b) < 10.0 * DBL_EPSILON) {
//         return 0.0;
//       } else if (std::fabs(a) > std::fabs(b)) {
//         const double primitive0 =
//             2. * b * pt0[0] +
//             ((a + 4. * a * (b * b) * (pt0[1] * pt0[1]) -
//               4. * (b * b * b) * (pt0[1] * pt0[1])) *
//              std::atan((2. * a * pt0[0]) /
//                        std::sqrt(1. + 4. * (b * b) * (pt0[1] * pt0[1])))) /
//                 (a * std::sqrt(1. + 4. * (b * b) * (pt0[1] * pt0[1])));
//         const double primitive1 =
//             2. * b * pt1[0] +
//             ((a + 4. * a * (b * b) * (pt1[1] * pt1[1]) -
//               4. * (b * b * b) * (pt1[1] * pt1[1])) *
//              std::atan((2. * a * pt1[0]) /
//                        std::sqrt(1. + 4. * (b * b) * (pt1[1] * pt1[1])))) /
//                 (a * std::sqrt(1. + 4. * (b * b) * (pt1[1] * pt1[1])));
//         return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
//       } else {
//         const double primitive0 =
//             -2. * a * pt0[1] -
//             ((b - 4. * (a * a * a) * (pt0[0] * pt0[0]) +
//               4. * (a * a) * b * (pt0[0] * pt0[0])) *
//              std::atan((2. * b * pt0[1]) /
//                        std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0])))) /
//                 (b * std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0])));
//         const double primitive1 =
//             -2. * a * pt1[1] -
//             ((b - 4. * (a * a * a) * (pt1[0] * pt1[0]) +
//               4. * (a * a) * b * (pt1[0] * pt1[0])) *
//              std::atan((2. * b * pt1[1]) /
//                        std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0])))) /
//                 (b * std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0])));
//         return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
//       }
//     } else {
//       const auto pt = arc_m.point(a_t);
//       const auto der = arc_m.derivative(a_t);
//       const double a = paraboloid_m.a();
//       const double b = paraboloid_m.b();
//       if (std::fabs(a) < 10.0 * DBL_EPSILON &&
//           std::fabs(b) < 10.0 * DBL_EPSILON) {
//         return 0.0;
//       } else if (std::fabs(a) > std::fabs(b)) {
//         const double primitive =
//             2. * b * pt[0] +
//             ((a + 4. * a * (b * b) * (pt[1] * pt[1]) -
//               4. * (b * b * b) * (pt[1] * pt[1])) *
//              std::atan((2. * a * pt[0]) /
//                        std::sqrt(1. + 4. * (b * b) * (pt[1] * pt[1])))) /
//                 (a * std::sqrt(1. + 4. * (b * b) * (pt[1] * pt[1])));
//         return primitive * der[1];
//       } else {
//         const double primitive =
//             -2. * a * pt[1] -
//             ((b - 4. * (a * a * a) * (pt[0] * pt[0]) +
//               4. * (a * a) * b * (pt[0] * pt[0])) *
//              std::atan((2. * b * pt[1]) /
//                        std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0])))) /
//                 (b * std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0])));
//         return primitive * der[0];
//       }
//     }
//   }

//  private:
//   const RationalBezierArc& arc_m;
//   const AlignedParaboloid& paraboloid_m;
// };  // namespace IRL

// class ArcContributionToGaussianCurvature_Functor {
//  public:
//   ArcContributionToGaussianCurvature_Functor(
//       const RationalBezierArc& a_arc, const AlignedParaboloid& a_paraboloid)
//       : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

//   double operator()(double a_t) const {
//     const auto weight = arc_m.weight();
//     if (weight > 1.0e15) {
//       const auto pt0 = arc_m.point(0.5 * a_t);
//       const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
//       const auto der0 = arc_m.derivative(0.5 * a_t);
//       const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
//       const double a = paraboloid_m.a();
//       const double b = paraboloid_m.b();
//       if (std::fabs(a) < 10.0 * DBL_EPSILON &&
//           std::fabs(b) < 10.0 * DBL_EPSILON) {
//         return 0.0;
//       } else {
//         const double primitive0 =
//             4.0 * a * b * pt0[0] /
//             ((1.0 + 4.0 * b * b * pt0[1] * pt0[1]) *
//              std::sqrt(1.0 + 4.0 * a * a * pt0[0] * pt0[0] +
//                        4.0 * b * b * pt0[1] * pt0[1]));
//         const double primitive1 =
//             4.0 * a * b * pt1[0] /
//             ((1.0 + 4.0 * b * b * pt1[1] * pt1[1]) *
//              std::sqrt(1.0 + 4.0 * a * a * pt1[0] * pt1[0] +
//                        4.0 * b * b * pt1[1] * pt1[1]));
//         return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
//       }
//     } else {
//       const auto pt = arc_m.point(a_t);
//       const auto der = arc_m.derivative(a_t);
//       const double a = paraboloid_m.a();
//       const double b = paraboloid_m.b();
//       if (std::fabs(a) < 10.0 * DBL_EPSILON &&
//           std::fabs(b) < 10.0 * DBL_EPSILON) {
//         return 0.0;
//       } else {
//         const double primitive = 4.0 * a * b * pt[0] /
//                                  ((1.0 + 4.0 * b * b * pt[1] * pt[1]) *
//                                   std::sqrt(1.0 + 4.0 * a * a * pt[0] * pt[0] +
//                                             4.0 * b * b * pt[1] * pt[1]));
//         return primitive * der[1];
//       }
//     }
//   }

//  private:
//   const RationalBezierArc& arc_m;
//   const AlignedParaboloid& paraboloid_m;
// };  // namespace IRL

inline double CylinderParametrizedSurfaceOutput::getSurfaceArea(void) {
  // if (!knows_surface_area_m) {
  //   const UnsignedIndex_t nArcs = this->size();
  //   surface_area_m = 0.0;
  //   size_t limit = 128;

  //   const double epsabs = 10.0 * DBL_EPSILON;
  //   const double epsrel = 0.0;
  //   auto& aligned_cylinder = cylinder_m.getAlignedCylinder();
  //   for (std::size_t t = 0; t < nArcs; ++t) {
  //     // Define the functor
  //     ArcContributionToCylinderSurfaceArea_Functor functor(arc_list_m[t],
  //                                                  aligned_cylinder);

  //     // Define the integrator.
  //     Eigen::Integrator<double> integrator(limit);

  //     // Define a quadrature rule.
  //     Eigen::Integrator<double>::QuadratureRule quadrature_rule =
  //         Eigen::Integrator<double>::GaussKronrod61;

  //     // Integrate.
  //     surface_area_m += integrator.quadratureAdaptive(functor, 0.0, 1.0, epsabs,
  //                                                     epsrel, quadrature_rule);
  //   }
  //   knows_surface_area_m = true;
  // }
  // return surface_area_m;
  std::cout << "calling a function that was not implemented : CylinderParametrizedSurfaceOutput::getSurfaceArea, 0 value returned" << std::endl;
  return 0;
}

inline Normal CylinderParametrizedSurfaceOutput::getAverageNormal(void) {
  if (!knows_avg_normal_m) {
    const UnsignedIndex_t nArcs = this->size();
    avg_normal_m = Normal();
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = cylinder_m.getAlignedCylinder();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToCylinderNormalX_Functor functorx(arc_list_m[t],
                                                aligned_paraboloid);
      ArcContributionToCylinderNormalY_Functor functory(arc_list_m[t],
                                                aligned_paraboloid);
      ArcContributionToCylinderNormalZ_Functor functorz(arc_list_m[t],
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
  std::cout << "calling a function that is wrong : CylinderParametrizedSurfaceOutput::getAverageNormal, wrong value returned" << std::endl;
  return avg_normal_m;
}

inline Normal CylinderParametrizedSurfaceOutput::getAverageNormalNonAligned(void) {
  auto aligned_normal = this->getAverageNormal();
  const auto& ref_frame = this->getCylinder().getReferenceFrame();
  auto normal = Normal();
  for (std::size_t d = 0; d < 3; ++d) {
    for (std::size_t n = 0; n < 3; ++n) {
      normal[n] += ref_frame[d][n] * aligned_normal[d];
    }
  }
  return normal;
}

inline double CylinderParametrizedSurfaceOutput::getMeanCurvatureIntegral(void) { return 0; }
// inline double ParaboloidParametrizedSurfaceOutput::getMeanCurvatureIntegral(void) {
//   if (!knows_int_mean_curv_m) {
//     const UnsignedIndex_t nArcs = this->size();
//     int_mean_curv_m = 0.0;
//     size_t limit = 128;

//     const double epsabs = 10.0 * DBL_EPSILON;
//     const double epsrel = 0.0;
//     auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
//     for (std::size_t t = 0; t < nArcs; ++t) {
//       // Define the functor
//       ArcContributionToMeanCurvature_Functor functor(arc_list_m[t],
//                                                      aligned_paraboloid);

//       // Define the integrator.
//       Eigen::Integrator<double> integrator(limit);

//       // Define a quadrature rule.
//       Eigen::Integrator<double>::QuadratureRule quadrature_rule =
//           Eigen::Integrator<double>::GaussKronrod61;

//       // Integrate.
//       int_mean_curv_m += integrator.quadratureAdaptive(
//           functor, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
//     }
//     knows_int_mean_curv_m = true;
//   }
//   return int_mean_curv_m;
// }

// inline double ParaboloidParametrizedSurfaceOutput::getGaussianCurvatureIntegral(void) {
//   if (!knows_int_gaussian_curv_m) {
//     const UnsignedIndex_t nArcs = this->size();
//     int_gaussian_curv_m = 0.0;
//     size_t limit = 128;

//     const double epsabs = 10.0 * DBL_EPSILON;
//     const double epsrel = 0.0;
//     auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
//     for (std::size_t t = 0; t < nArcs; ++t) {
//       // Define the functor
//       ArcContributionToGaussianCurvature_Functor functor(arc_list_m[t],
//                                                          aligned_paraboloid);

//       // Define the integrator.
//       Eigen::Integrator<double> integrator(limit);

//       // Define a quadrature rule.
//       Eigen::Integrator<double>::QuadratureRule quadrature_rule =
//           Eigen::Integrator<double>::GaussKronrod61;

//       // Integrate.
//       int_gaussian_curv_m += integrator.quadratureAdaptive(
//           functor, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
//     }
//     knows_int_gaussian_curv_m = true;
//   }
//   return int_gaussian_curv_m;
// }
inline double CylinderParametrizedSurfaceOutput::getGaussianCurvatureIntegral(void) {
  std::cout << "calling a function that was not implemented : CylinderParametrizedSurfaceOutput::getGaussianCurvatureIntegral, 0 value returned" << std::endl;
  return 0;
}

inline Normal CylinderParametrizedSurfaceOutput::getNormalAligned(const Pt a_pt) {
  auto& aligned_cylinder = this->getCylinder().getAlignedCylinder();
  auto aligned_normal = getCylinderSurfaceNormal(aligned_cylinder, a_pt);
  aligned_normal.normalize();
  return aligned_normal;
}

inline Normal CylinderParametrizedSurfaceOutput::getNormalNonAligned(const Pt a_pt) {
  const auto& datum = this->getCylinder().getDatum();
  const auto& ref_frame = this->getCylinder().getReferenceFrame();
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

// inline double CylinderParametrizedSurfaceOutput::getMeanCurvatureAligned(
//     const Pt a_pt) {
//   auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
//   return (2. * (aligned_paraboloid.a() + aligned_paraboloid.b() +
//                 4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
//                     aligned_paraboloid.b() * (a_pt[0] * a_pt[0]) +
//                 4. * aligned_paraboloid.a() *
//                     (aligned_paraboloid.b() * aligned_paraboloid.b()) *
//                     (a_pt[1] * a_pt[1]))) /
//          std::pow(1. +
//                       4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
//                           (a_pt[0] * a_pt[0]) +
//                       4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
//                           (a_pt[1] * a_pt[1]),
//                   1.5);
// }
inline double CylinderParametrizedSurfaceOutput::getMeanCurvatureAligned(
     const Pt a_pt) {
      std::cout << "calling a function that was not implemented : CylinderParametrizedSurfaceOutput::getMeanCurvatureAligned, 0 value returned" << std::endl;
      return 0;
     }

// inline double ParaboloidParametrizedSurfaceOutput::getMeanCurvatureNonAligned(
//     const Pt a_pt) {
//   const auto& datum = this->getParaboloid().getDatum();
//   const auto& ref_frame = this->getParaboloid().getReferenceFrame();
//   // assert(ref_frame.isOrthonormalBasis());
//   const Pt original_pt = a_pt - datum;
//   auto aligned_pt = a_pt;
//   for (std::size_t n = 0; n < 3; ++n) {
//     aligned_pt[n] = ref_frame[n] * original_pt;
//   }
//   return this->getMeanCurvatureAligned(aligned_pt);
// }
inline double CylinderParametrizedSurfaceOutput::getMeanCurvatureNonAligned(
     const Pt a_pt) {
      std::cout << "calling a function that was not implemented : CylinderParametrizedSurfaceOutput::getMeanCurvatureNonAligned, 0 value returned" << std::endl;
      return 0;
     }

// inline double ParaboloidParametrizedSurfaceOutput::getGaussianCurvatureAligned(
//     const Pt a_pt) {
//   auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
//   return 4. * aligned_paraboloid.a() * aligned_paraboloid.b() /
//          ((1. +
//            4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
//                (a_pt[0] * a_pt[0]) +
//            4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
//                (a_pt[1] * a_pt[1])) *
//           (1. +
//            4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
//                (a_pt[0] * a_pt[0]) +
//            4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
//                (a_pt[1] * a_pt[1])));
// }
inline double CylinderParametrizedSurfaceOutput::getGaussianCurvatureAligned(
     const Pt a_pt) {
      std::cout << "calling a function that was not implemented : CylinderParametrizedSurfaceOutput::getGaussianCurvatureAligned, 0 value returned" << std::endl;
      return 0;
     }

// inline double ParaboloidParametrizedSurfaceOutput::getGaussianCurvatureNonAligned(
//     const Pt a_pt) {
//   const auto& datum = this->getParaboloid().getDatum();
//   const auto& ref_frame = this->getParaboloid().getReferenceFrame();
//   // assert(ref_frame.isOrthonormalBasis());
//   const Pt original_pt = a_pt - datum;
//   auto aligned_pt = a_pt;
//   for (std::size_t n = 0; n < 3; ++n) {
//     aligned_pt[n] = ref_frame[n] * original_pt;
//   }
//   return this->getGaussianCurvatureAligned(aligned_pt);
// }
inline double CylinderParametrizedSurfaceOutput::getGaussianCurvatureNonAligned(
     const Pt a_pt) {
      std::cout << "calling a function that was not implemented : CylinderParametrizedSurfaceOutput::getGaussianCurvatureNonAligned, 0 value returned" << std::endl;
      return 0;
     }

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
  if (a_length_scale > 0.0) {
    length_scale_ref = a_length_scale;
  } else if (length_scale_ref <= 0.0) {
    auto surf = (*this);
    const double avg_length = std::sqrt(std::abs(surf.getSurfaceArea())) / 3.0;
    const double curv = std::fabs(surf.getAverageMeanCurvature());
    length_scale_ref = std::min(0.1 / curv, avg_length);
  }
  const auto& aligned_cylinder = cylinder_m.getAlignedCylinder();

  #ifdef VALDEBUG2
  std::cout << "triangulating a cylinder surface\nThere is " << nArcs << " arcs\nLet's find the close curves\n";
  #endif
  
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
  #ifdef VALDEBUG2
  std::cout << "starting a curve with an arc going from " << arc_list_m[t].start_point() << " to " << arc_list_m[t].end_point() << std::endl;
  #endif
    const std::uintptr_t start_id = arc_list_m[t].start_point_id();
    std::uintptr_t end_id = arc_list_m[t].end_point_id();
    int counter = 0;
    while (end_id != start_id) {
      for (std::size_t e = t + 1; e < nArcs; ++e) {
        if (arc_list_m[e].start_point_id() == end_id) {
          visited[e] = true;
  #ifdef VALDEBUG2
  std::cout << "next curve is an arc going from " << arc_list_m[e].start_point() << " to " << arc_list_m[e].end_point() << std::endl;
  #endif
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
  #ifdef VALDEBUG2
  std::cout << "end of that close curve\n";
  #endif
  }
  #ifdef VALDEBUG2
  std::cout << "in the end, there is " << list_of_closed_curves.size() << "close curves\n";
  #endif

  returned_surface->clearAll();

  if (valid_curves) {
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

  #ifdef VALDEBUG2
  std::cout << "let's define out polygons, there are " << nCurves << " of them\n";
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
  std::cout << "polygon number " << i << " has " << polygon[i].size() << " points\n";
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
        double z =
            sqrt(maximum(aligned_cylinder.r() - aligned_cylinder.b() * y * y, double(0)));
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
  std::cout << "after remesh, vlist has " << vlist.size() << " elements, and count has a value of : " << count << "\n";
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
    const auto& datum = cylinder_m.getDatum();
    const auto& ref_frame = cylinder_m.getReferenceFrame();
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
  }
}

inline std::ostream& operator<<(
    std::ostream& out,
    const CylinderParametrizedSurfaceOutput& a_parametrized_surface) {
  const auto& aligned_cylinder =
      a_parametrized_surface.getCylinder().getAlignedCylinder();
  out.precision(16);
  out << std::scientific << aligned_cylinder.b() << " "
      << aligned_cylinder.r() << std::endl;
  for (UnsignedIndex_t i = 0; i < a_parametrized_surface.size(); ++i) {
    out << a_parametrized_surface[i];
    if (i < a_parametrized_surface.size() - 1) out << std::endl;
  }
  return out;
}

}  // namespace IRL

#endif  // IRL_CYLINDER_RECONSTRUCTION_CYLINDER_PARAMETRIZED_SURFACE_TPP_
