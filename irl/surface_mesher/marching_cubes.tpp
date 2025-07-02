// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SURFACE_MESH_MARCHING_CUBES_TPP_
#define IRL_SURFACE_MESH_MARCHING_CUBES_TPP_

namespace IRL {

inline MarchingCubes::MarchingCubes(const RectangularCuboid a_domain,
                                    const ImplicitF a_function)
    : function_m(a_function), domain_m(a_domain) {}

inline void MarchingCubes::setDomain(const RectangularCuboid a_domain) {
  domain_m = a_domain;
}

inline void MarchingCubes::setFunction(const ImplicitF a_function) {
  function_m = a_function;
}

inline Pt VertexInterp(const MarchingCubes::ImplicitF function, const Pt p1,
                       const Pt p2, double valp1, double valp2) {
  return Pt(p1 + (p1 - p2) * (valp1 / (valp2 - valp1)));
}

// This is based off Paul Bourke's code at
// https://paulbourke.net/geometry/polygonise/
inline TriangulatedSurfaceOutput MarchingCubes::triangulate(
    const UnsignedIndex_t a_subdivision) {
  TriangulatedSurfaceOutput surface;
  surface.clearAll();
  const UnsignedIndex_t n = a_subdivision > 0 ? a_subdivision + 1 : 1;
  std::vector<long int> edgx((n - 1) * n * n, -1);
  std::vector<long int> edgy((n - 1) * n * n, -1);
  std::vector<long int> edgz((n - 1) * n * n, -1);
  std::vector<std::vector<std::vector<double>>> f(
      n, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));
  std::array<UnsignedIndex_t, 12> vertlist;
  std::array<double, 8> vals;
  std::array<Pt, 8> p;

  const Pt x0 = domain_m.getLowerLimits();
  const Pt x1 = domain_m.getUpperLimits();
  const Pt dx = (x1 - x0) * (1.0 / static_cast<double>(n - 1));
  for (UnsignedIndex_t i = 0; i < n; i++) {
    for (UnsignedIndex_t j = 0; j < n; j++) {
      for (UnsignedIndex_t k = 0; k < n; k++) {
        const double x = static_cast<double>(i) * dx[0];
        const double y = static_cast<double>(j) * dx[1];
        const double z = static_cast<double>(k) * dx[2];
        f[i][j][k] = function_m(Pt(x, y, z));
      }
    }
  }

  for (UnsignedIndex_t i = 0; i < n - 1; i++) {
    for (UnsignedIndex_t j = 0; j < n - 1; j++) {
      for (UnsignedIndex_t k = 0; k < n - 1; k++) {
        vals[0] = f[i][j][k];
        vals[1] = f[i + 1][j][k];
        vals[2] = f[i + 1][j][k + 1];
        vals[3] = f[i][j][k + 1];
        vals[4] = f[i][j + 1][k];
        vals[5] = f[i + 1][j + 1][k];
        vals[6] = f[i + 1][j + 1][k + 1];
        vals[7] = f[i][j + 1][k + 1];

        p[0].x() = p[3].x() = p[4].x() = p[7].x() =
            static_cast<double>(i) * dx[0];
        p[1].x() = p[2].x() = p[5].x() = p[6].x() =
            static_cast<double>(i + 1) * dx[0];
        p[0].y() = p[1].y() = p[2].y() = p[3].y() =
            static_cast<double>(j) * dx[1];
        p[4].y() = p[5].y() = p[6].y() = p[7].y() =
            static_cast<double>(j + 1) * dx[1];
        p[0].z() = p[1].z() = p[4].z() = p[5].z() =
            static_cast<double>(k) * dx[2];
        p[2].z() = p[3].z() = p[6].z() = p[7].z() =
            static_cast<double>(k + 1) * dx[2];

        UnsignedIndex_t cubeindex = 0;
        if (vals[0] < 0.0) cubeindex |= 1;
        if (vals[1] < 0.0) cubeindex |= 2;
        if (vals[2] < 0.0) cubeindex |= 4;
        if (vals[3] < 0.0) cubeindex |= 8;
        if (vals[4] < 0.0) cubeindex |= 16;
        if (vals[5] < 0.0) cubeindex |= 32;
        if (vals[6] < 0.0) cubeindex |= 64;
        if (vals[7] < 0.0) cubeindex |= 128;

        /* Cube is not entirely in/out of the surface */
        if (edgeTable[cubeindex] != 0) {
          /* Find the vertex where the surface intersects the cube
           * if vertice has not been generated before, then generates it */
          if (edgeTable[cubeindex] & 1) {
            /* Edge 0 follows axis X with coordinates {i, j, k} */
            UnsignedIndex_t edgeindex = i + j * (n - 1) + k * (n - 1) * n;
            if (edgx[edgeindex] == -1) {
              edgx[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[0], p[1], vals[0], vals[1]));
            }
            vertlist[0] = edgx[edgeindex];
          }
          if (edgeTable[cubeindex] & 2) {
            /* Edge 1 follows axis Z with coordinates {i + 1, j, k} */
            UnsignedIndex_t edgeindex = i + 1 + j * n + k * n * n;
            if (edgz[edgeindex] == -1) {
              edgz[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[1], p[2], vals[1], vals[2]));
            }
            vertlist[1] = edgz[edgeindex];
          }
          if (edgeTable[cubeindex] & 4) {
            /* Edge 2 follows axis X with coordinates {i, j, k + 1} */
            UnsignedIndex_t edgeindex = i + j * (n - 1) + (k + 1) * (n - 1) * n;
            if (edgx[edgeindex] == -1) {
              edgx[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[3], p[2], vals[3], vals[2]));
            }
            vertlist[2] = edgx[edgeindex];
          }
          if (edgeTable[cubeindex] & 8) {
            /* Edge 3 follows axis Z with coordinates {i, j, k} */
            UnsignedIndex_t edgeindex = i + j * n + k * n * n;
            if (edgz[edgeindex] == -1) {
              edgz[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[0], p[3], vals[0], vals[3]));
            }
            vertlist[3] = edgz[edgeindex];
          }
          if (edgeTable[cubeindex] & 16) {
            /* Edge 4 follows axis X with coordinates {i, j + 1, k} */
            UnsignedIndex_t edgeindex = i + (j + 1) * (n - 1) + k * (n - 1) * n;
            if (edgx[edgeindex] == -1) {
              edgx[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[4], p[5], vals[4], vals[5]));
            }
            vertlist[4] = edgx[edgeindex];
          }
          if (edgeTable[cubeindex] & 32) {
            /* Edge 5 follows axis Z with coordinates {i + 1, j + 1, k} */
            UnsignedIndex_t edgeindex = i + 1 + (j + 1) * n + k * n * n;
            if (edgz[edgeindex] == -1) {
              edgz[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[5], p[6], vals[5], vals[6]));
            }
            vertlist[5] = edgz[edgeindex];
          }
          if (edgeTable[cubeindex] & 64) {
            /* Edge 6 follows axis X with coordinates {i, j + 1, k + 1} */
            UnsignedIndex_t edgeindex =
                i + (j + 1) * (n - 1) + (k + 1) * (n - 1) * n;
            if (edgx[edgeindex] == -1) {
              edgx[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[7], p[6], vals[7], vals[6]));
            }
            vertlist[6] = edgx[edgeindex];
          }
          if (edgeTable[cubeindex] & 128) {
            /* Edge 7 follows axis Z with coordinates {i, j + 1, k} */
            UnsignedIndex_t edgeindex = i + (j + 1) * n + k * n * n;
            if (edgz[edgeindex] == -1) {
              edgz[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[4], p[7], vals[4], vals[7]));
            }
            vertlist[7] = edgz[edgeindex];
          }
          if (edgeTable[cubeindex] & 256) {
            /* Edge 8 follows axis Y with coordinates {i, j, k} */
            UnsignedIndex_t edgeindex = i + j * n + k * n * (n - 1);
            if (edgy[edgeindex] == -1) {
              edgy[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[0], p[4], vals[0], vals[4]));
            }
            vertlist[8] = edgy[edgeindex];
          }
          if (edgeTable[cubeindex] & 512) {
            /* Edge 9 follows axis Y with coordinates {i + 1, j, k} */
            UnsignedIndex_t edgeindex = i + 1 + j * n + k * n * (n - 1);
            if (edgy[edgeindex] == -1) {
              edgy[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[1], p[5], vals[1], vals[5]));
            }
            vertlist[9] = edgy[edgeindex];
          }
          if (edgeTable[cubeindex] & 1024) {
            /* Edge 10 follows axis Y with coordinates {i + 1, j, k + 1} */
            UnsignedIndex_t edgeindex = i + 1 + j * n + (k + 1) * n * (n - 1);
            if (edgy[edgeindex] == -1) {
              edgy[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[2], p[6], vals[2], vals[6]));
            }
            vertlist[10] = edgy[edgeindex];
          }
          if (edgeTable[cubeindex] & 2048) {
            /* Edge 11 follows axis Y with coordinates {i, j, k + 1} */
            UnsignedIndex_t edgeindex = i + j * n + (k + 1) * n * (n - 1);
            if (edgy[edgeindex] == -1) {
              edgy[edgeindex] = surface.nVertices();
              surface.addVertex(
                  VertexInterp(function_m, p[3], p[7], vals[3], vals[7]));
            }
            vertlist[11] = edgy[edgeindex];
          }

          /* Create the triangle */
          for (UnsignedIndex_t n = 0; triTable[cubeindex][n] != -1; n += 3) {
            /* Update vertex connectivity */
            surface.addTriangle(vertlist[triTable[cubeindex][n]],
                                vertlist[triTable[cubeindex][n + 1]],
                                vertlist[triTable[cubeindex][n + 2]]);
          }
        }
      }
    }
  }

  return surface;
}

}  // namespace IRL

#endif  // IRL_SURFACE_MESH_MARCHING_CUBES_TPP_
