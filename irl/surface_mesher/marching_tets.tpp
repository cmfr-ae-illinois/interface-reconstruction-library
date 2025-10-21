// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SURFACE_MESH_MARCHING_TETS_TPP_
#define IRL_SURFACE_MESH_MARCHING_TETS_TPP_

namespace IRL {

// ------------------------------------------------------------------
// Helper structs (private to this .tpp)
// ------------------------------------------------------------------

struct TetCell {
  std::array<Pt, 4> p;
  std::array<double, 4> val;
};

struct TetTriangle {
  std::array<Pt, 3> p;
};

inline Pt VertexInterp(double iso, const Pt &p1, const Pt &p2,
                       double val1, double val2) {
  double mu = (iso - val1) / (val2 - val1);
  return p1 + mu * (p2 - p1);
}


inline MarchingTets::MarchingTets(const RectangularCuboid a_domain,
                                  const ImplicitF a_function)
    : domain_m(a_domain), function_m(a_function) {}

// ------------------------------------------------------------------
// Polygonise a single tetrahedron — helper for marching tets
// ------------------------------------------------------------------
// Edge table for marching tets, since only 16 cases, can construct locally. Code constructed using Paul Bourke template
static inline int PolygoniseTri(const TetCell &g, double iso,
                                std::array<TetTriangle, 2> &tri) {
  int ntri = 0;
  int triindex = 0;

  if (g.val[0] < iso) triindex |= 1;
  if (g.val[1] < iso) triindex |= 2;
  if (g.val[2] < iso) triindex |= 4;
  if (g.val[3] < iso) triindex |= 8;

  switch (triindex) {
    case 0x00:
    case 0x0F:
      break;

    case 0x0E:
    case 0x01:
      tri[0].p[0] = VertexInterp(iso, g.p[0], g.p[1], g.val[0], g.val[1]);
      tri[0].p[1] = VertexInterp(iso, g.p[0], g.p[2], g.val[0], g.val[2]);
      tri[0].p[2] = VertexInterp(iso, g.p[0], g.p[3], g.val[0], g.val[3]);
      ntri++;
      break;

    case 0x0D:
    case 0x02:
      tri[0].p[0] = VertexInterp(iso, g.p[1], g.p[0], g.val[1], g.val[0]);
      tri[0].p[1] = VertexInterp(iso, g.p[1], g.p[3], g.val[1], g.val[3]);
      tri[0].p[2] = VertexInterp(iso, g.p[1], g.p[2], g.val[1], g.val[2]);
      ntri++;
      break;

    case 0x0C:
    case 0x03:
      tri[0].p[0] = VertexInterp(iso, g.p[0], g.p[3], g.val[0], g.val[3]);
      tri[0].p[1] = VertexInterp(iso, g.p[0], g.p[2], g.val[0], g.val[2]);
      tri[0].p[2] = VertexInterp(iso, g.p[1], g.p[3], g.val[1], g.val[3]);
      ntri++;
      tri[1].p[0] = tri[0].p[2];
      tri[1].p[1] = VertexInterp(iso, g.p[1], g.p[2], g.val[1], g.val[2]);
      tri[1].p[2] = tri[0].p[1];
      ntri++;
      break;

    case 0x0B:
    case 0x04:
      tri[0].p[0] = VertexInterp(iso, g.p[2], g.p[0], g.val[2], g.val[0]);
      tri[0].p[1] = VertexInterp(iso, g.p[2], g.p[1], g.val[2], g.val[1]);
      tri[0].p[2] = VertexInterp(iso, g.p[2], g.p[3], g.val[2], g.val[3]);
      ntri++;
      break;

    case 0x0A:
    case 0x05:
      tri[0].p[0] = VertexInterp(iso, g.p[0], g.p[1], g.val[0], g.val[1]);
      tri[0].p[1] = VertexInterp(iso, g.p[2], g.p[3], g.val[2], g.val[3]);
      tri[0].p[2] = VertexInterp(iso, g.p[0], g.p[3], g.val[0], g.val[3]);
      ntri++;
      tri[1].p[0] = tri[0].p[0];
      tri[1].p[1] = VertexInterp(iso, g.p[1], g.p[2], g.val[1], g.val[2]);
      tri[1].p[2] = tri[0].p[1];
      ntri++;
      break;

    case 0x09:
    case 0x06:
      tri[0].p[0] = VertexInterp(iso, g.p[0], g.p[1], g.val[0], g.val[1]);
      tri[0].p[1] = VertexInterp(iso, g.p[1], g.p[3], g.val[1], g.val[3]);
      tri[0].p[2] = VertexInterp(iso, g.p[2], g.p[3], g.val[2], g.val[3]);
      ntri++;
      tri[1].p[0] = tri[0].p[0];
      tri[1].p[1] = VertexInterp(iso, g.p[0], g.p[2], g.val[0], g.val[2]);
      tri[1].p[2] = tri[0].p[2];
      ntri++;
      break;

    case 0x07:
    case 0x08:
      tri[0].p[0] = VertexInterp(iso, g.p[3], g.p[0], g.val[3], g.val[0]);
      tri[0].p[1] = VertexInterp(iso, g.p[3], g.p[2], g.val[3], g.val[2]);
      tri[0].p[2] = VertexInterp(iso, g.p[3], g.p[1], g.val[3], g.val[1]);
      ntri++;
      break;
  }

  return ntri;
}

inline TriangulatedSurfaceOutput MarchingTets::triangulate(
    const UnsignedIndex_t a_subdivision) {
  TriangulatedSurfaceOutput surface;
  surface.clearAll();

  const UnsignedIndex_t n = a_subdivision > 0 ? a_subdivision + 1 : 1;
  std::vector<std::vector<std::vector<double>>> f(
      n, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));

  // Sample scalar field
  const Pt x0 = domain_m.getLowerLimits();
  const Pt x1 = domain_m.getUpperLimits();
  const Pt dx = (x1 - x0) * (1.0 / static_cast<double>(n - 1));

  for (UnsignedIndex_t i = 0; i < n; ++i) {
    for (UnsignedIndex_t j = 0; j < n; ++j) {
      for (UnsignedIndex_t k = 0; k < n; ++k) {
        const double x = static_cast<double>(i) * dx[0] + x0[0];
        const double y = static_cast<double>(j) * dx[1] + x0[1];
        const double z = static_cast<double>(k) * dx[2] + x0[2];
        f[i][j][k] = function_m(Pt(x, y, z));
      }
    }
  }

  // Tetrahedral decomposition of each cube (6 tets per cube)
  // crack-free for this vertex numbering:
const int tetrahedra[6][4] = {
  {0, 1, 3, 7},
  {0, 3, 2, 7},
  {0, 2, 6, 7},
  {0, 6, 4, 7},
  {0, 4, 5, 7},
  {0, 5, 1, 7}
};


  std::array<Pt, 8> p;
  std::array<double, 8> val;

  for (UnsignedIndex_t i = 0; i < n - 1; ++i) {
    for (UnsignedIndex_t j = 0; j < n - 1; ++j) {
      for (UnsignedIndex_t k = 0; k < n - 1; ++k) {

        // Define cube vertices and scalar values
        for (int c = 0; c < 8; ++c) {
          const int cz = (c & 1) ? 1 : 0;
          const int cy = (c & 2) ? 1 : 0;
          const int cx = (c & 4) ? 1 : 0;
          const double x = (i + cx) * dx[0] + x0[0];
          const double y = (j + cy) * dx[1] + x0[1];
          const double z = (k + cz) * dx[2] + x0[2];
          p[c] = Pt(x, y, z);
          val[c] = f[i + cx][j + cy][k + cz];
        }

        // Process 6 tetrahedra per cube
        for (const auto &tet : tetrahedra) {
          // Collect vertices for this tetrahedron
          TetCell g;
          for (int m = 0; m < 4; ++m) {
            g.p[m] = p[tet[m]];
            g.val[m] = val[tet[m]];
          }

          // Polygonize this tetrahedron
          std::array<TetTriangle, 2> tris;
            int ntris = PolygoniseTri(g, 0.0, tris);

          // Add triangles to output surface
          for (int t = 0; t < ntris; ++t) {
            const auto &tri = tris[t];
            UnsignedIndex_t v0 = surface.nVertices();
            surface.addVertex(tri.p[0]);
            surface.addVertex(tri.p[1]);
            surface.addVertex(tri.p[2]);
            surface.addTriangle(v0, v0 + 1, v0 + 2);
          }
        }
      }
    }
  }

  return surface;
}

} // namespace IRL

#endif // IRL_SURFACE_MESH_MARCHING_TETS_TPP_
