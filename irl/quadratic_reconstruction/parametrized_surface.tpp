// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_QUADRATIC_RECONSTRUCTION_PARAMETRIZED_SURFACE_TPP_
#define IRL_QUADRATIC_RECONSTRUCTION_PARAMETRIZED_SURFACE_TPP_

#include <fstream>
#include <iomanip>

// #define IRL_USE_EARCUT
// #define IRL_USE_TRIANGLE
#define IRL_USE_CGAL
// #define IRL_USE_GEOGRAM

#include "external/NumericalIntegration/NumericalIntegration.h"

namespace IRL {

template <class List>
bool noDuplicates(List& list) {
  if (list.size() == 0) {
    return true;
  }
  for (UnsignedIndex_t i = 0; i < list.size() - 1; i++) {
    for (UnsignedIndex_t j = i + 1; j < list.size(); j++) {
      if (list[i] == list[j]) {
        return false;
      }
    }
  }
  return true;
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void splitEdge(const UnsignedIndex_t edge_id, VertexList& vertices,
               TriList& triangles, EdgeList& edges, TriConnectivity& tri_neigh,
               TriEdges& tri_edges, VertConnectivity& vert_tri,
               VertexValence& vert_neigh) {
  // Start and end point
  const int v0 = edges[edge_id][0];
  const int v1 = edges[edge_id][1];
  assert(v0 != v1);
  assert(v0 >= 0 && v0 < vertices.size());
  assert(v1 >= 0 && v1 < vertices.size());

  // Find opposite triagnles and vertices
  const int tp = edges[edge_id][2];
  assert(tp >= 0);
  int vp = -1;
  int dp = -1;
  int ep = -1;
  int tpneigh = -1;
  const int tm = edges[edge_id][3];
  int vm = -1;
  int dm = -1;
  int em = -1;
  int tmneigh = -1;
  assert(tp != tm);

  if (tp >= 0) {
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] != v0 && triangles[3 * tp + d] != v1) {
        vp = triangles[3 * tp + d];
        break;
      }
    }
    for (int d = 0; d < 3; d++) {
      if ((triangles[3 * tp + d] == v1 &&
           triangles[3 * tp + (d + 1) % 3] == vp) ||
          (triangles[3 * tp + d] == vp &&
           triangles[3 * tp + (d + 1) % 3] == v1)) {
        dp = d;
        tpneigh = tri_neigh[tp][d];
        ep = tri_edges[tp][d];
        assert(ep >= 0);
        assert((edges[ep][0] == v1 && edges[ep][1] == vp) ||
               (edges[ep][1] == v1 && edges[ep][0] == vp));
        assert((edges[ep][2] == tp && edges[ep][3] == tpneigh) ||
               (edges[ep][3] == tp && edges[ep][2] == tpneigh));
        break;
      }
    }
    assert(vp >= 0);
  }
  if (tm >= 0) {
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] != v0 && triangles[3 * tm + d] != v1) {
        vm = triangles[3 * tm + d];
        break;
      }
    }
    for (int d = 0; d < 3; d++) {
      if ((triangles[3 * tm + d] == v1 &&
           triangles[3 * tm + (d + 1) % 3] == vm) ||
          (triangles[3 * tm + d] == vm &&
           triangles[3 * tm + (d + 1) % 3] == v1)) {
        dm = d;
        tmneigh = tri_neigh[tm][d];
        em = tri_edges[tm][d];
        assert(em >= 0);
        assert((edges[em][0] == v1 && edges[em][1] == vm) ||
               (edges[em][1] == v1 && edges[em][0] == vm));
        assert((edges[em][2] == tm && edges[em][3] == tmneigh) ||
               (edges[em][3] == tm && edges[em][2] == tmneigh));
        break;
      }
    }

    assert(vm >= 0);
  }
  assert(vm != vp);

  // Create new mid-point
  const int vnew = vertices.size();
  vertices.push_back(Pt(0.5 * (vertices[v0] + vertices[v1])));
  vert_tri.resize(vert_tri.size() + 1);
  vert_neigh.resize(vert_neigh.size() + 1);
  vert_tri.back().resize(0);
  vert_neigh.back().resize(0);

  // Create new triangle(s)
  int tpnew = -1;
  if (tp >= 0) {
    tpnew = triangles.size() / 3;
    tri_neigh.push_back({-1, -1, -1});
    tri_edges.push_back({-1, -1, -1});
    triangles.push_back(vnew);
    triangles.push_back(v1);
    triangles.push_back(vp);
    // std::cout << "adding triangle " << tpnew << " = " << vnew << ", " << v1
    //           << ", " << vp << std::endl;
  }
  int tmnew = -1;
  if (tm >= 0) {
    tmnew = triangles.size() / 3;
    tri_neigh.push_back({-1, -1, -1});
    tri_edges.push_back({-1, -1, -1});
    triangles.push_back(vnew);
    triangles.push_back(v1);
    triangles.push_back(vm);
    // std::cout << "adding triangle " << tmnew << " = " << vnew << ", " << v1
    //           << ", " << vm << std::endl;
  }
  assert(tpnew != tmnew);

  // Create new edges
  edges[edge_id][1] = vnew;
  const int enew = edges.size();
  if (tp >= 0) {
    edges.push_back({vnew, v1, tpnew, tmnew});
  } else {
    edges.push_back({vnew, v1, tmnew, tpnew});
  }
  int epnew = -1;
  if (tp >= 0) {
    epnew = edges.size();
    edges.push_back({vp, vnew, tp, tpnew});
  }
  int emnew = -1;
  if (tm >= 0) {
    emnew = edges.size();
    edges.push_back({vm, vnew, tm, tmnew});
  }

  // Update old and new triangles
  if (tp >= 0) {
    bool found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] == v1) {
        triangles[3 * tp + d] = vnew;
        found_vertex = true;
        break;
      }
    }
    assert(found_vertex);
    tri_neigh[tp][dp] = tpnew;
    tri_edges[tp][dp] = epnew;
    tri_neigh[tpnew][0] = tmnew;
    tri_neigh[tpnew][1] = tpneigh;
    tri_neigh[tpnew][2] = tp;
    tri_edges[tpnew][0] = enew;
    tri_edges[tpnew][1] = ep;
    tri_edges[tpnew][2] = epnew;
    if (edges[ep][2] == tp) {
      edges[ep][2] = tpnew;
    } else if (edges[ep][3] == tp) {
      edges[ep][3] = tpnew;
    } else {
      assert(false);
    }
    if (tpneigh >= 0) {
      bool found_neigh = false;
      for (int d = 0; d < 3; d++) {
        if (tri_neigh[tpneigh][d] == tp) {
          tri_neigh[tpneigh][d] = tpnew;
          found_neigh = true;
          break;
        }
      }
      assert(found_neigh);
    }
  }
  if (tm >= 0) {
    bool found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] == v1) {
        triangles[3 * tm + d] = vnew;
        found_vertex = true;
        break;
      }
    }
    assert(found_vertex);
    tri_neigh[tm][dm] = tmnew;
    tri_edges[tm][dm] = emnew;
    tri_neigh[tmnew][0] = tpnew;
    tri_neigh[tmnew][1] = tmneigh;
    tri_neigh[tmnew][2] = tm;
    tri_edges[tmnew][0] = enew;
    tri_edges[tmnew][1] = em;
    tri_edges[tmnew][2] = emnew;
    if (edges[em][2] == tm) {
      edges[em][2] = tmnew;
    } else if (edges[em][3] == tm) {
      edges[em][3] = tmnew;
    } else {
      assert(false);
    }
    if (tmneigh >= 0) {
      bool found_neigh = false;
      for (int d = 0; d < 3; d++) {
        if (tri_neigh[tmneigh][d] == tm) {
          tri_neigh[tmneigh][d] = tmnew;
          found_neigh = true;
          break;
        }
      }
      assert(found_neigh);
    }
  }

  // Update vertice neighbours
  if (vp >= 0) {
    vert_neigh[vp].push_back(vnew);
  }
  if (vm >= 0) {
    vert_neigh[vm].push_back(vnew);
  }
  bool found_neigh = false;
  for (int d = 0; d < vert_neigh[v0].size(); d++) {
    if (vert_neigh[v0][d] == v1) {
      vert_neigh[v0][d] = vnew;
      found_neigh = true;
      break;
    }
  }
  assert(found_neigh);
  found_neigh = false;
  for (int d = 0; d < vert_neigh[v1].size(); d++) {
    if (vert_neigh[v1][d] == v0) {
      vert_neigh[v1][d] = vnew;
      found_neigh = true;
      break;
    }
  }
  assert(found_neigh);
  vert_neigh[vnew].push_back(v0);
  vert_neigh[vnew].push_back(v1);
  if (vp >= 0) {
    vert_neigh[vnew].push_back(vp);
  }
  if (vm >= 0) {
    vert_neigh[vnew].push_back(vm);
  }

  // Update vertex->tri connectivity
  if (tp >= 0) {
    bool found_tri = false;
    for (int d = 0; d < vert_tri[v1].size(); d++) {
      if (vert_tri[v1][d] == tp) {
        vert_tri[v1][d] = tpnew;
        found_tri = true;
        break;
      }
    }
    assert(found_tri);
    vert_tri[vp].push_back(tpnew);
    vert_tri[vnew].push_back(tp);
    vert_tri[vnew].push_back(tpnew);
  }
  if (tm >= 0) {
    bool found_tri = false;
    for (int d = 0; d < vert_tri[v1].size(); d++) {
      if (vert_tri[v1][d] == tm) {
        vert_tri[v1][d] = tmnew;
        found_tri = true;
        break;
      }
    }
    assert(found_tri);
    vert_tri[vm].push_back(tmnew);
    vert_tri[vnew].push_back(tm);
    vert_tri[vnew].push_back(tmnew);
  }

  if (vp >= 0) {
    assert(noDuplicates(vert_tri[vp]));
    assert(noDuplicates(vert_neigh[vp]));
  }
  if (vm >= 0) {
    assert(noDuplicates(vert_tri[vm]));
    assert(noDuplicates(vert_neigh[vm]));
  }
  assert(noDuplicates(vert_tri[v0]));
  assert(noDuplicates(vert_tri[v1]));
  assert(noDuplicates(vert_tri[vnew]));
  assert(noDuplicates(vert_neigh[v0]));
  assert(noDuplicates(vert_neigh[v1]));
  assert(noDuplicates(vert_neigh[vnew]));
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void collapseEdge(const UnsignedIndex_t edge_id, VertexList& vertices,
                  TriList& triangles, EdgeList& edges,
                  TriConnectivity& tri_neigh, TriEdges& tri_edges,
                  VertConnectivity& vert_tri, VertexValence& vert_neigh) {
  // Start and end point
  const int v0 = edges[edge_id][0];
  const int v1 = edges[edge_id][1];
  assert(v0 != v1);
  assert(v0 >= 0 && v0 < vertices.size());
  assert(v1 >= 0 && v1 < vertices.size());

  // Find opposite triagnles and vertices
  const int tp = edges[edge_id][2];
  assert(tp >= 0);
  int vp = -1;
  int dp0 = -1;
  int dp1 = -1;
  int ep0 = -1;
  int ep1 = -1;
  int tpneigh0 = -1;
  int tpneigh1 = -1;
  const int tm = edges[edge_id][3];
  int vm = -1;
  int dm0 = -1;
  int dm1 = -1;
  int em0 = -1;
  int em1 = -1;
  int tmneigh0 = -1;
  int tmneigh1 = -1;
  assert(tp != tm);

  if (tp >= 0) {
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] != v0 && triangles[3 * tp + d] != v1) {
        vp = triangles[3 * tp + d];
        break;
      }
    }
    for (int d = 0; d < 3; d++) {
      if ((triangles[3 * tp + d] == v1 &&
           triangles[3 * tp + (d + 1) % 3] == vp) ||
          (triangles[3 * tp + d] == vp &&
           triangles[3 * tp + (d + 1) % 3] == v1)) {
        dp1 = d;
        tpneigh1 = tri_neigh[tp][d];
        ep1 = tri_edges[tp][d];
        assert(ep1 >= 0);
        assert((edges[ep1][2] == tp && edges[ep1][3] == tpneigh1) ||
               (edges[ep1][3] == tp && edges[ep1][2] == tpneigh1));
      } else if ((triangles[3 * tp + d] == v0 &&
                  triangles[3 * tp + (d + 1) % 3] == vp) ||
                 (triangles[3 * tp + d] == vp &&
                  triangles[3 * tp + (d + 1) % 3] == v0)) {
        dp0 = d;
        tpneigh0 = tri_neigh[tp][d];
        ep0 = tri_edges[tp][d];
        assert(ep0 >= 0);
        assert((edges[ep0][2] == tp && edges[ep0][3] == tpneigh0) ||
               (edges[ep0][3] == tp && edges[ep0][2] == tpneigh0));
      }
    }
    assert(vp >= 0);
  }
  if (tm >= 0) {
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] != v0 && triangles[3 * tm + d] != v1) {
        vm = triangles[3 * tm + d];
        break;
      }
    }
    for (int d = 0; d < 3; d++) {
      if ((triangles[3 * tm + d] == v1 &&
           triangles[3 * tm + (d + 1) % 3] == vm) ||
          (triangles[3 * tm + d] == vm &&
           triangles[3 * tm + (d + 1) % 3] == v1)) {
        dm1 = d;
        tmneigh1 = tri_neigh[tm][d];
        em1 = tri_edges[tm][d];
        assert(em1 >= 0);
        assert((edges[em1][2] == tm && edges[em1][3] == tmneigh1) ||
               (edges[em1][3] == tm && edges[em1][2] == tmneigh1));
      } else if ((triangles[3 * tm + d] == v0 &&
                  triangles[3 * tm + (d + 1) % 3] == vm) ||
                 (triangles[3 * tm + d] == vm &&
                  triangles[3 * tm + (d + 1) % 3] == v0)) {
        dm0 = d;
        tmneigh0 = tri_neigh[tm][d];
        em0 = tri_edges[tm][d];
        assert(em0 >= 0);
        assert((edges[em0][2] == tm && edges[em0][3] == tmneigh0) ||
               (edges[em0][3] == tm && edges[em0][2] == tmneigh0));
      }
    }
    assert(vm >= 0);
  }
  assert(vm != vp);

  // Update start triangles
  for (int i = 0; i < vert_tri[v0].size(); i++) {
    const int tri = vert_tri[v0][i];
    bool found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tri + d] == v0) {
        triangles[3 * tri + d] = v1;
        found_vertex = true;
      }
      const int edge = tri_edges[tri][d];
      if (edges[edge][0] == v0) {
        edges[edge][0] = v1;
      } else if (edges[edge][1] == v0) {
        edges[edge][1] = v1;
      }
    }
    assert(found_vertex);
    bool tri_not_here = true;
    for (int j = 0; j < vert_tri[v1].size(); j++) {
      if (vert_tri[v1][j] == tri) {
        tri_not_here = false;
        break;
      }
    }
    if (tri_not_here) {
      // if (tri != tp && tri != tm) {
      vert_tri[v1].push_back(tri);
    }
  }
  vert_tri[v0].resize(0);

  for (int i = 0; i < vert_neigh[v0].size(); i++) {
    const int neigh = vert_neigh[v0][i];
    if (neigh != v1) {
      bool neigh_has_v1 = false;
      int index_v0_neigh = -1;
      int count_v0 = 0;
      for (int j = 0; j < vert_neigh[neigh].size(); j++) {
        if (vert_neigh[neigh][j] == v1) {
          neigh_has_v1 = true;
        } else if (vert_neigh[neigh][j] == v0) {
          index_v0_neigh = j;
          count_v0++;
        }
      }
      assert(count_v0 == 1);
      assert(index_v0_neigh >= 0);
      if (!neigh_has_v1) {
        vert_neigh[neigh][index_v0_neigh] = v1;
        bool v1_had_neigh = false;
        for (int k = 0; k < vert_neigh[v1].size(); k++) {
          if (vert_neigh[v1][k] == neigh) {
            v1_had_neigh = true;
            break;
          }
        }
        assert(!v1_had_neigh);
        vert_neigh[v1].push_back(neigh);
      } else {
        vert_neigh[neigh].erase(vert_neigh[neigh].begin() + index_v0_neigh);
      }
    } else {
      bool found_vertex = false;
      for (int j = 0; j < vert_neigh[v1].size(); j++) {
        if (vert_neigh[v1][j] == v0) {
          vert_neigh[v1].erase(vert_neigh[v1].begin() + j);
          found_vertex = true;
          break;
        }
      }
      assert(found_vertex);
    }
  }
  vert_neigh[v0].resize(0);

  if (tpneigh0 >= 0) {
    bool found_edge = false;
    bool found_neigh = false;
    for (int d = 0; d < 3; d++) {
      if (tri_edges[tpneigh0][d] == ep0) {
        tri_edges[tpneigh0][d] = ep1;
        found_edge = true;
      }
      if (tri_neigh[tpneigh0][d] == tp) {
        tri_neigh[tpneigh0][d] = tpneigh1;
        found_neigh = true;
      }
    }
    assert(found_edge);
    assert(found_neigh);
  }
  if (tpneigh1 >= 0) {
    bool found_neigh = false;
    for (int d = 0; d < 3; d++) {
      if (tri_neigh[tpneigh1][d] == tp) {
        tri_neigh[tpneigh1][d] = tpneigh0;
        found_neigh = true;
        break;
      }
    }
    assert(found_neigh);
  }
  if (tmneigh0 >= 0) {
    bool found_edge = false;
    bool found_neigh = false;
    for (int d = 0; d < 3; d++) {
      if (tri_edges[tmneigh0][d] == em0) {
        found_edge = true;
        tri_edges[tmneigh0][d] = em1;
      }
      if (tri_neigh[tmneigh0][d] == tm) {
        tri_neigh[tmneigh0][d] = tmneigh1;
        found_neigh = true;
      }
    }
    assert(found_edge);
    assert(found_neigh);
  }
  if (tmneigh1 >= 0) {
    bool found_neigh = false;
    for (int d = 0; d < 3; d++) {
      if (tri_neigh[tmneigh1][d] == tm) {
        tri_neigh[tmneigh1][d] = tmneigh0;
        found_neigh = true;
        break;
      }
    }
    assert(found_neigh);
  }
  if (ep1 >= 0) {
    if (edges[ep1][2] == tp) {
      assert(edges[ep1][3] == tpneigh1);
      edges[ep1][2] = tpneigh0;
    } else if (edges[ep1][3] == tp) {
      assert(edges[ep1][2] == tpneigh1);
      edges[ep1][3] = tpneigh0;
    } else {
      assert(false);
    }
  }
  if (em1 >= 0) {
    if (edges[em1][2] == tm) {
      assert(edges[em1][3] == tmneigh1);
      edges[em1][2] = tmneigh0;
    } else if (edges[em1][3] == tm) {
      assert(edges[em1][2] == tmneigh1);
      edges[em1][3] = tmneigh0;
    } else {
      assert(false);
    }
  }

  // Remove triangles adjacent to edge
  if (tp >= 0) {
    for (int d = 0; d < 3; d++) {
      triangles[3 * tp + d] = -1;
    }
  }
  if (tm >= 0) {
    for (int d = 0; d < 3; d++) {
      triangles[3 * tm + d] = -1;
    }
  }
  if (tp >= 0) {
    bool found = false;
    for (int i = 0; i < vert_tri[v1].size(); i++) {
      if (vert_tri[v1][i] == tp) {
        vert_tri[v1].erase(vert_tri[v1].begin() + i);
        found = true;
        break;
      }
    }
    assert(found);
    found = false;
    for (int i = 0; i < vert_tri[vp].size(); i++) {
      if (vert_tri[vp][i] == tp) {
        vert_tri[vp].erase(vert_tri[vp].begin() + i);
        found = true;
        break;
      }
    }
    assert(found);
  }
  if (tm >= 0) {
    bool found = false;
    for (int i = 0; i < vert_tri[v1].size(); i++) {
      if (vert_tri[v1][i] == tm) {
        vert_tri[v1].erase(vert_tri[v1].begin() + i);
        found = true;
        break;
      }
    }
    assert(found);
    found = false;
    for (int i = 0; i < vert_tri[vm].size(); i++) {
      if (vert_tri[vm][i] == tm) {
        vert_tri[vm].erase(vert_tri[vm].begin() + i);
        found = true;
        break;
      }
    }
    assert(found);
  }

  // Remove edge
  edges[edge_id] = {-1, -1, -1, -1};
  if (ep0 >= 0) {
    edges[ep0] = {-1, -1, -1, -1};
  }
  if (em0 >= 0) {
    edges[em0] = {-1, -1, -1, -1};
  }
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void flipEdge(const UnsignedIndex_t edge_id, VertexList& vertices,
              TriList& triangles, EdgeList& edges, TriConnectivity& tri_neigh,
              TriEdges& tri_edges, VertConnectivity& vert_tri,
              VertexValence& vert_neigh) {
  // Start and end point
  const int v0 = edges[edge_id][0];
  const int v1 = edges[edge_id][1];
  assert(v0 != v1);
  assert(v0 >= 0 && v0 < vertices.size());
  assert(v1 >= 0 && v1 < vertices.size());

  // Find opposite triagnles and vertices
  const int tp = edges[edge_id][2];
  assert(tp >= 0);
  int vp = -1;
  int dp0 = -1;
  int dp1 = -1;
  int ep0 = -1;
  int ep1 = -1;
  int tpneigh0 = -1;
  int tpneigh1 = -1;
  const int tm = edges[edge_id][3];
  assert(tm >= 0);
  int vm = -1;
  int dm0 = -1;
  int dm1 = -1;
  int em0 = -1;
  int em1 = -1;
  int tmneigh0 = -1;
  int tmneigh1 = -1;
  assert(tp != tm);

  // std::cout << "Flipping egde " << edge_id << " = " << v0 << " " << v1 << " "
  //           << tp << " " << tm << std::endl;
  for (int d = 0; d < 3; d++) {
    if (triangles[3 * tp + d] != v0 && triangles[3 * tp + d] != v1) {
      vp = triangles[3 * tp + d];
      break;
    }
  }
  for (int d = 0; d < 3; d++) {
    if ((triangles[3 * tp + d] == v1 &&
         triangles[3 * tp + (d + 1) % 3] == vp) ||
        (triangles[3 * tp + d] == vp &&
         triangles[3 * tp + (d + 1) % 3] == v1)) {
      dp1 = d;
      tpneigh1 = tri_neigh[tp][d];
      ep1 = tri_edges[tp][d];
      assert(ep1 >= 0);
      assert((edges[ep1][0] == v1 && edges[ep1][1] == vp) ||
             (edges[ep1][1] == v1 && edges[ep1][0] == vp));
      assert((edges[ep1][2] == tp && edges[ep1][3] == tpneigh1) ||
             (edges[ep1][3] == tp && edges[ep1][2] == tpneigh1));
    } else if ((triangles[3 * tp + d] == v0 &&
                triangles[3 * tp + (d + 1) % 3] == vp) ||
               (triangles[3 * tp + d] == vp &&
                triangles[3 * tp + (d + 1) % 3] == v0)) {
      dp0 = d;
      tpneigh0 = tri_neigh[tp][d];
      ep0 = tri_edges[tp][d];
      assert(ep0 >= 0);
      assert((edges[ep0][0] == v0 && edges[ep0][1] == vp) ||
             (edges[ep0][1] == v0 && edges[ep0][0] == vp));
      assert((edges[ep0][2] == tp && edges[ep0][3] == tpneigh0) ||
             (edges[ep0][3] == tp && edges[ep0][2] == tpneigh0));
    }
  }
  assert(vp >= 0);
  for (int d = 0; d < 3; d++) {
    if (triangles[3 * tm + d] != v0 && triangles[3 * tm + d] != v1) {
      vm = triangles[3 * tm + d];
      break;
    }
  }
  for (int d = 0; d < 3; d++) {
    if ((triangles[3 * tm + d] == v1 &&
         triangles[3 * tm + (d + 1) % 3] == vm) ||
        (triangles[3 * tm + d] == vm &&
         triangles[3 * tm + (d + 1) % 3] == v1)) {
      dm1 = d;
      tmneigh1 = tri_neigh[tm][d];
      em1 = tri_edges[tm][d];
      assert(em1 >= 0);
      assert((edges[em1][0] == v1 && edges[em1][1] == vm) ||
             (edges[em1][1] == v1 && edges[em1][0] == vm));
      assert((edges[em1][2] == tm && edges[em1][3] == tmneigh1) ||
             (edges[em1][3] == tm && edges[em1][2] == tmneigh1));
    } else if ((triangles[3 * tm + d] == v0 &&
                triangles[3 * tm + (d + 1) % 3] == vm) ||
               (triangles[3 * tm + d] == vm &&
                triangles[3 * tm + (d + 1) % 3] == v0)) {
      dm0 = d;
      tmneigh0 = tri_neigh[tm][d];
      em0 = tri_edges[tm][d];
      assert(em0 >= 0);
      assert((edges[em0][0] == v0 && edges[em0][1] == vm) ||
             (edges[em0][1] == v0 && edges[em0][0] == vm));
      assert((edges[em0][2] == tm && edges[em0][3] == tmneigh0) ||
             (edges[em0][3] == tm && edges[em0][2] == tmneigh0));
    }
  }
  assert(vm >= 0);
  assert(vm != vp);

  bool do_flip = true;
  if (tpneigh0 >= 0 && tmneigh0 >= 0) {
    if (tpneigh0 == tmneigh0) {
      do_flip = false;
    }
  }
  if (tpneigh1 >= 0 && tmneigh1 >= 0) {
    if (tpneigh1 == tmneigh1) {
      do_flip = false;
    }
  }

  if (do_flip) {
    // Flip edge
    edges[edge_id][0] = vp;
    edges[edge_id][1] = vm;

    bool found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] == v0) {
        triangles[3 * tp + (d + 1) % 3] = vm;
        triangles[3 * tp + (d + 2) % 3] = vp;
        tri_neigh[tp][d] = tmneigh0;
        tri_neigh[tp][(d + 1) % 3] = tm;
        tri_neigh[tp][(d + 2) % 3] = tpneigh0;
        tri_edges[tp][d] = em0;
        tri_edges[tp][(d + 1) % 3] = edge_id;
        tri_edges[tp][(d + 2) % 3] = ep0;
        found_vertex = true;
        break;
      }
    }
    found_vertex = false;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] == v1) {
        triangles[3 * tm + (d + 1) % 3] = vm;
        triangles[3 * tm + (d + 2) % 3] = vp;
        tri_neigh[tm][d] = tmneigh1;
        tri_neigh[tm][(d + 1) % 3] = tp;
        tri_neigh[tm][(d + 2) % 3] = tpneigh1;
        tri_edges[tm][d] = em1;
        tri_edges[tm][(d + 1) % 3] = edge_id;
        tri_edges[tm][(d + 2) % 3] = ep1;
        found_vertex = true;
        break;
      }
    }
    assert(found_vertex);
    if (tpneigh1 >= 0) {
      bool found_tri = false;
      for (int d = 0; d < 3; d++) {
        if (tri_neigh[tpneigh1][d] == tp) {
          tri_neigh[tpneigh1][d] = tm;
          found_tri = true;
          break;
        }
      }
      assert(found_tri);
    }
    assert(found_vertex);
    if (tmneigh0 >= 0) {
      bool found_tri = false;
      for (int d = 0; d < 3; d++) {
        if (tri_neigh[tmneigh0][d] == tm) {
          tri_neigh[tmneigh0][d] = tp;
          found_tri = true;
          break;
        }
      }
      assert(found_tri);
    }
    if (edges[ep1][2] == tp) {
      edges[ep1][2] = tm;
    } else if (edges[ep1][3] == tp) {
      edges[ep1][3] = tm;
    } else {
      assert(false);
    }
    if (edges[em0][2] == tm) {
      edges[em0][2] = tp;
    } else if (edges[em0][3] == tm) {
      edges[em0][3] = tp;
    } else {
      assert(false);
    }

    vert_tri[vp].push_back(tm);
    vert_tri[vm].push_back(tp);
    for (int i = 0; i < vert_tri[v0].size(); i++) {
      if (vert_tri[v0][i] == tm) {
        vert_tri[v0].erase(vert_tri[v0].begin() + i);
        break;
      }
    }
    for (int i = 0; i < vert_tri[v1].size(); i++) {
      if (vert_tri[v1][i] == tp) {
        vert_tri[v1].erase(vert_tri[v1].begin() + i);
        break;
      }
    }

    bool has_vm = false;
    for (int i = 0; i < vert_neigh[vp].size(); i++) {
      if (vert_neigh[vp][i] == vm) {
        has_vm = true;
        break;
      }
    }
    if (!has_vm) {
      vert_neigh[vp].push_back(vm);
    }
    bool has_vp = false;
    for (int i = 0; i < vert_neigh[vm].size(); i++) {
      if (vert_neigh[vm][i] == vp) {
        has_vp = true;
        break;
      }
    }
    if (!has_vp) {
      vert_neigh[vm].push_back(vp);
    }
    for (int i = 0; i < vert_neigh[v0].size(); i++) {
      if (vert_neigh[v0][i] == v1) {
        vert_neigh[v0].erase(vert_neigh[v0].begin() + i);
        break;
      }
    }
    for (int i = 0; i < vert_neigh[v1].size(); i++) {
      if (vert_neigh[v1][i] == v0) {
        vert_neigh[v1].erase(vert_neigh[v1].begin() + i);
        break;
      }
    }

    assert(noDuplicates(vert_tri[vp]));
    assert(noDuplicates(vert_tri[vm]));
    assert(noDuplicates(vert_tri[v0]));
    assert(noDuplicates(vert_tri[v1]));
    assert(noDuplicates(vert_neigh[vp]));
    assert(noDuplicates(vert_neigh[vm]));
    assert(noDuplicates(vert_neigh[v0]));
    assert(noDuplicates(vert_neigh[v1]));
    assert((edges[em0][2] == tmneigh0 && edges[em0][3] == tp) ||
           (edges[em0][3] == tmneigh0 && edges[em0][2] == tp));
    assert((edges[ep0][2] == tpneigh0 && edges[ep0][3] == tp) ||
           (edges[ep0][3] == tpneigh0 && edges[ep0][2] == tp));
    assert((edges[em1][2] == tmneigh1 && edges[em1][3] == tm) ||
           (edges[em1][3] == tmneigh1 && edges[em1][2] == tm));
    assert((edges[ep1][2] == tpneigh1 && edges[ep1][3] == tm) ||
           (edges[ep1][3] == tpneigh1 && edges[ep1][2] == tm));
  }
}

template <class VertexList, class EdgeList, class VertexValence>
void relaxVertices(VertexList& vertices, EdgeList& edges,
                   VertexValence& vert_neigh,
                   const UnsignedIndex_t fixed_vertices,
                   const UnsignedIndex_t iterations,
                   const double relax_factor) {
  if (vertices.size() > fixed_vertices) {
    std::vector<bool> fixed;
    fixed.resize(vertices.size(), false);
    for (UnsignedIndex_t i = 0; i < edges.size(); i++) {
      if ((edges[i][0] >= 0 || edges[i][1] >= 0) &&
          (edges[i][2] < 0 || edges[i][3] < 0)) {
        fixed[edges[i][0]] = true;
        fixed[edges[i][1]] = true;
      }
    }
    assert(vert_neigh.size() == vertices.size());
    std::vector<Pt> barycenters;
    barycenters.resize(vertices.size() - fixed_vertices);
    for (UnsignedIndex_t it = 0; it < iterations; it++) {
      for (UnsignedIndex_t i = fixed_vertices; i < vertices.size(); i++) {
        barycenters[i - fixed_vertices] = Pt(0, 0, 0);
        if (!fixed[i] && vert_neigh[i].size() > 0) {
          assert(noDuplicates(vert_neigh[i]));
          for (UnsignedIndex_t j = 0; j < vert_neigh[i].size(); j++) {
            assert(vert_neigh[i][j] >= 0);
            assert(vert_neigh[i][j] < vertices.size());
            barycenters[i - fixed_vertices] += vertices[vert_neigh[i][j]];
          }
          barycenters[i - fixed_vertices] /=
              static_cast<double>(vert_neigh[i].size());
        }
      }
      for (UnsignedIndex_t i = fixed_vertices; i < vertices.size(); i++) {
        if (!fixed[i] && vert_neigh[i].size() > 0) {
          vertices[i] = (1.0 - relax_factor) * vertices[i] +
                        relax_factor * barycenters[i - fixed_vertices];
        }
      }
    }
  }
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void splitLongEdges(VertexList& vertices, TriList& triangles, EdgeList& edges,
                    TriConnectivity& tri_neigh, TriEdges& tri_edges,
                    VertConnectivity& vert_tri, VertexValence& vert_neigh,
                    const double high) {
  const double high_sq = high * high;
  bool remains_long_edges = true;
  UnsignedIndex_t it = 0;
  while (it < 1 && remains_long_edges) {
    it++;
    const int n_edges = edges.size();
    remains_long_edges = false;
    for (int i = 0; i < n_edges; ++i) {
      const int v0 = edges[i][0];
      const int v1 = edges[i][1];
      if (v0 < 0 || v1 < 0) {
        continue;
      }
      const double edge_length_sq =
          squaredMagnitude(vertices[v0] - vertices[v1]);
      if (edge_length_sq > high_sq) {
        remains_long_edges = true;
        splitEdge(i, vertices, triangles, edges, tri_neigh, tri_edges, vert_tri,
                  vert_neigh);
      }
    }
  }
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void collapseShortEdges(VertexList& vertices, TriList& triangles,
                        EdgeList& edges, TriConnectivity& tri_neigh,
                        TriEdges& tri_edges, VertConnectivity& vert_tri,
                        VertexValence& vert_neigh, const double low,
                        const double high,
                        const UnsignedIndex_t fixed_vertices) {
  const double low_sq = low * low;
  const double high_sq = high * high;
  bool remains_short_edges = true;
  UnsignedIndex_t it = 0;
  while (it < 1 && remains_short_edges) {
    it++;
    const int n_edges = edges.size();
    remains_short_edges = false;
    for (int i = 0; i < n_edges; ++i) {
      const int v0 = edges[i][0];
      const int v1 = edges[i][1];
      if (v0 < fixed_vertices || v1 < 0 || vert_tri[v0].size() < 5 ||
          vert_tri[v1].size() < 5) {
        continue;
      }
      const double edge_length_sq =
          squaredMagnitude(vertices[v0] - vertices[v1]);
      if (edge_length_sq < low_sq) {
        bool collapse = true;
        for (int j = 0; j < vert_neigh[v0].size(); j++) {
          const int neigh = vert_neigh[v0][j];
          const double test_edge_length_sq =
              squaredMagnitude(vertices[neigh] - vertices[v1]);
          if (test_edge_length_sq > high_sq) {
            collapse = false;
            break;
          }
        }
        if (collapse) {
          std::cout << "Collapsing edge " << i << " " << v0 << " -- " << v1
                    << " -- " << edges[i][2] << " -- " << edges[i][3]
                    << std::endl;
          remains_short_edges = true;
          collapseEdge(i, vertices, triangles, edges, tri_neigh, tri_edges,
                       vert_tri, vert_neigh);
        }
        // else {
        //   break;
        // }
      }
    }
  }
}

template <class VertexList, class TriList, class EdgeList,
          class TriConnectivity, class TriEdges, class VertConnectivity,
          class VertexValence>
void equalizeValence(VertexList& vertices, TriList& triangles, EdgeList& edges,
                     TriConnectivity& tri_neigh, TriEdges& tri_edges,
                     VertConnectivity& vert_tri, VertexValence& vert_neigh,
                     const UnsignedIndex_t fixed_vertices) {
  const int n_edges = edges.size();
  for (int i = 0; i < n_edges; ++i) {
    const int v0 = edges[i][0];
    const int v1 = edges[i][1];
    const int tp = edges[i][2];
    const int tm = edges[i][3];
    if (v0 < 0 || v1 < 0 || vert_neigh[v0].size() == 0 ||
        vert_neigh[v1].size() == 0 || tp < 0 || tm < 0) {
      continue;
    }
    int vp = -1;
    int vm = -1;
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tp + d] != v0 && triangles[3 * tp + d] != v1) {
        vp = triangles[3 * tp + d];
        break;
      }
    }
    assert(vp >= 0);
    for (int d = 0; d < 3; d++) {
      if (triangles[3 * tm + d] != v0 && triangles[3 * tm + d] != v1) {
        vm = triangles[3 * tm + d];
        break;
      }
    }
    assert(vm >= 0);
    int val0 = vert_neigh[v0].size();
    int val1 = vert_neigh[v1].size();
    int valp = vert_neigh[vp].size();
    int valm = vert_neigh[vm].size();
    int target0 = v0 < fixed_vertices ? 4 : 6;
    int target1 = v1 < fixed_vertices ? 4 : 6;
    int targetp = vp < fixed_vertices ? 4 : 6;
    int targetm = vm < fixed_vertices ? 4 : 6;
    int deviation_pre = std::abs(val0 - target0) + std::abs(val1 - target1) +
                        std::abs(valp - targetp) + std::abs(valm - targetm);
    flipEdge(i, vertices, triangles, edges, tri_neigh, tri_edges, vert_tri,
             vert_neigh);
    val0 = vert_neigh[v0].size();
    val1 = vert_neigh[v1].size();
    valp = vert_neigh[vp].size();
    valm = vert_neigh[vm].size();
    int deviation_post = std::abs(val0 - target0) + std::abs(val1 - target1) +
                         std::abs(valp - targetp) + std::abs(valm - targetm);
    if (deviation_pre <= deviation_post) {
      flipEdge(i, vertices, triangles, edges, tri_neigh, tri_edges, vert_tri,
               vert_neigh);
    }
  }
}

template <class MomentType, class SurfaceType>
MomentType& AddSurfaceOutput<MomentType, SurfaceType>::getMoments(void) {
  return volume_moments_m;
}

template <class MomentType, class SurfaceType>
const MomentType& AddSurfaceOutput<MomentType, SurfaceType>::getMoments(
    void) const {
  return volume_moments_m;
}

template <class MomentType, class SurfaceType>
SurfaceType& AddSurfaceOutput<MomentType, SurfaceType>::getSurface(void) {
  return surface_m;
}

template <class MomentType, class SurfaceType>
const SurfaceType& AddSurfaceOutput<MomentType, SurfaceType>::getSurface(
    void) const {
  return surface_m;
}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput()
    : knows_surface_area_m{false},
      knows_avg_normal_m{false},
      knows_int_mean_curv_m{false},
      knows_int_gaussian_curv_m{false},
      length_scale_m{-1.0} {}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    ParametrizedSurfaceOutput&& a_rhs)
    : arc_list_m(std::move(a_rhs.arc_list_m)),
      knows_surface_area_m(a_rhs.knows_surface_area_m),
      surface_area_m(a_rhs.surface_area_m),
      knows_avg_normal_m(a_rhs.knows_avg_normal_m),
      avg_normal_m(a_rhs.avg_normal_m),
      knows_int_mean_curv_m(a_rhs.knows_int_mean_curv_m),
      int_mean_curv_m(a_rhs.int_mean_curv_m),
      knows_int_gaussian_curv_m(a_rhs.knows_int_gaussian_curv_m),
      int_gaussian_curv_m(a_rhs.int_gaussian_curv_m),
      length_scale_m{a_rhs.length_scale_m},
      pt_from_bezier_split_m(std::move(a_rhs.pt_from_bezier_split_m)) {}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    const ParametrizedSurfaceOutput& a_rhs)
    : arc_list_m(a_rhs.arc_list_m),
      knows_surface_area_m(a_rhs.knows_surface_area_m),
      surface_area_m(a_rhs.surface_area_m),
      knows_avg_normal_m(a_rhs.knows_avg_normal_m),
      avg_normal_m(a_rhs.avg_normal_m),
      knows_int_mean_curv_m(a_rhs.knows_int_mean_curv_m),
      int_mean_curv_m(a_rhs.int_mean_curv_m),
      knows_int_gaussian_curv_m(a_rhs.knows_int_gaussian_curv_m),
      int_gaussian_curv_m(a_rhs.int_gaussian_curv_m),
      length_scale_m{a_rhs.length_scale_m} {}

inline ParametrizedSurfaceOutput& ParametrizedSurfaceOutput::operator=(
    ParametrizedSurfaceOutput&& a_rhs) {
  if (this != &a_rhs) {
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

inline ParametrizedSurfaceOutput& ParametrizedSurfaceOutput::operator=(
    const ParametrizedSurfaceOutput& a_rhs) {
  if (this != &a_rhs) {
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

inline void ParametrizedSurfaceOutput::setLengthScale(
    const double a_length_scale) {
  length_scale_m = a_length_scale;
}

inline RationalBezierArc& ParametrizedSurfaceOutput::operator[](
    const UnsignedIndex_t a_index) {
  return arc_list_m[a_index];
}

inline const RationalBezierArc& ParametrizedSurfaceOutput::operator[](
    const UnsignedIndex_t a_index) const {
  return arc_list_m[a_index];
}

inline std::vector<RationalBezierArc>& ParametrizedSurfaceOutput::getArcs(
    void) {
  return arc_list_m;
}

inline std::vector<Pt*>& ParametrizedSurfaceOutput::getPts(void) {
  return pt_from_bezier_split_m;
}

inline void ParametrizedSurfaceOutput::addArc(
    const RationalBezierArc& a_rational_bezier_arc) {
  #ifdef VALDEBUG2
  std::cout << "adding the following arc to the output : " << a_rational_bezier_arc << std::endl;
  #endif
  arc_list_m.push_back(a_rational_bezier_arc);
}

inline void ParametrizedSurfaceOutput::addPt(Pt* a_pt) {
  pt_from_bezier_split_m.push_back(a_pt);
}

inline const std::vector<RationalBezierArc>::size_type
ParametrizedSurfaceOutput::size(void) const {
  return arc_list_m.size();
}

inline void ParametrizedSurfaceOutput::clearArcs(void) { arc_list_m.clear(); }

inline void ParametrizedSurfaceOutput::clearPts(void) {
  for (auto& elem : pt_from_bezier_split_m) {
    delete elem;
  }
  pt_from_bezier_split_m.clear();
}

inline void ParametrizedSurfaceOutput::clear(void) {
  this->clearArcs();
  this->clearPts();
}

inline ParametrizedSurfaceOutput::~ParametrizedSurfaceOutput(void) {
  for (auto elem : pt_from_bezier_split_m) {
    delete elem;
  }
}

inline double ParametrizedSurfaceOutput::getAverageMeanCurvature(void) {
  return this->getMeanCurvatureIntegral() / safelyTiny(this->getSurfaceArea());
}

inline double ParametrizedSurfaceOutput::getAverageGaussianCurvature(void) {
  return this->getGaussianCurvatureIntegral() /
         safelyTiny(this->getSurfaceArea());
}

inline std::ostream& operator<<(
    std::ostream& out,
    const ParametrizedSurfaceOutput& a_parametrized_surface) {
  out.precision(16);
  for (UnsignedIndex_t i = 0; i < a_parametrized_surface.size(); ++i) {
    out << a_parametrized_surface[i];
    if (i < a_parametrized_surface.size() - 1) out << std::endl;
  }
  return out;
}

}  // namespace IRL

#endif  // IRL_QUADRATIC_RECONSTRUCTION_PARAMETRIZED_SURFACE_TPP_
