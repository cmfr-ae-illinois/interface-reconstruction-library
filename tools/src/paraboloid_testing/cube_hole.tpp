// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_CUBE_HOLE_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_CUBE_HOLE_TPP_

#include <cassert>

#include "irl/geometry/general/moment_calculation_through_simplices.h"
#include "irl/geometry/half_edge_structures/half_edge.h"
#include "irl/geometry/half_edge_structures/half_edge.h"

 namespace IRL {

namespace cube_hole_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 20> face_triangle_decomposition{{
{4, 5, 6}, 
{4, 6, 7}, 
{6, 2, 3}, 
{6, 3, 7}, 
{9, 8, 1}, 
{9, 1, 2}, 
{9, 2, 6}, 
{9, 6, 10}, 
{11, 10, 6}, 
{11, 6, 5}, 
{11, 5, 1}, 
{11, 1, 8}, 
{8, 9, 13}, 
{8, 13, 12}, 
{9, 10, 14}, 
{9, 14, 13}, 
{10, 11, 15}, 
{10, 15, 14}, 
{11, 8, 12}, 
{11, 12, 15}}}; 
} // namespace cube_hole_triangulation 

template<class Derived, class VertexType>
HalfEdgePolyhedron<VertexType> CubeHoleSpecialization<Derived, VertexType>::generateHalfEdgeVersion(void) const{
HalfEdgePolyhedron<VertexType> half_edge_version;
this->setHalfEdgeVersion(&half_edge_version);
return half_edge_version;
}

template<class Derived, class VertexType>
template<class HalfEdgePolyhedronType> void CubeHoleSpecialization<Derived, VertexType>::setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const{
using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

a_half_edge_version->resize(56, 16, 12);

for (UnsignedIndex_t v = 0; v < 16; ++v){
a_half_edge_version->getVertex(v).setLocation((*this)[v]);
}

static constexpr std::array<UnsignedIndex_t, 56> ending_vertex_mapping{{2, 1, 0, 3, 1, 5, 4, 0, 5, 6, 7, 4, 2, 3, 7, 6, 8, 1, 2, 6, 10, 9, 10, 6, 5, 1, 8, 11, 12, 0, 4, 7, 14, 15, 14, 7, 3, 0, 12, 13, 9, 13, 12, 8, 10, 14, 13, 9, 11, 15, 14, 10, 8, 12, 15, 11}};
static constexpr std::array<UnsignedIndex_t, 56> previous_half_edge_mapping{{3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10, 15, 12, 13, 14, 21, 16, 17, 18, 19, 20, 27, 22, 23, 24, 25, 26, 33, 28, 29, 30, 31, 32, 39, 34, 35, 36, 37, 38, 43, 40, 41, 42, 47, 44, 45, 46, 51, 48, 49, 50, 55, 52, 53, 54}};
static constexpr std::array<UnsignedIndex_t, 56> next_half_edge_mapping{{1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 17, 18, 19, 20, 21, 16, 23, 24, 25, 26, 27, 22, 29, 30, 31, 32, 33, 28, 35, 36, 37, 38, 39, 34, 41, 42, 43, 40, 45, 46, 47, 44, 49, 50, 51, 48, 53, 54, 55, 52}};
static constexpr std::array<UnsignedIndex_t, 56> face_mapping{{0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11}};
static constexpr std::array<UnsignedIndex_t, 56> opposite_half_edge_mapping{{13, 18, 4, 37, 2, 25, 8, 30, 6, 24, 15, 31, 19, 0, 36, 10, 40, 26, 1, 12, 23, 44, 48, 20, 9, 5, 17, 52, 54, 38, 7, 11, 35, 50, 46, 32, 14, 3, 29, 42, 16, 47, 39, 53, 21, 51, 34, 41, 22, 55, 33, 45, 27, 43, 28, 49}};
for (UnsignedIndex_t n = 0; n < static_cast<UnsignedIndex_t>(ending_vertex_mapping.size()); ++n){
HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
current_half_edge = HalfEdgeType(&a_half_edge_version->getVertex(ending_vertex_mapping[n]), &a_half_edge_version->getHalfEdge(previous_half_edge_mapping[n]), &a_half_edge_version->getHalfEdge(next_half_edge_mapping[n]), &a_half_edge_version->getFace(face_mapping[n]));
current_half_edge.setOppositeHalfEdge(&a_half_edge_version->getHalfEdge(opposite_half_edge_mapping[n]));
current_half_edge.getFace()->setStartingHalfEdge(&current_half_edge);
current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
}

}

template<class Derived, class VertexType>
constexpr UnsignedIndex_t CubeHoleSpecialization<Derived, VertexType>::getNumberOfSimplicesInDecomposition(void)  {
return static_cast<UnsignedIndex_t>(cube_hole_triangulation::face_triangle_decomposition.size());
}

template<class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4> CubeHoleSpecialization<Derived, VertexType>::getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
assert(a_tet < CubeHoleSpecialization::getNumberOfSimplicesInDecomposition());
return {cube_hole_triangulation::face_triangle_decomposition[a_tet][0],
cube_hole_triangulation::face_triangle_decomposition[a_tet][1],
cube_hole_triangulation::face_triangle_decomposition[a_tet][2],
cube_hole_triangulation::datum_index};
}

template<class Derived, class VertexType>
ProxyTet<Derived> CubeHoleSpecialization<Derived, VertexType>::getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const{
assert(a_tet < this->getNumberOfSimplicesInDecomposition());
return { static_cast<const Derived&>(*this),this->getSimplexIndicesFromDecomposition(a_tet)};
}


 } // namespace IRL
#endif //SRC_GEOMETRY_POLYHEDRONS_CUBE_HOLE_TPP_
