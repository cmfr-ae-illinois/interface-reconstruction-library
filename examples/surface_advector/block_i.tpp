// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_BLOCK_I_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_BLOCK_I_TPP_

#include <cassert>

#include "irl/geometry/general/moment_calculation_through_simplices.h"
#include "irl/geometry/half_edge_structures/half_edge.h"
#include "irl/geometry/half_edge_structures/half_edge.h"

 namespace IRL {

namespace block_i_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 30> face_triangle_decomposition{{
{12, 23, 22}, 
{12, 22, 21}, 
{12, 21, 20}, 
{12, 20, 19}, 
{12, 19, 18}, 
{12, 18, 17}, 
{12, 17, 16}, 
{12, 16, 15}, 
{12, 15, 14}, 
{12, 14, 13}, 
{1, 13, 14}, 
{1, 14, 2}, 
{2, 14, 15}, 
{2, 15, 3}, 
{3, 15, 16}, 
{3, 16, 4}, 
{4, 16, 17}, 
{4, 17, 5}, 
{5, 17, 18}, 
{5, 18, 6}, 
{6, 18, 19}, 
{6, 19, 7}, 
{7, 19, 20}, 
{7, 20, 8}, 
{8, 20, 21}, 
{8, 21, 9}, 
{9, 21, 22}, 
{9, 22, 10}, 
{10, 22, 23}, 
{10, 23, 11}}}; 
} // namespace block_i_triangulation 

template<class Derived, class VertexType>
HalfEdgePolyhedron<VertexType> BlockISpecialization<Derived, VertexType>::generateHalfEdgeVersion(void) const{
HalfEdgePolyhedron<VertexType> half_edge_version;
this->setHalfEdgeVersion(&half_edge_version);
return half_edge_version;
}

template<class Derived, class VertexType>
template<class HalfEdgePolyhedronType> void BlockISpecialization<Derived, VertexType>::setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const{
using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

a_half_edge_version->resize(72, 24, 14);

for (UnsignedIndex_t v = 0; v < 24; ++v){
a_half_edge_version->getVertex(v).setLocation((*this)[v]);
}

static constexpr std::array<UnsignedIndex_t, 72> ending_vertex_mapping{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 12, 13, 1, 0, 13, 14, 2, 1, 14, 15, 3, 2, 15, 16, 4, 3, 16, 17, 5, 4, 17, 18, 6, 5, 18, 19, 7, 6, 19, 20, 8, 7, 20, 21, 9, 8, 21, 22, 10, 9, 22, 23, 11, 10, 23, 12, 0, 11}};
static constexpr std::array<UnsignedIndex_t, 72> previous_half_edge_mapping{{11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 23, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 27, 24, 25, 26, 31, 28, 29, 30, 35, 32, 33, 34, 39, 36, 37, 38, 43, 40, 41, 42, 47, 44, 45, 46, 51, 48, 49, 50, 55, 52, 53, 54, 59, 56, 57, 58, 63, 60, 61, 62, 67, 64, 65, 66, 71, 68, 69, 70}};
static constexpr std::array<UnsignedIndex_t, 72> next_half_edge_mapping{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 12, 25, 26, 27, 24, 29, 30, 31, 28, 33, 34, 35, 32, 37, 38, 39, 36, 41, 42, 43, 40, 45, 46, 47, 44, 49, 50, 51, 48, 53, 54, 55, 52, 57, 58, 59, 56, 61, 62, 63, 60, 65, 66, 67, 64, 69, 70, 71, 68}};
static constexpr std::array<UnsignedIndex_t, 72> face_mapping{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13}};
static constexpr std::array<UnsignedIndex_t, 72> opposite_half_edge_mapping{{27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 69, 65, 61, 57, 53, 49, 45, 41, 37, 33, 29, 25, 70, 23, 28, 0, 26, 22, 32, 1, 30, 21, 36, 2, 34, 20, 40, 3, 38, 19, 44, 4, 42, 18, 48, 5, 46, 17, 52, 6, 50, 16, 56, 7, 54, 15, 60, 8, 58, 14, 64, 9, 62, 13, 68, 10, 66, 12, 24, 11}};
for (UnsignedIndex_t n = 0; n < static_cast<UnsignedIndex_t>(ending_vertex_mapping.size()); ++n){
HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
current_half_edge = HalfEdgeType(&a_half_edge_version->getVertex(ending_vertex_mapping[n]), &a_half_edge_version->getHalfEdge(previous_half_edge_mapping[n]), &a_half_edge_version->getHalfEdge(next_half_edge_mapping[n]), &a_half_edge_version->getFace(face_mapping[n]));
current_half_edge.setOppositeHalfEdge(&a_half_edge_version->getHalfEdge(opposite_half_edge_mapping[n]));
current_half_edge.getFace()->setStartingHalfEdge(&current_half_edge);
current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
}

}

template<class Derived, class VertexType>
constexpr UnsignedIndex_t BlockISpecialization<Derived, VertexType>::getNumberOfSimplicesInDecomposition(void)  {
return static_cast<UnsignedIndex_t>(block_i_triangulation::face_triangle_decomposition.size());
}

template<class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4> BlockISpecialization<Derived, VertexType>::getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
assert(a_tet < BlockISpecialization::getNumberOfSimplicesInDecomposition());
return {block_i_triangulation::face_triangle_decomposition[a_tet][0],
block_i_triangulation::face_triangle_decomposition[a_tet][1],
block_i_triangulation::face_triangle_decomposition[a_tet][2],
block_i_triangulation::datum_index};
}

template<class Derived, class VertexType>
ProxyTet<Derived> BlockISpecialization<Derived, VertexType>::getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const{
assert(a_tet < this->getNumberOfSimplicesInDecomposition());
return { static_cast<const Derived&>(*this),this->getSimplexIndicesFromDecomposition(a_tet)};
}


 } // namespace IRL
#endif //SRC_GEOMETRY_POLYHEDRONS_BLOCK_I_TPP_
