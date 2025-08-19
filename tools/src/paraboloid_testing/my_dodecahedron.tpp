// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_MY_DODECAHEDRON_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_MY_DODECAHEDRON_TPP_

#include <cassert>

#include "irl/geometry/general/moment_calculation_through_simplices.h"
#include "irl/geometry/half_edge_structures/half_edge.h"
#include "irl/geometry/half_edge_structures/half_edge.h"

 namespace IRL {

namespace my_dodecahedron_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 27> face_triangle_decomposition{{
{7, 1, 2}, 
{7, 2, 9}, 
{7, 9, 8}, 
{9, 2, 3}, 
{9, 3, 11}, 
{9, 11, 10}, 
{11, 3, 4}, 
{11, 4, 13}, 
{11, 13, 12}, 
{19, 15, 16}, 
{19, 16, 17}, 
{19, 17, 18}, 
{8, 15, 19}, 
{8, 19, 6}, 
{8, 6, 7}, 
{15, 8, 9}, 
{15, 9, 10}, 
{15, 10, 16}, 
{16, 10, 11}, 
{16, 11, 12}, 
{16, 12, 17}, 
{14, 5, 6}, 
{14, 6, 19}, 
{14, 19, 18}, 
{17, 12, 13}, 
{17, 13, 14}, 
{17, 14, 18}}}; 
} // namespace my_dodecahedron_triangulation 

template<class Derived, class VertexType>
HalfEdgePolyhedron<VertexType> MyDodecahedronSpecialization<Derived, VertexType>::generateHalfEdgeVersion(void) const{
HalfEdgePolyhedron<VertexType> half_edge_version;
this->setHalfEdgeVersion(&half_edge_version);
return half_edge_version;
}

template<class Derived, class VertexType>
template<class HalfEdgePolyhedronType> void MyDodecahedronSpecialization<Derived, VertexType>::setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const{
using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

a_half_edge_version->resize(60, 20, 12);

for (UnsignedIndex_t v = 0; v < 20; ++v){
a_half_edge_version->getVertex(v).setLocation((*this)[v]);
}

static constexpr std::array<UnsignedIndex_t, 60> ending_vertex_mapping{{4, 3, 2, 1, 0, 0, 1, 7, 6, 5, 1, 2, 9, 8, 7, 2, 3, 11, 10, 9, 3, 4, 13, 12, 11, 4, 0, 5, 14, 13, 15, 16, 17, 18, 19, 15, 19, 6, 7, 8, 8, 9, 10, 16, 15, 10, 11, 12, 17, 16, 5, 6, 19, 18, 14, 12, 13, 14, 18, 17}};
static constexpr std::array<UnsignedIndex_t, 60> previous_half_edge_mapping{{4, 0, 1, 2, 3, 9, 5, 6, 7, 8, 14, 10, 11, 12, 13, 19, 15, 16, 17, 18, 24, 20, 21, 22, 23, 29, 25, 26, 27, 28, 34, 30, 31, 32, 33, 39, 35, 36, 37, 38, 44, 40, 41, 42, 43, 49, 45, 46, 47, 48, 54, 50, 51, 52, 53, 59, 55, 56, 57, 58}};
static constexpr std::array<UnsignedIndex_t, 60> next_half_edge_mapping{{1, 2, 3, 4, 0, 6, 7, 8, 9, 5, 11, 12, 13, 14, 10, 16, 17, 18, 19, 15, 21, 22, 23, 24, 20, 26, 27, 28, 29, 25, 31, 32, 33, 34, 30, 36, 37, 38, 39, 35, 41, 42, 43, 44, 40, 46, 47, 48, 49, 45, 51, 52, 53, 54, 50, 56, 57, 58, 59, 55}};
static constexpr std::array<UnsignedIndex_t, 60> face_mapping{{0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11}};
static constexpr std::array<UnsignedIndex_t, 60> opposite_half_edge_mapping{{26, 21, 16, 11, 6, 27, 4, 10, 38, 51, 7, 3, 15, 41, 39, 12, 2, 20, 46, 42, 17, 1, 25, 56, 47, 22, 0, 5, 50, 57, 36, 44, 49, 59, 53, 40, 30, 52, 8, 14, 35, 13, 19, 45, 31, 43, 18, 24, 55, 32, 28, 9, 37, 34, 58, 48, 23, 29, 54, 33}};
for (UnsignedIndex_t n = 0; n < static_cast<UnsignedIndex_t>(ending_vertex_mapping.size()); ++n){
HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
current_half_edge = HalfEdgeType(&a_half_edge_version->getVertex(ending_vertex_mapping[n]), &a_half_edge_version->getHalfEdge(previous_half_edge_mapping[n]), &a_half_edge_version->getHalfEdge(next_half_edge_mapping[n]), &a_half_edge_version->getFace(face_mapping[n]));
current_half_edge.setOppositeHalfEdge(&a_half_edge_version->getHalfEdge(opposite_half_edge_mapping[n]));
current_half_edge.getFace()->setStartingHalfEdge(&current_half_edge);
current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
}

}

template<class Derived, class VertexType>
constexpr UnsignedIndex_t MyDodecahedronSpecialization<Derived, VertexType>::getNumberOfSimplicesInDecomposition(void)  {
return static_cast<UnsignedIndex_t>(my_dodecahedron_triangulation::face_triangle_decomposition.size());
}

template<class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4> MyDodecahedronSpecialization<Derived, VertexType>::getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
assert(a_tet < MyDodecahedronSpecialization::getNumberOfSimplicesInDecomposition());
return {my_dodecahedron_triangulation::face_triangle_decomposition[a_tet][0],
my_dodecahedron_triangulation::face_triangle_decomposition[a_tet][1],
my_dodecahedron_triangulation::face_triangle_decomposition[a_tet][2],
my_dodecahedron_triangulation::datum_index};
}

template<class Derived, class VertexType>
ProxyTet<Derived> MyDodecahedronSpecialization<Derived, VertexType>::getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const{
assert(a_tet < this->getNumberOfSimplicesInDecomposition());
return { static_cast<const Derived&>(*this),this->getSimplexIndicesFromDecomposition(a_tet)};
}


 } // namespace IRL
#endif //SRC_GEOMETRY_POLYHEDRONS_MY_DODECAHEDRON_TPP_
