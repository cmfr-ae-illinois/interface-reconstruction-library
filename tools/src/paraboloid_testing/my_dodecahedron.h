// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_MY_DODECAHEDRON_H_
#define SRC_GEOMETRY_POLYHEDRONS_MY_DODECAHEDRON_H_

#include "irl/geometry/general/geometry_type_traits.h"
#include "irl/geometry/general/stored_vertex_access.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/tet.h"
#include "irl/parameters/defined_types.h"

namespace IRL{

template <class Derived, class VertexType>
class MyDodecahedronSpecialization: public BasePolyhedron<Derived, VertexType, ProxyTet<Derived>> {
public: 

HalfEdgePolyhedron<VertexType> generateHalfEdgeVersion(void) const; 

template<class HalfEdgePolyhedronType>
void setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const;

static constexpr UnsignedIndex_t getNumberOfSimplicesInDecomposition(void);

static constexpr std::array<UnsignedIndex_t, 4> getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet);

ProxyTet<Derived> getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const;

}; 

template <class VertexType>
class StoredMyDodecahedron: public StoredVertexAccess<StoredMyDodecahedron<VertexType>,VertexType,20>, public MyDodecahedronSpecialization<StoredMyDodecahedron<VertexType>, VertexType>{
friend StoredVertexAccess<StoredMyDodecahedron<VertexType>,VertexType,20>;

public: 

using StoredVertexAccess<StoredMyDodecahedron<VertexType>,VertexType,20>::StoredVertexAccess;

StoredMyDodecahedron(void) = default;
}; 


// Predefined types 
using MyDodecahedron= StoredMyDodecahedron<Pt>; 


template<class VertexType>
struct is_polyhedron<StoredMyDodecahedron<VertexType>> : std::true_type {};
} // namespace IRL 

#include "my_dodecahedron.tpp"
#endif //SRC_GEOMETRY_POLYHEDRONS_MY_DODECAHEDRON_H_
