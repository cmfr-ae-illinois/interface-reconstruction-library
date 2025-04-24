// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Valentin Wasquel <wasquel.valentin@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_CYLINDER_RECONSTRUCTION_CYLINDER_RECONSTRUCTION_TYPE_TRAITS_H_
#define IRL_CYLINDER_RECONSTRUCTION_CYLINDER_RECONSTRUCTION_TYPE_TRAITS_H_

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/planar_reconstruction/planar_reconstruction_type_traits.h"

namespace IRL {
//******************************************************************* //
//   Is a Cylinder based reconstruction                               //
//******************************************************************* //
template <class C>
struct is_cylinder_reconstruction : std::false_type {};

template <class C>
struct is_cylinder_reconstruction<const C> : is_cylinder_reconstruction<C> {
};

template <class ScalarType>
struct is_cylinder_reconstruction<CylinderBase<ScalarType>>
    : std::true_type {};

//******************************************************************* //
//   Contains a Cylinder based reconstruction                         //
//******************************************************************* //
template <class C>
struct has_cylinder_reconstruction : std::false_type {};

template <class C>
struct has_cylinder_reconstruction<const C>
    : has_cylinder_reconstruction<C> {};

template <class ScalarType>
struct has_cylinder_reconstruction<CylinderBase<ScalarType>>
    : std::true_type {};

template <class ScalarType>
struct has_cylinder_reconstruction<LocalizedCylinder<ScalarType>>
    : std::true_type {};

template <class ScalarType>
struct has_cylinder_reconstruction<LocalizedCylinderLink<ScalarType>>
    : std::true_type {};

//********************************************************************* //
//  Contains a localizer, base from planar_reconstruction_type_traits.h //
//********************************************************************* //
template <class ScalarType>
struct has_localizer<LocalizedCylinder<ScalarType>> : std::true_type {};

template <class ScalarType>
struct has_localizer<LocalizedCylinderLink<ScalarType>> : std::true_type {};

//******************************************************************* //
//        Linked Planar Reconstructions
//        True for anything that inherits from UnDirectedGraphNode<Self>
//******************************************************************* //
template <class ScalarType>
struct is_reconstruction_link<LocalizedCylinderLink<ScalarType>>
    : std::true_type {};

}  // namespace IRL

#endif  // IRL_CYLINDER_RECONSTRUCTION_CYLINDER_RECONSTRUCTION_TYPE_TRAITS_H_
