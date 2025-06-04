// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_H_
#define IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_H_

#include <variant>
#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/planar_reconstruction/planar_separator.h"

namespace IRL {
using SeparatorVariant = std::variant<PlanarSeparator, Paraboloid, Cylinder>;
using LocalizedSeparatorVariant =
    JoinedReconstructions<PlanarLocalizer, SeparatorVariant>;
using LocalizedSeparatorVariantLink =
    ReconstructionLink<LocalizedSeparatorVariant, UnDirectedGraphNode>;
template <>
struct is_reconstruction_link<LocalizedSeparatorVariantLink> : std::true_type {
};

}  // namespace IRL

#endif  // IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_H_
