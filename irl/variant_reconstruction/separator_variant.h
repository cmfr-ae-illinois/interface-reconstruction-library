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
#include <tuple>


#include "irl/generic_cutting/general/class_classifications.h"

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/planar_reconstruction/planar_separator.h"

namespace IRL {

class SeparatorVariant
    : public std::variant<PlanarSeparator, Paraboloid, Cylinder> {
 public:
  using base = std::variant<PlanarSeparator, Paraboloid, Cylinder>;
  using base::base;
  using base::operator=;

  void setToPlanarSeparator(void);
  void setToParaboloid(void);
  void setToCylinder(void);

  void serialize(ByteBuffer* a_buffer) const;
  void unpackSerialized(ByteBuffer* a_buffer);

  std::tuple<double,Eigen::Vector3d,Eigen::Matrix3d>
    getSignedDistanceAndGradAndHessianSep(const Pt& a_pt, const Pt& a_centroid) const;
};

inline std::ostream& operator<<(std::ostream& out,
                                const SeparatorVariant& a_reconstruction);

using LocalizedSeparatorVariant =
    JoinedReconstructions<PlanarLocalizer, SeparatorVariant>;
using LocalizedSeparatorVariantLink =
    ReconstructionLink<LocalizedSeparatorVariant, UnDirectedGraphNode>;
template <>
struct is_reconstruction_link<LocalizedSeparatorVariantLink> : std::true_type {
};

}  // namespace IRL

#include "irl/variant_reconstruction/separator_variant.tpp"

#endif  // IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_H_
