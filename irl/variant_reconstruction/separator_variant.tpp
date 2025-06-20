// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_TPP_
#define IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_TPP_

namespace IRL {

inline void SeparatorVariant::serialize(ByteBuffer* a_buffer) const {
  const std::size_t index = this->index();
  a_buffer->pack(&index, 1);
  if (const auto separator = std::get_if<PlanarSeparator>(this)) {
    separator->serialize(a_buffer);
  } else if (const auto separator = std::get_if<Paraboloid>(this)) {
    separator->serialize(a_buffer);
  } else if (const auto separator = std::get_if<Cylinder>(this)) {
    separator->serialize(a_buffer);
  } else {
    throw std::runtime_error("Variant type cannot be serialized");
  }
}

inline void SeparatorVariant::unpackSerialized(ByteBuffer* a_buffer) {
  std::size_t index;
  a_buffer->unpack(&index, 1);
  if (index == 0) {
    PlanarSeparator separator;
    separator.unpackSerialized(a_buffer);
    (*this) = separator;
  } else if (index == 1) {
    Paraboloid separator;
    separator.unpackSerialized(a_buffer);
    (*this) = separator;
  } else if (index == 2) {
    Cylinder separator;
    separator.unpackSerialized(a_buffer);
    (*this) = separator;
  } else {
    throw std::runtime_error("Variant type cannot be unpacked");
  }
}

}  // namespace IRL

#endif  // IRL_PLANAR_RECONSTRUCTION_SEPARATOR_VARIANT_TPP_
