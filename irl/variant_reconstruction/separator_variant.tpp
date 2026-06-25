// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SEPARATOR_VARIANT_TPP_
#define IRL_SEPARATOR_VARIANT_TPP_

namespace IRL {

inline SeparatorVariant::SeparatorVariant(
    const SeparatorUnion& a_separator_union) {
  (*this) = a_separator_union;
}

inline SeparatorVariant& SeparatorVariant::operator=(
    const SeparatorUnion& a_separator_union) {
  switch (a_separator_union.type()) {
    case SeparatorUnion::SeparatorType::OnePlane:
      (*this) = PlanarSeparator::fromOnePlane(a_separator_union.getPlane());
      break;

    case SeparatorUnion::SeparatorType::Paraboloid:
      (*this) = a_separator_union.getParaboloid();
      break;

    case SeparatorUnion::SeparatorType::Cylinder:
      (*this) = a_separator_union.getCylinder();
      break;

    default:
      throw std::runtime_error(
          "Unsupported SeparatorUnion type in "
          "SeparatorVariant conversion");
  }

  return *this;
}

inline void SeparatorVariant::setToPlanarSeparator(void) {
  if (not std::holds_alternative<PlanarSeparator>(*this)) {
    (*this) = PlanarSeparator();
  }
}

inline void SeparatorVariant::setToParaboloid(void) {
  if (not std::holds_alternative<Paraboloid>(*this)) {
    (*this) = Paraboloid();
  }
}

inline void SeparatorVariant::setToCylinder(void) {
  if (not std::holds_alternative<Cylinder>(*this)) {
    (*this) = Cylinder();
  }
}

inline std::pair<double, double> SeparatorVariant::getPrincipalCurvatures(
    void) {
  if (const auto separator = std::get_if<PlanarSeparator>(this)) {
    return std::make_pair(0.0, 0.0);
  } else if (const auto separator = std::get_if<Paraboloid>(this)) {
    const double a = separator->getAlignedParaboloid().a();
    const double b = separator->getAlignedParaboloid().b();
    return std::make_pair(2.0 * a, 2.0 * b);
  } else if (const auto separator = std::get_if<Cylinder>(this)) {
    const double r = separator->getAlignedCylinder().r();
    return std::make_pair(std::sqrt(r), 0.0);
  } else {
    throw std::runtime_error("Variant type cannot return principal curvatures");
  }
}

inline void SeparatorVariant::setPrincipalCurvatures(const double k1,
                                                     const double k2) {
  if (const auto separator = std::get_if<Paraboloid>(this)) {
    const auto& datum = separator->getDatum();
    const auto& frame = separator->getReferenceFrame();
    if (k1 == 0.0 && k2 == 0.0) {  // Replace by plane
      const auto plane = Plane(frame[2], frame[2] * datum);
      (*this) = PlanarSeparator::fromOnePlane(plane);
    } else {
      separator->setAlignedParaboloid(AlignedParaboloid({0.5 * k1, 0.5 * k2}));
    }
  } else if (const auto separator = std::get_if<Cylinder>(this)) {
    throw std::runtime_error(
        "Principal curvatures cannot be set for variant type");
  } else {
    throw std::runtime_error(
        "Principal curvatures cannot be set for variant type");
  }
}

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

inline std::ostream& operator<<(std::ostream& out,
                                const SeparatorVariant& a_reconstruction) {
  if (const auto separator = std::get_if<PlanarSeparator>(&a_reconstruction)) {
    out << (*separator);
  } else if (const auto separator =
                 std::get_if<Paraboloid>(&a_reconstruction)) {
    out << (*separator);
  } else if (const auto separator = std::get_if<Cylinder>(&a_reconstruction)) {
    out << (*separator);
  } else {
    throw std::runtime_error("Variant type cannot be printed");
  }
  return out;
}

}  // namespace IRL

#endif  // IRL_SEPARATOR_VARIANT_TPP_
