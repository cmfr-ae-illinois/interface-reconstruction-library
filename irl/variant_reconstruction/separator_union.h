// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SEPARATOR_UNION_H_
#define IRL_SEPARATOR_UNION_H_

#include <cstdint>

#include "irl/generic_cutting/general/class_classifications.h"
#include "irl/parameters/defined_types.h"

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/geometry/general/plane.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace IRL {

// This class is meant to be trivially copyable for use in AMReX
// For application outside AMReX, use SeparatorVariant instead
class SeparatorUnion {
 public:
  enum class SeparatorType : std::uint8_t {
    OnePlane = 1,
    TwoPlanes = 2,
    Paraboloid = 3,
    Cylinder = 4
  };

  SeparatorUnion() {
    type_m = SeparatorType::OnePlane;
    planes_m[0] = Plane();
  };

  SeparatorUnion(const Plane& a_plane) {
    type_m = SeparatorType::OnePlane;
    planes_m[0] = a_plane;
  };

  SeparatorUnion(const Plane& a_plane_0, const Plane& a_plane_1) {
    type_m = SeparatorType::TwoPlanes;
    planes_m[0] = a_plane_0;
    planes_m[1] = a_plane_1;
  };

  SeparatorUnion(const Paraboloid& a_paraboloid) {
    type_m = SeparatorType::Paraboloid;
    paraboloid_m = a_paraboloid;
  };

  SeparatorUnion(const Cylinder& a_cylinder) {
    type_m = SeparatorType::Cylinder;
    cylinder_m = a_cylinder;
  };

  SeparatorUnion& operator=(const SeparatorUnion& other) {
    type_m = other.type();
    switch (type_m) {
      case SeparatorType::OnePlane:
        planes_m[0] = other.getPlanes()[0];
        break;
      case SeparatorType::TwoPlanes:
        planes_m = other.getPlanes();
        break;
      case SeparatorType::Paraboloid:
        paraboloid_m = other.getParaboloid();
        break;
      case SeparatorType::Cylinder:
        cylinder_m = other.getCylinder();
        break;
      default:
        std::runtime_error(
            "Unrecognized reconstruction type in SeparatorUnion");
    }
    return *this;
  };

  // This just replaces the current separator union
  // Only defined because it is required by AMReX
  SeparatorUnion& operator+=(const SeparatorUnion& other) {
    type_m = other.type();
    switch (type_m) {
      case SeparatorType::OnePlane:
        planes_m[0] = other.getPlanes()[0];
        break;
      case SeparatorType::TwoPlanes:
        planes_m = other.getPlanes();
        break;
      case SeparatorType::Paraboloid:
        paraboloid_m = other.getParaboloid();
        break;
      case SeparatorType::Cylinder:
        cylinder_m = other.getCylinder();
        break;
      default:
        std::runtime_error(
            "Unrecognized reconstruction type in SeparatorUnion");
    }
    return *this;
  };

  SeparatorUnion& operator=(const Plane& other) {
    type_m = SeparatorType::OnePlane;
    planes_m[0] = other;
    return *this;
  };

  SeparatorUnion& operator=(const PlanarSeparator& other) {
    switch (other.getNumberOfPlanes()) {
      case 0:
        type_m = SeparatorType::OnePlane;
        planes_m[0] = Plane();
        break;
      case 1:
        type_m = SeparatorType::OnePlane;
        planes_m[0] = other[0];
        break;
      case 2:
        type_m = SeparatorType::TwoPlanes;
        planes_m[0] = other[0];
        planes_m[1] = other[0];
        break;
      default:
        std::runtime_error("SeparatorUnion cannot contain more than 2 planes");
    }
    return *this;
  };

  SeparatorUnion& operator=(const Paraboloid& other) {
    type_m = SeparatorType::Paraboloid;
    paraboloid_m = other;
    return *this;
  };

  SeparatorUnion& operator=(const Cylinder& other) {
    type_m = SeparatorType::Cylinder;
    cylinder_m = other;
    return *this;
  };

  const SeparatorType& type(void) const { return type_m; };
  SeparatorType& type(void) { return type_m; };
  const Plane& getPlane(void) const { return planes_m[0]; };
  Plane& getPlane(void) { return planes_m[0]; };
  const Plane& getPlane(const UnsignedIndex_t a_index) const {
    return planes_m[a_index];
  };
  Plane& getPlane(const UnsignedIndex_t a_index) { return planes_m[a_index]; };
  const std::array<Plane, 2>& getPlanes(void) const { return planes_m; };
  std::array<Plane, 2>& getPlanes(void) { return planes_m; };
  const Paraboloid& getParaboloid(void) const { return paraboloid_m; };
  Paraboloid& getParaboloid(void) { return paraboloid_m; };
  const Cylinder& getCylinder(void) const { return cylinder_m; };
  Cylinder& getCylinder(void) { return cylinder_m; };

  void setToPlane(void) {
    type_m = SeparatorType::OnePlane;
    planes_m[0] = Plane();
  };
  void setToPlane(const Plane& a_plane) {
    type_m = SeparatorType::OnePlane;
    planes_m[0] = a_plane;
  };
  void setToPlanes(const Plane& a_plane0, const Plane& a_plane1) {
    type_m = SeparatorType::TwoPlanes;
    planes_m[0] = a_plane0;
    planes_m[1] = a_plane1;
  };
  void setToParaboloid(void) { type_m = SeparatorType::Paraboloid; };
  void setToParaboloid(const Paraboloid& a_paraboloid) {
    type_m = SeparatorType::Paraboloid;
    paraboloid_m = a_paraboloid;
  };
  void setToCylinder(void) { type_m = SeparatorType::Cylinder; };
  void setToCylinder(const Cylinder& a_cylinder) {
    type_m = SeparatorType::Cylinder;
    cylinder_m = a_cylinder;
  };

  void serialize(ByteBuffer* a_buffer) const {
    a_buffer->pack(&this->type(), 1);
    switch (this->type()) {
      case SeparatorType::OnePlane:
        this->getPlane().serialize(a_buffer);
        break;
      case SeparatorType::TwoPlanes:
        this->getPlanes()[0].serialize(a_buffer);
        this->getPlanes()[1].serialize(a_buffer);
        break;
      case SeparatorType::Paraboloid:
        this->getParaboloid().serialize(a_buffer);
        break;
      case SeparatorType::Cylinder:
        this->getCylinder().serialize(a_buffer);
        break;
      default:
        std::runtime_error("SeparatorUnion type cannot be serialized");
    }
  };

  void unpackSerialized(ByteBuffer* a_buffer) {
    a_buffer->unpack(&this->type(), 1);
    switch (this->type()) {
      case SeparatorType::OnePlane:
        this->getPlane().unpackSerialized(a_buffer);
        break;
      case SeparatorType::TwoPlanes:
        this->getPlanes()[0].unpackSerialized(a_buffer);
        this->getPlanes()[1].unpackSerialized(a_buffer);
        break;
      case SeparatorType::Paraboloid:
        this->getParaboloid().unpackSerialized(a_buffer);
        break;
      case SeparatorType::Cylinder:
        this->getCylinder().unpackSerialized(a_buffer);
        break;
      default:
        std::runtime_error("SeparatorUnion type cannot be unpacked");
    }
  };

  void shift(const Pt a_shift) {
    switch (type_m) {
      case SeparatorType::OnePlane:
        planes_m[0].distance() += planes_m[0].normal() * a_shift;
        break;
      case SeparatorType::TwoPlanes:
        planes_m[0].distance() += planes_m[0].normal() * a_shift;
        planes_m[1].distance() += planes_m[1].normal() * a_shift;
        break;
      case SeparatorType::Paraboloid:
        paraboloid_m.setDatum(paraboloid_m.getDatum() + a_shift);
        break;
      case SeparatorType::Cylinder:
        cylinder_m.setDatum(cylinder_m.getDatum() + a_shift);
        break;
      default:
        std::runtime_error("SeparatorUnion type cannot shift datum");
    }
  };

 private:
  SeparatorType type_m;
  union {
    std::array<Plane, 2> planes_m;
    Paraboloid paraboloid_m;
    Cylinder cylinder_m;
  };
};

inline std::ostream& operator<<(std::ostream& out,
                                const SeparatorUnion& a_reconstruction);

inline SeparatorUnion operator*(const SeparatorUnion& a_sep1,
                                const SeparatorUnion& a_sep2) {
  return a_sep1;
}

inline SeparatorUnion operator*(const double a_double,
                                const SeparatorUnion& a_sep) {
  return a_sep;
}

inline SeparatorUnion operator*(const SeparatorUnion& a_sep,
                                const double a_double) {
  return a_double * a_sep;
}

using LocalizedSeparatorUnion =
    JoinedReconstructions<PlanarLocalizer, SeparatorUnion>;
using LocalizedSeparatorUnionLink =
    ReconstructionLink<LocalizedSeparatorUnion, UnDirectedGraphNode>;
template <>
struct is_reconstruction_link<LocalizedSeparatorUnionLink> : std::true_type {};

}  // namespace IRL

#include "irl/variant_reconstruction/separator_union.tpp"

#endif  // IRL_SEPARATOR_UNION_H_
