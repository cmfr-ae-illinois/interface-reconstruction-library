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
class alignas(16) SeparatorUnion {
 public:
  enum class SeparatorType : std::uint8_t {
    OnePlane = 1,
    TwoPlanes = 2,
    Paraboloid = 3,
    Cylinder = 4
  };

  // Default constructor
  SeparatorUnion(void);

  // Constructors from 1/2 planes, 1 paraboloid, 1 cylinder
  SeparatorUnion(const Plane& a_plane);
  SeparatorUnion(const Plane& a_plane_0, const Plane& a_plane_1);
  SeparatorUnion(const Paraboloid& a_paraboloid);
  SeparatorUnion(const Cylinder& a_cylinder);

  // Overloaded copy
  SeparatorUnion& operator=(const SeparatorUnion& other);

  // This just replaces the current separator union
  // Only defined because it is required by AMReX
  SeparatorUnion& operator+=(const SeparatorUnion& other);

  // Overloaded Plane copy (sets type to 1 plane)
  SeparatorUnion& operator=(const Plane& other);

  // Overloaded PlanarSeparator copy (sets type to 1 or 2 planes)
  SeparatorUnion& operator=(const PlanarSeparator& other);

  // Overloaded Paraboloid copy (sets type to paraboloid)
  SeparatorUnion& operator=(const Paraboloid& other);

  // Overloaded Cylinder copy (sets type to cylinder)
  SeparatorUnion& operator=(const Cylinder& other);

  // Returns type of stored separator
  const SeparatorType& type(void) const;
  SeparatorType& type(void);

  // Returns first stored plane
  const Plane& getPlane(void) const;
  Plane& getPlane(void);

  // Returns stored plane corresponding to index
  const Plane& getPlane(const UnsignedIndex_t a_index) const;
  Plane& getPlane(const UnsignedIndex_t a_index);

  // Returns stored array containing 2 planes
  const std::array<Plane, 2>& getPlanes(void) const;
  std::array<Plane, 2>& getPlanes(void);

  // Returns stored Paraboloid
  const Paraboloid& getParaboloid(void) const;
  Paraboloid& getParaboloid(void);

  // Returns stored Cylinder
  const Cylinder& getCylinder(void) const;
  Cylinder& getCylinder(void);

  // Sets plane with index (only does something if type is 1/2 planes)
  void setPlane(const int a_index, const Plane& a_plane);

  // Sets type to 1 plane
  void setToPlane(void);

  // Sets type to 1 plane and initializes plane
  void setToPlane(const Plane& a_plane);

  // Sets type to 2 planes and initializes planes
  void setToPlanes(const Plane& a_plane0, const Plane& a_plane1);

  // Sets type to Paraboloid
  void setToParaboloid(void);

  // Sets type to Paraboloid and initializes object
  void setToParaboloid(const Paraboloid& a_paraboloid);

  // Sets type to Cylinder
  void setToCylinder(void);

  // Sets type to Cylinder and initializes object
  void setToCylinder(const Cylinder& a_cylinder);

  // Sets to full Paraboloid
  void setToFull(void);

  // Sets to empty Paraboloid
  void setToEmpty(void);

  // Returns boolean for full/empy shortcuts
  const bool isFull(void) const;
  const bool isEmpty(void) const;

  // Serializes object to byte buffer
  void serialize(ByteBuffer* a_buffer) const;

  // Unpacks object from byte buffer
  void unpackSerialized(ByteBuffer* a_buffer);

  // Shifts origin of geometric object in 3D
  void shift(const Pt a_shift);

  // Reflects object about location in x, y, or z direction
  void reflect(const SeparatorUnion& a_ref, const int a_dir,
               const double a_loc);

 private:
  SeparatorType type_m;
  union {
    std::array<Plane, 2> planes_m;
    Paraboloid paraboloid_m;
    Cylinder cylinder_m;
  };
};

static_assert(
    sizeof(SeparatorUnion) == 128,
    "SeparatorUnion size must be 128 for its corresponding Fortran type "
    "to be correct");

inline std::ostream& operator<<(std::ostream& out,
                                const SeparatorUnion& a_reconstruction);

inline SeparatorUnion operator*(const SeparatorUnion& a_sep1,
                                const SeparatorUnion& a_sep2);

inline SeparatorUnion operator*(const double a_double,
                                const SeparatorUnion& a_sep);

inline SeparatorUnion operator*(const SeparatorUnion& a_sep,
                                const double a_double);

using LocalizedSeparatorUnion =
    JoinedReconstructions<PlanarLocalizer, SeparatorUnion>;
using LocalizedSeparatorUnionLink =
    ReconstructionLink<LocalizedSeparatorUnion, UnDirectedGraphNode>;
template <>
struct is_reconstruction_link<LocalizedSeparatorUnionLink> : std::true_type {};

}  // namespace IRL

#include "irl/variant_reconstruction/separator_union.tpp"

#endif  // IRL_SEPARATOR_UNION_H_
