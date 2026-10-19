// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_AREA_H_
#define IRL_MOMENTS_AREA_H_

#include <ostream>

#include "irl/moments/volume.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief An area class which is just a double
/// with special properties. Allows more general
/// writing of functions.
template <class ScalarType>
class AreaBase {
 public:
  using value_type = ScalarType;

  /// \brief Default constructor.
  AreaBase(void);

  /// \brief Constructor that initializes area to value,
  /// want to allow implicit casting of double to Area here.
  constexpr AreaBase(const ScalarType a_value);

  /// \brief Overload assignment to convert precision.
  constexpr AreaBase(const AreaBase<double>& a_area);
  constexpr AreaBase(const AreaBase<Quad_t>& a_area);
  constexpr AreaBase(const VolumeBase<ScalarType>& a_area);

  static AreaBase fromScalarConstant(const ScalarType a_value);

  /// \brief Obtain Area from the supplied geometry.
  template <class GeometryType>
  static AreaBase calculateMoments(GeometryType* a_geometry);

  /// \brief Overload += operator to update area.
  AreaBase& operator+=(const AreaBase& a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  AreaBase& operator*=(const ScalarType a_rhs);

  /// \brief Overload /= operator to divide by constant double
  AreaBase& operator/=(const ScalarType a_rhs);

  /// \brief Allow implicit conversion to double.
  operator ScalarType() const;

  /// \brief Return value of stored area.
  ScalarType& area(void);
  /// \brief Return const reference to stored area.
  const ScalarType& area(void) const;

  /// \brief Overload assignment to assign constant value to moments.
  AreaBase& operator=(const ScalarType a_value);

  /// \brief Overload assignment
  AreaBase& operator=(const VolumeBase<ScalarType>& a_value);

  // Unary - operator
  AreaBase operator-(void) const;

  /// \brief Default destructor.
  ~AreaBase(void) = default;

 private:
  ScalarType area_m;  ///< \brief Area of something
};

template <class ScalarType>
std::ostream& operator<<(std::ostream& out, const AreaBase<ScalarType>& a_area);

/// \brief Overload + operator to add two geometric moments together
template <class ScalarType>
inline AreaBase<ScalarType> operator+(const AreaBase<ScalarType>& a_am1,
                                      const AreaBase<ScalarType>& a_am2);
/// \brief Overload - operator to subtract one
/// geometric moment object from another.
template <class ScalarType>
inline AreaBase<ScalarType> operator-(const AreaBase<ScalarType>& a_am1,
                                      const AreaBase<ScalarType>& a_am2);
/// \brief Overload * operator to multiply area/centroid
template <class ScalarType>
inline AreaBase<ScalarType> operator*(const ScalarType a_multiplier,
                                      const AreaBase<ScalarType>& a_am);
/// \brief Overload * operator to multiply area/centroid
template <class ScalarType>
inline AreaBase<ScalarType> operator*(const AreaBase<ScalarType>& a_am,
                                      const ScalarType a_multiplier);

using Area = AreaBase<double>;

}  // namespace IRL

#include "irl/moments/area.tpp"

#endif  // IRL_MOMENTS_AREA_H_ */
