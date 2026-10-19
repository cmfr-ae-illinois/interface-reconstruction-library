// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_AREA_TPP_
#define IRL_MOMENTS_AREA_TPP_

namespace IRL {

template <class ScalarType>
inline AreaBase<ScalarType>::AreaBase(void) : area_m{ScalarType(0)} {}

template <class ScalarType>
inline constexpr AreaBase<ScalarType>::AreaBase(const ScalarType a_value)
    : area_m{a_value} {}

template <class ScalarType>
inline constexpr AreaBase<ScalarType>::AreaBase(const AreaBase<double>& a_value)
    : area_m{ScalarType(a_value.area())} {}

template <class ScalarType>
inline constexpr AreaBase<ScalarType>::AreaBase(const AreaBase<Quad_t>& a_value)
    : area_m{ScalarType(a_value.area())} {}

template <class ScalarType>
inline constexpr AreaBase<ScalarType>::AreaBase(
    const VolumeBase<ScalarType>& a_value)
    : area_m{ScalarType(a_value.volume())} {}

template <class ScalarType>
inline AreaBase<ScalarType> AreaBase<ScalarType>::fromScalarConstant(
    const ScalarType a_value) {
  return AreaBase<ScalarType>(a_value);
}

template <class ScalarType>
template <class GeometryType>
inline AreaBase<ScalarType> AreaBase<ScalarType>::calculateMoments(
    GeometryType* a_geometry) {
  return a_geometry->calculateVolume();
}

template <class ScalarType>
inline ScalarType& AreaBase<ScalarType>::area(void) {
  return area_m;
}

template <class ScalarType>
inline const ScalarType& AreaBase<ScalarType>::area(void) const {
  return area_m;
}

template <class ScalarType>
inline AreaBase<ScalarType>& AreaBase<ScalarType>::operator+=(
    const AreaBase<ScalarType>& a_rhs) {
  area_m += a_rhs.area_m;
  return (*this);
}

template <class ScalarType>
inline AreaBase<ScalarType>& AreaBase<ScalarType>::operator*=(
    const ScalarType a_rhs) {
  area_m *= a_rhs;
  return (*this);
}

template <class ScalarType>
inline AreaBase<ScalarType>& AreaBase<ScalarType>::operator/=(
    const ScalarType a_rhs) {
  area_m /= a_rhs;
  return (*this);
}

template <class ScalarType>
inline AreaBase<ScalarType>::operator ScalarType() const {
  return area_m;
}

template <class ScalarType>
inline AreaBase<ScalarType>& AreaBase<ScalarType>::operator=(
    const ScalarType a_value) {
  area_m = a_value;
  return (*this);
}

template <class ScalarType>
inline AreaBase<ScalarType>& AreaBase<ScalarType>::operator=(
    const VolumeBase<ScalarType>& a_value) {
  area_m = a_value;
  return (*this);
}

template <class ScalarType>
inline AreaBase<ScalarType> AreaBase<ScalarType>::operator-(void) const {
  return AreaBase<ScalarType>(-area_m);
}

template <class ScalarType>
inline std::ostream& operator<<(std::ostream& out,
                                const AreaBase<ScalarType>& a_area) {
  out << static_cast<double>(a_area);
  return out;
}

template <class ScalarType>
inline AreaBase<ScalarType> operator+(const AreaBase<ScalarType>& a_am1,
                                      const AreaBase<ScalarType>& a_am2) {
  return AreaBase<ScalarType>(a_am1.area() + a_am2.area());
}

template <class ScalarType>
inline AreaBase<ScalarType> operator-(const AreaBase<ScalarType>& a_am1,
                                      const AreaBase<ScalarType>& a_am2) {
  return AreaBase<ScalarType>(a_am1.area() - a_am2.area());
}

template <class ScalarType>
inline AreaBase<ScalarType> operator*(const ScalarType a_multiplier,
                                      const AreaBase<ScalarType>& a_am) {
  return AreaBase<ScalarType>(a_multiplier * a_am.area());
}

template <class ScalarType>
inline AreaBase<ScalarType> operator*(const AreaBase<ScalarType>& a_am,
                                      const ScalarType a_multiplier) {
  return a_multiplier * a_am;
}

}  // namespace IRL

#endif  // IRL_MOMENTS_AREA_TPP_
