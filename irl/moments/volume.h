// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_VOLUME_H_
#define IRL_MOMENTS_VOLUME_H_

#include <ostream>

#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief A volume class which is just a double
/// with special properties. Allows more general
/// writing of functions.
template <class ScalarType>
class VolumeBase {
 public:
  using gradient_type = VolumeBase<ScalarType>;
  using value_type = ScalarType;

  /// \brief Default constructor.
  VolumeBase(void);

  /// \brief Constructor that initializes volume to value,
  /// want to allow implicit casting of double to Volume here.
  constexpr VolumeBase(const ScalarType a_value);

  /// \brief Overload assignment to convert precision.
  constexpr VolumeBase(const VolumeBase<double>& a_volume);
  constexpr VolumeBase(const VolumeBase<Quad_t>& a_volume);

  static VolumeBase fromScalarConstant(const ScalarType a_value);

  /// \brief Obtain Volume from the supplied geometry.
  template <class GeometryType>
  static VolumeBase calculateMoments(GeometryType* a_geometry);

  /// \brief Dummy function to allow general use along
  /// with VolumeMoments.
  void multiplyByVolume(void);

  /// \brief Dummy function to allow general use along
  /// with VolumeMoments.
  void normalizeByVolume(void);

  /// \brief Overload += operator to update volume.
  VolumeBase& operator+=(const VolumeBase& a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  VolumeBase& operator*=(const ScalarType a_rhs);

  /// \brief Overload /= operator to divide by constant double
  VolumeBase& operator/=(const ScalarType a_rhs);

  /// \brief Allow implicit conversion to double.
  operator ScalarType() const;

  /// \brief Return value of stored volume.
  ScalarType& volume(void);
  /// \brief Return const reference to stored volume.
  const ScalarType& volume(void) const;

  /// \brief Overload assignment to assign constant value to moments.
  VolumeBase& operator=(const ScalarType a_value);

  // Unary - operator
  VolumeBase operator-(void) const;

  /// \brief Default destructor.
  ~VolumeBase(void) = default;

 private:
  ScalarType volume_m;  ///< \brief Volume of something
};

template <class ScalarType>
std::ostream& operator<<(std::ostream& out,
                         const VolumeBase<ScalarType>& a_volume);

/// \brief Overload + operator to add two geometric moments together
template <class ScalarType>
inline VolumeBase<ScalarType> operator+(const VolumeBase<ScalarType>& a_vm1,
                                        const VolumeBase<ScalarType>& a_vm2);
/// \brief Overload - operator to subtract one
/// geometric moment object from another.
template <class ScalarType>
inline VolumeBase<ScalarType> operator-(const VolumeBase<ScalarType>& a_vm1,
                                        const VolumeBase<ScalarType>& a_vm2);
/// \brief Overload * operator to multiply volume/centroid
template <class ScalarType>
inline VolumeBase<ScalarType> operator*(const ScalarType a_multiplier,
                                        const VolumeBase<ScalarType>& a_vm);
/// \brief Overload * operator to multiply volume/centroid
template <class ScalarType>
inline VolumeBase<ScalarType> operator*(const VolumeBase<ScalarType>& a_vm,
                                        const ScalarType a_multiplier);

using Volume = VolumeBase<double>;

}  // namespace IRL

#include "irl/moments/volume.tpp"

#endif  // IRL_MOMENTS_VOLUME_H_ */
