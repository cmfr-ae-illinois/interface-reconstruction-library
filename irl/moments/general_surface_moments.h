// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_GENERAL_SURFACE_MOMENTS_H_
#define IRL_MOMENTS_GENERAL_SURFACE_MOMENTS_H_

#include <array>
#include <ostream>

#include "irl/moments/volume.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief A general surface moments class
/// that stores surface moments in row-major order.
/// For 3D, this is:
/// 1, x, y, z, x^2, xy, xz, y^2, yz, z^2, x^3, x^2 y, ...
/// For 2D, this is
/// 1, x, y, x^2, xy, y^2, x^3, x^2 y, ...
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
class GeneralSurfaceMomentsBase {
 public:
  using value_type = ScalarType;
  static constexpr std::size_t linear_length =
      DIM == 3 ? (ORDER + 1) * (ORDER + 2) * (ORDER + 3) / 6
               : (ORDER + 1) * (ORDER + 2) / 2;
  using storage = std::array<ScalarType, linear_length>;

  /// \brief Default constructor.
  GeneralSurfaceMomentsBase(void);

  /// \brief Constructor from DP general moment.
  GeneralSurfaceMomentsBase(
      const GeneralSurfaceMomentsBase<ORDER, DIM, double>& a_moments);

  /// @brief  \breif Fill all moments with `a_value`
  static GeneralSurfaceMomentsBase fromScalarConstant(const ScalarType a_value);

  /// \brief Obtain GeneralSurfaceMoments from the supplied geometry.
  template <class GeometryType>
  static GeneralSurfaceMomentsBase calculateMoments(GeometryType* a_geometry);

  /// \brief Return array of central moments.
  storage calculateCentralMoments(void) const;

  /// \brief Return scale and translation invariant moments,
  /// properly normalized by M^0
  storage calculateInvariantMoments(void) const;

  /// \brief Normalize currently stored moments by normalization used for
  /// scale invariance.
  void normalizeAsInvariant(void);

  /// \brief Return reference to stored area.
  ScalarType& area(void);

  /// \brief Return const reference to stored area.
  const ScalarType area(void) const;

  /// \brief Multiply all moments by the area (zeroeth moments)
  void multiplyByArea(void);

  /// \brief Divide all moments by the area (zeroeth moments)
  void normalizeByArea(void);

  /// \brief change the datum of stored moments
  void moveMoments(const PtBase<ScalarType>& datum);

  /// \brief change the frame of reference of stored moments
  void moveAndRotateMoments(const PtBase<ScalarType>& datum,
                            const ReferenceFrameBase<ScalarType>& frame);

  /// \brief Overload += operator to update moments.
  GeneralSurfaceMomentsBase& operator+=(const GeneralSurfaceMomentsBase& a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  GeneralSurfaceMomentsBase& operator*=(const ScalarType a_rhs);

  /// \brief Overload /= operator to divide by constant double
  GeneralSurfaceMomentsBase& operator/=(const ScalarType a_rhs);

  /// \brief Overload assignment to assign constant value to moments.
  GeneralSurfaceMomentsBase& operator=(const ScalarType a_value);

  /// \brief Overload assignment to assign DP moments to moments.
  GeneralSurfaceMomentsBase& operator=(
      const GeneralSurfaceMomentsBase<ORDER, DIM, double>& a_rhs);

  /// \brief Total number of moments.
  UnsignedIndex_t size(void) const;

  /// \brief Provide mutable access to underlying moments
  ScalarType& operator[](const UnsignedIndex_t index);

  /// \brief Provide const access to underlying moments
  ScalarType operator[](const UnsignedIndex_t index) const;

  /// \brief Default destructor.
  ~GeneralSurfaceMomentsBase(void) = default;

 private:
  storage moments_m;  // Stored moments in row-major order
};

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM>
using GeneralSurfaceMoments = GeneralSurfaceMomentsBase<ORDER, DIM, double>;

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
std::ostream& operator<<(
    std::ostream& out,
    const GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType>& a_volume);

/// \brief Overload + operator to add two geometric moments together
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
inline GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType> operator+(
    const GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType>& a_vm1,
    const GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType>& a_vm2);
/// \brief Overload - operator to subtract one
/// geometric moment object from another.
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
inline GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType> operator-(
    const GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType>& a_vm1,
    const GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType>& a_vm2);
/// \brief Overload * operator to multiply moments
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
inline GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType> operator*(
    const ScalarType a_multiplier,
    const GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType>& a_vm);
/// \brief Overload * operator to multiply moments
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM, class ScalarType>
inline GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType> operator*(
    const GeneralSurfaceMomentsBase<ORDER, DIM, ScalarType>& a_vm,
    const ScalarType a_multiplier);

template <UnsignedIndex_t ORDER>
using GeneralSurfaceMoments3D = GeneralSurfaceMoments<ORDER, 3>;

template <UnsignedIndex_t ORDER>
using GeneralSurfaceMoments2D = GeneralSurfaceMoments<ORDER, 2>;

}  // namespace IRL

#include "irl/moments/general_surface_moments.tpp"

#endif  // IRL_MOMENTS_GENERAL_SURFACE_MOMENTS_H_ */
