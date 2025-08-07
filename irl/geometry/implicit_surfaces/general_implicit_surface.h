// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_IMPLICIT_SURFACES_GENERAL_IMPLICIT_SURFACE_H_
#define IRL_GEOMETRY_IMPLICIT_SURFACES_GENERAL_IMPLICIT_SURFACE_H_

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "irl/generic_cutting/generic_cutting.h"

namespace IRL {

template <class ScalarType, size_t MAX_REFINE_LEVEL = 5>
class GeneralImplicitSurface {
 public:
  using Vec3 = Eigen::Matrix<ScalarType, 3, 1>;
  using Mat3 = Eigen::Matrix<ScalarType, 3, 3>;

  virtual ~GeneralImplicitSurface(void) = default;

  static constexpr size_t getMaxRefineLevel() { return MAX_REFINE_LEVEL; }

  // Implicit function
  virtual ScalarType F(const ScalarType& x, const ScalarType& y,
                       const ScalarType& z) const = 0;

  // Gradient
  virtual Vec3 gradF(const ScalarType& x, const ScalarType& y,
                     const ScalarType& z) const = 0;

  // Hessian
  virtual Mat3 hessF(const ScalarType& x, const ScalarType& y,
                     const ScalarType& z) const = 0;

  // Volume (optional)
  virtual bool hasVolume() const { return false; }
  virtual ScalarType volume() const {
    throw std::runtime_error("volume() not implemented for this shape.");
  }

  // Surface area (optional)
  virtual bool hasSurfaceArea() const { return false; }
  virtual ScalarType surfaceArea() const {
    throw std::runtime_error("surfaceArea() not implemented for this shape.");
  }
};

}  // namespace IRL

#endif  // IRL_GEOMETRY_IMPLICIT_SURFACES_GENERAL_IMPLICIT_SURFACE_H_