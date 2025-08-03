// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_QUADRATIC_RECONSTRUCTION_COONS_PATCH_H_
#define IRL_QUADRATIC_RECONSTRUCTION_COONS_PATCH_H_

#include <cstdint>
#include <utility>

#include <Eigen/Dense>

#include "irl/geometry/general/pt.h"
#include "irl/quadratic_reconstruction/rational_bezier_arc.h"

namespace IRL {

/// \brief Coons patch bounded by quadratic rational bezier arcs
template <class ScalarType>
class CoonsPatchBase {
 public:
  /// \brief Default constructor.
  CoonsPatchBase(void);

  /// \brief Constructor that initializes the coons patch
  CoonsPatchBase(const RationalBezierArcBase<ScalarType>& arc_0,
                 const RationalBezierArcBase<ScalarType>& arc_1,
                 const RationalBezierArcBase<ScalarType>& arc_2,
                 const RationalBezierArcBase<ScalarType>& arc_3);

  /// \brief Return array of arcs.
  std::array<RationalBezierArcBase<ScalarType>, 4> getBoundaryArcs() const;
  /// \brief Return specific arc
  RationalBezierArcBase<ScalarType> getArc(const int edgeIndex) const;
  /// \brief Evaluates point on patch X(u,v)
  PtBase<ScalarType> evaluate(const ScalarType u, const ScalarType v) const;
  /// \brief Jacobian matrix
  Eigen::Matrix<ScalarType, 2, 2> jacobian(const ScalarType u,
                                           const ScalarType v) const;
  /// \bried determinant of jacobian matrix
  ScalarType detJacobian(const ScalarType u, const ScalarType v) const;

  /// \brief Default destructor.
  ~CoonsPatchBase(void) = default;

 private:
  RationalBezierArcBase<ScalarType> c0, c1, c2, c3;
};

using CoonsPatch = CoonsPatchBase<double>;

}  // namespace IRL

#include "irl/quadratic_reconstruction/coons_patch.tpp"

#endif  // IRL_QUADRATIC_RECONSTRUCTION_COONS_PATCH_H_
