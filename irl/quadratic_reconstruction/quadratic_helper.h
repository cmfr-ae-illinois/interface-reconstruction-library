// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com> ? (I don't know if this is correct)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_QUADRATIC_RECONSTRUCTION_QUADRATIC_HELPER_H_
#define IRL_QUADRATIC_RECONSTRUCTION_QUADRATIC_HELPER_H_

#include <math.h>
#include "irl/data_structures/stack_vector.h"

namespace IRL {

template <class ScalarType>
inline StackVector<ScalarType, 2> solveQuadratic(const ScalarType a,
                                                 const ScalarType b,
                                                 const ScalarType c);

template <class GradientType, class ScalarType>
inline StackVector<std::pair<ScalarType, GradientType>, 2>
solveQuadraticWithGradient(const ScalarType a, const ScalarType b,
                           const ScalarType c, const GradientType& a_grad,
                           const GradientType& b_grad,
                           const GradientType& c_grad);

}  // namespace IRL

#include "irl/quadratic_reconstruction/quadratic_helper.tpp"

#endif // IRL_QUADRATIC_RECONSTRUCTION_QUADRATIC_HELPER_H_