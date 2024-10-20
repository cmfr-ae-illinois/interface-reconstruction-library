// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_VOF_ADVECTION_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_VOF_ADVECTION_H_

#include <string>

#include "examples/2d_advector/data.h"
#include "examples/2d_advector/irl2d.h"

void advectVOF(const std::string& a_simulation_type,
               const std::string& a_advection_method,
               const std::string& a_reconstruction_method, const double a_dt,
               const double a_time, const Data<double>& a_U,
               const Data<double>& a_V, Data<IRL2D::Moments>* a_liquid_moments,
               Data<IRL2D::Moments>* a_gas_moments,
               Data<IRL2D::Parabola>* a_interface);

void correctMomentLocations(Data<IRL2D::Moments>* a_liquid_moments);

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_VOF_ADVECTION_H_
