// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_RECONSTRUCTION_TYPES_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_RECONSTRUCTION_TYPES_H_

#include <string>

#include "examples/2d_advector_debug/data.h"
#include "examples/2d_advector_debug/irl2d.h"

void getReconstruction(const std::string& a_reconstruction_method,
                       const Data<IRL2D::Moments>& a_liquid_moments,
                       const Data<IRL2D::Moments>& a_gas_moments,
                       const double a_dt, const Data<double>& a_U,
                       const Data<double>& a_V,
                       Data<IRL2D::Parabola>* a_interface);
struct ELVIRA {
  static void getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                                const Data<IRL2D::Moments>& a_gas_moments,
                                const double a_dt, const Data<double>& a_U,
                                const Data<double>& a_V,
                                Data<IRL2D::Parabola>* a_interface);
};

struct LVIRA {
  static void getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                                const Data<IRL2D::Moments>& a_gas_moments,
                                const double a_dt, const Data<double>& a_U,
                                const Data<double>& a_V,
                                Data<IRL2D::Parabola>* a_interface);
};

struct LVIRAQ {
  static void getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                                const Data<IRL2D::Moments>& a_gas_moments,
                                const double a_dt, const Data<double>& a_U,
                                const Data<double>& a_V,
                                Data<IRL2D::Parabola>* a_interface);
};

struct MOF {
  static void getReconstruction(const Data<IRL2D::Moments>& a_liquid_moments,
                                const Data<IRL2D::Moments>& a_gas_moments,
                                const double a_dt, const Data<double>& a_U,
                                const Data<double>& a_V,
                                Data<IRL2D::Parabola>* a_interface);
};

void RotateReferenceFrame(IRL2D::Parabola& interface, double angle);

void RecenterMoments(IRL2D::Moments* moments, const IRL2D::Vec& center);

void correctInterfaceBorders(Data<IRL2D::Parabola>* a_interface);

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_RECONSTRUCTION_TYPES_H_
