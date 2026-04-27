// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_RECONSTRUCTION_METRICS_INTERFACE_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_RECONSTRUCTION_METRICS_INTERFACE_H_

#include "irl/c_interface/interface_reconstruction_methods/c_jibben_neighborhood.h"
#include "irl/c_interface/interface_reconstruction_methods/c_pu_neighborhood.h"
#include "irl/interface_reconstruction_methods/reconstruction_metrics_interface.h"

extern "C" {
/// \file c_localizers.h
///
/// These C-style funcions are
/// mapped to functions available in src/reconstruction_interface.h.
///
/// This file includes functions to place PlanarSeparator objects
/// in geometries. These methods differ in what they require. For
/// the individual needs of each reconstruction method,
/// it is best to constult its specific documentation.

double c_reconstructionMetricWithJibben3D(
    const c_JibbenNeigh* a_jibben_neighborhood);

double c_angularVarianceMetricWithJibben3D(
    const c_JibbenNeigh* a_jibben_neighborhood);

double c_volumeErrorSquaredWithJibben3D(
    const c_JibbenNeigh* a_jibben_neighborhood, const double* a_dx);

}  // end extern C

#endif  // IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_RECONSTRUCTION_METRICS_INTERFACE_H_
