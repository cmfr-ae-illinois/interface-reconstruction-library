// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_METRICS_INTERFACE_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_METRICS_INTERFACE_H_

#include "irl/interface_reconstruction_methods/jibben.h"
#include "irl/interface_reconstruction_methods/jibben_neighborhood.h"
#include "irl/interface_reconstruction_methods/pu.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"

namespace IRL {

/// \brief Reconstruction metric for Jibben 3D (resolved or not)
inline double reconstructionMetricWithJibben3D(
    const JibbenNeighborhood& a_neighborhood_geometry);

/// \brief Angular variance reconstruction metric for Jibben 3D
inline double angularVarianceMetricWithJibben3D(
    const JibbenNeighborhood& a_neighborhood_geometry);

inline double volumeErrorSquaredWithJibben3D(
    const JibbenNeighborhood& a_neighborhood_geometry, const double& a_dx);

}  // namespace IRL

#include "irl/interface_reconstruction_methods/reconstruction_metrics_interface.tpp"

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_METRICS_INTERFACE_H_
