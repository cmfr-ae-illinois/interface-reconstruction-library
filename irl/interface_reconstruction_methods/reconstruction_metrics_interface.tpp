// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_METRICS_INTERFACE_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_METRICS_INTERFACE_TPP_

namespace IRL {

double reconstructionMetricWithJibben3D(
    const JibbenNeighborhood& a_neighborhood_geometry) {
  Jibben_3D jibben_solver(&a_neighborhood_geometry);
  return jibben_solver.getNormalScatterMetric();
}

double angularVarianceMetricWithJibben3D(
    const JibbenNeighborhood& a_neighborhood_geometry) {
  Jibben_3D jibben_solver(&a_neighborhood_geometry);
  return jibben_solver.getAngularVariance();
}

double volumeErrorSquaredWithJibben3D(
    const JibbenNeighborhood& a_neighborhood_geometry, const double& a_dx) {
  Jibben_3D jibben_solver(&a_neighborhood_geometry);
  jibben_solver.solve2(&a_neighborhood_geometry, -1.0);
  return jibben_solver.getVolumeErrorSquared(a_dx);
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_METRICS_INTERFACE_TPP_
