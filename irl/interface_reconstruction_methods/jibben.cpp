// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/interface_reconstruction_methods/jibben.h"

namespace IRL {

Paraboloid Jibben_3D::solve(const JibbenNeighborhood* a_neighborhood_pointer,
                            const double a_delta) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_m = a_neighborhood_pointer;
  delta_m = a_delta;
  return this->solve();
}

Paraboloid Jibben_3D::solve(void) {
  // Allocate least-squares system (using Eigen)
  const UnsignedIndex_t npolygons = neighborhood_m->size();
  Eigen::MatrixXd Amat = Eigen::MatrixXd::Zero(npolygons, 6);
  Eigen::VectorXd bvec = Eigen::VectorXd::Zero(npolygons);
  Eigen::Vector3d plane = Eigen::Vector3d::Zero();
  Eigen::VectorXd integrals = Eigen::VectorXd::Zero(6);

  const double delta = neighborhood_m->getDelta();
  if (delta_m < 0.0) {
    delta_m = delta;
  }

  for (UnsignedIndex_t n = 0; n < npolygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const double input_weight = neighborhood_m->getWeight(n);
    // Skip empty polygons
    const UnsignedIndex_t nvertices = polygon.getNumberOfVertices();
    if (nvertices == 0) {
      continue;
    }
    // Local polygon normal and centroid
    const Pt centroid = polygon.calculateCentroid();
    const Normal& normal = polygon.getPlaneOfExistence().normal();
    // Ignore polygons oriented more than 90 degrees from the local normal
    if (normal[2] < DBL_EPSILON) {
      continue;
    }
    // Compute distance
    const double distance = IRL::magnitude(centroid);
    if (distance > delta_m) {
      continue;
    }
    const double r = distance / delta_m;
    const double distance_weight =
        (1.0 + 4.0 * r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
    const double weight = input_weight * distance_weight;
    // Compute momonial integrals and feed to LS system
    integrals = Eigen::VectorXd::Zero(6);
    double b_dot_sum = 0.0;
    for (UnsignedIndex_t v = 0; v < nvertices; ++v) {
      const UnsignedIndex_t vn = (v + 1) % nvertices;
      const double xv = polygon[v][0], yv = polygon[v][1];
      const double xvn = polygon[vn][0], yvn = polygon[vn][1];
      integrals(0) += weight * (xv * yvn - xvn * yv) / 2.0;
      integrals(1) += weight * (xv + xvn) * (xv * yvn - xvn * yv) / 6.0;
      integrals(2) += weight * (yv + yvn) * (xv * yvn - xvn * yv) / 6.0;
      integrals(3) +=
          weight * (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0;
      integrals(4) +=
          weight * (yvn - yv) *
          (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
           2.0 * xv * xvn * yvn + xvn * xvn * yv + 3.0 * xvn * xvn * yvn) /
          24.0;
      integrals(5) +=
          weight * (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;
    }
    Amat.row(n) = integrals;
    // Compute polygon plane coefficients
    plane(0) = (centroid * normal) / normal[2];
    plane(1) = -normal[0] / normal[2];
    plane(2) = -normal[1] / normal[2];
    bvec(n) = integrals.head(3).dot(plane);
  }
  // Unconstrained LS solution
  const Eigen::VectorXd svec =
      Amat.completeOrthogonalDecomposition().pseudoInverse() * bvec;

  error_m = static_cast<double>((Amat * svec - bvec).norm()) /
            (delta_m * delta_m * delta_m);

  const double a = svec[0], b = svec[1], c = svec[2], d = svec[3], e = svec[4],
               f = svec[5];
  const double theta = 0.5 * std::atan2(e, (safelyTiny(d - f)));
  const double cos_t = std::cos(theta), sin_t = std::sin(theta);
  const double A = -(d * cos_t * cos_t + f * sin_t * sin_t + e * cos_t * sin_t);
  const double B = -(f * cos_t * cos_t + d * sin_t * sin_t - e * cos_t * sin_t);
  const double denom_inv = 1.0 / safelyTiny(4.0 * d * f - e * e);
  const double u = (2.0 * b * f - c * e) * denom_inv;
  const double v = -(b * e - 2.0 * d * c) * denom_inv;
  const double w = -(a + (-b * b * f + b * c * e - c * c * d) * denom_inv);

  const Pt datum = neighborhood_m->getDatum();
  const ReferenceFrame frame = neighborhood_m->getReferenceFrame();
  const Pt paraboloid_datum =
      datum - u * frame[0] - v * frame[1] - w * frame[2];
  const UnitQuaternion rotation(theta, frame[2]);
  const ReferenceFrame paraboloid_frame = rotation * frame;
  Paraboloid paraboloid(paraboloid_datum, paraboloid_frame, A, B);
  // Recenter paraboloid close to plic centroid
  paraboloid.regenerateAtLocation(datum);

  return paraboloid;
}

const double Jibben_3D::getFittingError(void) const { return error_m; }

}  // namespace IRL
