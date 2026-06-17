// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/interface_reconstruction_methods/iterative_jibben.h"
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace IRL {

iJibben_3D::iJibben_3D(const JibbenNeighborhood* a_neighborhood_pointer) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_m = a_neighborhood_pointer;
}

void iJibben_3D::getParaboloidCoefficients(void) {
  // Allocate least-squares system (using Eigen)
  const UnsignedIndex_t npolygons = neighborhood_m->size();
  Eigen::MatrixXd Amat = Eigen::MatrixXd::Zero(npolygons, 6);
  Eigen::VectorXd bvec = Eigen::VectorXd::Zero(npolygons);
  Eigen::Vector3d plane = Eigen::Vector3d::Zero();
  Eigen::VectorXd integrals = Eigen::VectorXd::Zero(6);

  for (UnsignedIndex_t n = 0; n < npolygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
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
    // Compute momonial integrals and feed to LS system
    integrals = Eigen::VectorXd::Zero(6);
    double b_dot_sum = 0.0;
    for (UnsignedIndex_t v = 0; v < nvertices; ++v) {
      const UnsignedIndex_t vn = (v + 1) % nvertices;
      const double xv = polygon[v][0], yv = polygon[v][1];
      const double xvn = polygon[vn][0], yvn = polygon[vn][1];
      integrals(0) += (xv * yvn - xvn * yv) / 2.0;
      integrals(1) += (xv + xvn) * (xv * yvn - xvn * yv) / 6.0;
      integrals(2) += (yv + yvn) * (xv * yvn - xvn * yv) / 6.0;
      integrals(3) += (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0;
      integrals(4) +=
          (yvn - yv) *
          (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
           2.0 * xv * xvn * yvn + xvn * xvn * yv + 3.0 * xvn * xvn * yvn) /
          24.0;
      integrals(5) += (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;
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

  // storing paraboloid coefficients
  for (UnsignedIndex_t i = 0; i < 6; i++) {
    coefficients_m[i] = svec(i);
  }
}

double iJibben_3D::computeMeanCurvature(const double& xi,
                                        const double& eta) const {
  // zeta = c0 + c1*xi + c2*eta + c3*xi^2 + c4*xi*eta + c5*eta^2
  const double zeta_xi = coefficients_m[1] + 2.0 * coefficients_m[3] * xi +
                         coefficients_m[4] * eta;

  const double zeta_eta = coefficients_m[2] + coefficients_m[4] * xi +
                          2.0 * coefficients_m[5] * eta;

  const double zeta_xixi = 2.0 * coefficients_m[3];
  const double zeta_etaeta = 2.0 * coefficients_m[5];
  const double zeta_xieta = coefficients_m[4];

  const double numerator =
      zeta_xixi + zeta_etaeta + zeta_xixi * zeta_eta * zeta_eta +
      zeta_etaeta * zeta_xi * zeta_xi - 2.0 * zeta_xieta * zeta_xi * zeta_eta;

  const double denominator =
      std::pow(1.0 + zeta_xi * zeta_xi + zeta_eta * zeta_eta, 1.5);

  return -numerator / denominator;
}

IRL::Normal iJibben_3D::computeNormal(const double& xi,
                                      const double& eta) const {
  const double zeta_xi = coefficients_m[1] + 2.0 * coefficients_m[3] * xi +
                         coefficients_m[4] * eta;

  const double zeta_eta = coefficients_m[2] + coefficients_m[4] * xi +
                          2.0 * coefficients_m[5] * eta;

  Normal normal(-zeta_xi, -zeta_eta, 1.0);
  normal.normalize();
  return normal;
}

// Returns a paraboloid in the global principal frame
Paraboloid iJibben_3D::makeParaboloid(
    const IRL::Pt& a_local_datum,
    const IRL::ReferenceFrame& a_local_frame) const {
  const double c0 = coefficients_m[0];
  const double c1 = coefficients_m[1];
  const double c2 = coefficients_m[2];
  const double c3 = coefficients_m[3];
  const double c4 = coefficients_m[4];
  const double c5 = coefficients_m[5];

  Eigen::Vector2d d(c1, c2);

  Eigen::Matrix2d Q;
  Q << c3, 0.5 * c4, 0.5 * c4, c5;

  const double detQ = Q.determinant();
  if (std::abs(detQ) < 1.0e-14) {
    throw std::runtime_error(
        "Cannot form principal paraboloid: quadratic matrix is singular.");
  }

  // Vertex in the original local coordinates.
  //
  // grad(zeta) = d + 2 Q s = 0
  // s_v = -1/2 Q^{-1} d
  const Eigen::Vector2d s_v = -0.5 * Q.inverse() * d;

  const double xi_v = s_v(0);
  const double eta_v = s_v(1);

  const double zeta_v = c0 + d.dot(s_v) + s_v.transpose() * Q * s_v;

  // Diagonalize the quadratic form
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(Q);

  if (eig.info() != Eigen::Success) {
    throw std::runtime_error("Eigen decomposition failed.");
  }

  const double lambda1 = eig.eigenvalues()(0);
  const double lambda2 = eig.eigenvalues()(1);

  const Eigen::Vector2d q1 = eig.eigenvectors().col(0);
  Eigen::Vector2d q2 = eig.eigenvectors().col(1);

  const double a = -lambda1;
  const double b = -lambda2;

  // Original local basis expressed globally
  const IRL::Normal& e_xi = a_local_frame[0];
  const IRL::Normal& e_eta = a_local_frame[1];
  const IRL::Normal& e_zeta = a_local_frame[2];

  // Vertex in global coordinates:
  IRL::Pt vertex_global(
      a_local_datum[0] + xi_v * e_xi[0] + eta_v * e_eta[0] + zeta_v * e_zeta[0],
      a_local_datum[1] + xi_v * e_xi[1] + eta_v * e_eta[1] + zeta_v * e_zeta[1],
      a_local_datum[2] + xi_v * e_xi[2] + eta_v * e_eta[2] +
          zeta_v * e_zeta[2]);

  // Principal tangent directions in global coordinates.
  IRL::Normal p1(q1(0) * e_xi[0] + q1(1) * e_eta[0],
                 q1(0) * e_xi[1] + q1(1) * e_eta[1],
                 q1(0) * e_xi[2] + q1(1) * e_eta[2]);

  IRL::Normal p2(q2(0) * e_xi[0] + q2(1) * e_eta[0],
                 q2(0) * e_xi[1] + q2(1) * e_eta[1],
                 q2(0) * e_xi[2] + q2(1) * e_eta[2]);

  IRL::Normal n(e_zeta[0], e_zeta[1], e_zeta[2]);

  p1.normalize();
  p2.normalize();
  n.normalize();

  // Ensure right-handed orientation
  const double cx = p1[1] * p2[2] - p1[2] * p2[1];
  const double cy = p1[2] * p2[0] - p1[0] * p2[2];
  const double cz = p1[0] * p2[1] - p1[1] * p2[0];

  const double handedness = cx * n[0] + cy * n[1] + cz * n[2];

  if (handedness < 0.0) {
    p2 = IRL::Normal(-p2[0], -p2[1], -p2[2]);
  }

  IRL::ReferenceFrame principal_frame;
  principal_frame[0] = p1;
  principal_frame[1] = p2;
  principal_frame[2] = n;

  return Paraboloid(vertex_global, principal_frame, a, b);
}

Paraboloid iJibben_3D::makeParaboloid2(
    const JibbenNeighborhood* a_neighborhood_pointer) const {
  const double a = coefficients_m[0], b = coefficients_m[1],
               c = coefficients_m[2], d = coefficients_m[3],
               e = coefficients_m[4], f = coefficients_m[5];
  const double theta = 0.5 * std::atan2(e, (safelyTiny(d - f)));
  const double cos_t = std::cos(theta), sin_t = std::sin(theta);
  const double A = -(d * cos_t * cos_t + f * sin_t * sin_t + e * cos_t * sin_t);
  const double B = -(f * cos_t * cos_t + d * sin_t * sin_t - e * cos_t * sin_t);
  const double denom_inv = 1.0 / safelyTiny(4.0 * d * f - e * e);
  const double u = (2.0 * b * f - c * e) * denom_inv;
  const double v = -(b * e - 2.0 * d * c) * denom_inv;
  const double w = -(a + (-b * b * f + b * c * e - c * c * d) * denom_inv);

  const Pt datum = a_neighborhood_pointer->getDatum();
  const ReferenceFrame frame = a_neighborhood_pointer->getReferenceFrame();
  const Pt paraboloid_datum =
      datum - u * frame[0] - v * frame[1] - w * frame[2];
  const UnitQuaternion rotation(theta, frame[2]);
  const ReferenceFrame paraboloid_frame = rotation * frame;
  Paraboloid paraboloid(paraboloid_datum, paraboloid_frame, A, B);
  return paraboloid;
}

std::pair<double, Normal> iJibben_3D::computeAveragedCurvatureAndNormal()
    const {
  const auto& polygon = neighborhood_m->getCenterPolygon();

  const UnsignedIndex_t nvertices = polygon.getNumberOfVertices();

  double curvature = 0.0;
  Normal normal(0.0, 0.0, 1.0);

  if (nvertices == 0) {
    return std::make_pair(curvature, normal);
  }

  double curvature_sum = 0.0;
  IRL::Normal normal_sum(0.0, 0.0, 0.0);

  for (UnsignedIndex_t v = 0; v < nvertices; ++v) {
    const double xi = polygon[v][0];
    const double eta = polygon[v][1];

    curvature_sum += computeMeanCurvature(xi, eta);

    const double zeta_xi = coefficients_m[1] + 2.0 * coefficients_m[3] * xi +
                           coefficients_m[4] * eta;

    const double zeta_eta = coefficients_m[2] + coefficients_m[4] * xi +
                            2.0 * coefficients_m[5] * eta;

    normal_sum[0] += -zeta_xi;
    normal_sum[1] += -zeta_eta;
    normal_sum[2] += 1.0;
  }

  normal_sum.normalize();

  curvature = curvature_sum / static_cast<double>(nvertices);
  normal = normal_sum;

  return std::make_pair(curvature, normal);
}

}  // namespace IRL
