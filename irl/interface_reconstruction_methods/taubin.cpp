// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/interface_reconstruction_methods/taubin.h"

namespace IRL {

Paraboloid Taubin_3D::solve(const JibbenNeighborhood* a_neighborhood_pointer,
                            const double a_delta) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_m = a_neighborhood_pointer;
  delta_m = a_delta;
  return this->solve();
}

void Taubin_3D::getIntersectionPts(
    const IRL::Polygon& a_polygon, const IRL::Plane& a_cutting_plane,
    IRL::StackVector<IRL::Pt, 2>* a_intersection_pts) {
  a_intersection_pts->resize(0);
  double distance = a_cutting_plane.signedDistanceToPoint(a_polygon[0]);
  for (int n = 0; n < a_polygon.getNumberOfVertices(); ++n) {
    const int next_id = (n + 1) % a_polygon.getNumberOfVertices();
    double next_distance =
        a_cutting_plane.signedDistanceToPoint(a_polygon[next_id]);
    if (distance * next_distance < 0.0) {
      a_intersection_pts->push_back(IRL::Pt::fromEdgeIntersection(
          a_polygon[n], distance, a_polygon[next_id], next_distance));
      if (a_intersection_pts->size() == 2) {
        break;
      }
    }
    distance = next_distance;
  }
}

IRL::Normal Taubin_3D::calculatePolygonNormal(const IRL::Polygon& a_polygon) {
  const int nverts = a_polygon.getNumberOfVertices();
  if (nverts < 3) {
    return IRL::Normal(0, 0, 0);
  }
  for (int n = 0; n < nverts; ++n) {
    const IRL::Pt p0 = a_polygon[n];
    const IRL::Pt p1 = a_polygon[(n + 1) % nverts];
    const IRL::Pt p2 = a_polygon[(n + 2) % nverts];
    const IRL::Normal normal = IRL::crossProduct(p1 - p0, p2 - p0);
    const double normal_magnitude = IRL::magnitude(normal);
    if (normal_magnitude > 100.0 * DBL_EPSILON) {
      return normal / normal_magnitude;
    }
  }
  return IRL::Normal(0, 0, 0);
}

IRL::ReferenceFrame Taubin_3D::referenceFrameFromNormal(
    const IRL::Normal a_normal) {
  IRL::ReferenceFrame frame;
  int largest_dir = 0;
  if (std::fabs(a_normal[largest_dir]) < std::fabs(a_normal[1]))
    largest_dir = 1;
  if (std::fabs(a_normal[largest_dir]) < std::fabs(a_normal[2]))
    largest_dir = 2;
  if (largest_dir == 0)
    frame[0] = IRL::crossProduct(a_normal, IRL::Normal(0.0, 1.0, 0.0));
  else if (largest_dir == 1)
    frame[0] = IRL::crossProduct(a_normal, IRL::Normal(0.0, 0.0, 1.0));
  else
    frame[0] = IRL::crossProduct(a_normal, IRL::Normal(1.0, 0.0, 0.0));
  frame[0].normalize();
  frame[1] = IRL::crossProduct(a_normal, frame[0]);
  frame[2] = a_normal;
  return frame;
}

double Taubin_3D::getNormalWeight(const IRL::Normal& a_nref,
                                  const IRL::Normal& a_nloc) {
  const double n_dot =
      a_nref[0] * a_nloc[0] + a_nref[1] * a_nloc[1] + a_nref[2] * a_nloc[2];
  return std::max(0.0, n_dot);
}

double Taubin_3D::getDistanceWeight(const IRL::Pt& a_pref,
                                    const IRL::Pt& a_ploc) {
  delta_m = neighborhood_m->getDelta();
  const double distance = IRL::magnitude(a_ploc - a_pref);
  const double r = distance / delta_m;
  const double distance_weight =
      r >= 1.0 ? 0.0 : (1.0 + 4.0 * r) * std::pow(1.0 - r, 4.0);
  return distance_weight;
}

void Taubin_3D::sampleLocalPoints(
    const IRL::ReferenceFrame& local_frame, const IRL::Pt& local_datum,
    const std::vector<std::pair<IRL::Pt, IRL::Pt>>& end_points_list,
    const int& num_samples_per_segment,
    std::vector<Eigen::Vector2d>* points_local_frame) {
  // rotation matrix and origin
  Eigen::Vector3d e1(local_frame[0][0], local_frame[0][1], local_frame[0][2]);
  Eigen::Vector3d e2(local_frame[2][0], local_frame[2][1], local_frame[2][2]);
  Eigen::Vector3d e3(local_frame[1][0], local_frame[1][1], local_frame[1][2]);
  Eigen::Matrix3d R;
  R.col(0) = e1;
  R.col(1) = e2;
  R.col(2) = e3;
  Eigen::Vector3d o(local_datum[0], local_datum[1], local_datum[2]);

  // sampling points
  for (int i = 0; i < end_points_list.size(); i++) {
    // starting and ending points
    Eigen::Vector3d start(end_points_list[i].first[0],
                          end_points_list[i].first[1],
                          end_points_list[i].first[2]);
    Eigen::Vector3d end(end_points_list[i].second[0],
                        end_points_list[i].second[1],
                        end_points_list[i].second[2]);

    // start and end in local rotated frame
    Eigen::Vector3d start_local = R.transpose() * (start - o);
    Eigen::Vector3d end_local = R.transpose() * (end - o);

    // sampling points
    for (int j = 0; j < num_samples_per_segment; j++) {
      double t = static_cast<double>(j) /
                 (static_cast<double>(num_samples_per_segment) - 1);
      Eigen::Vector3d pt_local = start_local * (1.0 - t) + end_local * t;
      // std::cout << pt_local.transpose() << std::endl;
      Eigen::Vector2d pt_local_2d(pt_local[0], pt_local[1]);
      points_local_frame->push_back(pt_local_2d);
    }
  }
}

double Taubin_3D::getSignedTaubinCurvature(
    const std::vector<Eigen::Vector2d>& a_points,
    const std::vector<double>& a_weights) {
  // number of points
  const int n = a_points.size();  // in local frame
  if (n < 3)
    return std::numeric_limits<double>::quiet_NaN();  // not enough points

  // moment matrix
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n; i++) {
    double xi = a_points[i].x(), yi = a_points[i].y();
    double zi = xi * xi + yi * yi;
    double w = a_weights[i];
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += w * u * u.transpose();  // weighted outer product
  }

  // constraint matrix
  Eigen::Matrix4d C;
  C.setZero();
  C(0, 0) = 4.0 * M(0, 3);
  C(0, 1) = 2.0 * M(1, 3);
  C(0, 2) = 2.0 * M(2, 3);
  C(1, 0) = C(0, 1);
  C(1, 1) = static_cast<double>(n);
  C(2, 0) = C(0, 2);
  C(2, 2) = static_cast<double>(n);

  // collinearity check for safety (again)
  const double Mxx = M(1, 1) / static_cast<double>(n);
  const double Myy = M(2, 2) / static_cast<double>(n);
  const double Mxy = M(1, 2) / static_cast<double>(n);
  const double cov_det = Mxx * Myy - Mxy * Mxy;
  if (std::abs(cov_det) < 1e-14) {
    return 0.0;  // zero curvature if collinear
  }

  // solving generalized eigenvalue problem
  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, C);

  // eigen values/vectors
  const auto& evals = ges.eigenvalues();
  const auto& evecs = ges.eigenvectors();

  // need to find eigenvector with smallest non-negative real eigenvalue
  int bestIndex = -1;
  double bestLambda = std::numeric_limits<double>::infinity();
  for (int i = 0; i < 4; i++) {
    const auto lam = evals(i);
    if (std::abs(lam.imag()) > 1e-9) continue;
    const double lambdaReal = lam.real();
    if (lambdaReal <= 0.0) continue;
    if (lambdaReal < bestLambda) {
      bestLambda = lambdaReal;
      bestIndex = i;
    }
  }
  if (bestIndex < 0) {  // no eigenvalue found
    return std::numeric_limits<double>::quiet_NaN();
  }

  // extracting eigenvector components
  Eigen::Vector4cd v_c = evecs.col(bestIndex);
  Eigen::Vector4d a = v_c.real();  // imaginary parts should be tiny
  double A = a(0);
  double B = a(1);
  double Cc = a(2);
  double D = a(3);

  if (std::abs(A) < 1e-14) {
    return 0.0;  // very small A is infinite radius again
  }

  // finding radius of curvature
  const double num = B * B + Cc * Cc - 4.0 * A * D;
  if (num <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double R = 0.5 * std::sqrt(num) / std::abs(A);

  // signed curvature
  double yc_local = -Cc / (2.0 * A);
  double sign = (yc_local >= 0.0) ? -1.0 : +1.0;
  // std::cout << 1.0 / R << std::endl;
  return sign * (1.0 / R);
}

std::tuple<double, double, double, bool> Taubin_3D::getPrincipalCurvatures(
    const std::vector<double>& a_theta, const std::vector<double>& a_k_theta) {
  const std::size_t m = a_theta.size();
  // if (m < 3 || m != a_k_theta.size()) {
  //   return {std::numeric_limits<double>::quiet_NaN(),
  //           std::numeric_limits<double>::quiet_NaN(),
  //           std::numeric_limits<double>::quiet_NaN(), false};
  // }

  // building least squares system
  Eigen::MatrixXd X(m, 3);
  Eigen::VectorXd y(m);

  for (std::size_t s = 0; s < m; ++s) {
    const double th = a_theta[s];
    const double ks = a_k_theta[s];

    X(s, 0) = 1.0;
    X(s, 1) = std::cos(2.0 * th);
    X(s, 2) = std::sin(2.0 * th);
    y(s) = ks;
  }

  // solving minimization problem
  Eigen::Vector3d a = X.colPivHouseholderQr().solve(y);
  const double res_norm = (X * a - y).norm();
  // if (std::isnan(res_norm)) {  // NaN check
  //   return {std::numeric_limits<double>::quiet_NaN(),
  //           std::numeric_limits<double>::quiet_NaN(),
  //           std::numeric_limits<double>::quiet_NaN(), false};
  // }

  const double a0 = a(0);
  const double a1 = a(1);
  const double a2 = a(2);

  const double Delta = std::sqrt(a1 * a1 + a2 * a2);
  double k1, k2, phi;

  const double eps = 1e-14;
  if (Delta < eps) {
    k1 = a0;
    k2 = a0;
    phi = 0.0;  // arbitrary
    return {k1, k2, phi, true};
  }

  k1 = a0 + Delta;
  k2 = a0 - Delta;
  phi = 0.5 * std::atan2(a2, a1);

  return {k1, k2, phi, true};
}

Paraboloid Taubin_3D::solve(void) {
  const UnsignedIndex_t npolygons = neighborhood_m->size();

  // params for slicing and circle fit
  const UnsignedIndex_t nsamples_per_segment = 10;
  const UnsignedIndex_t nslices = 18;
  IRL::StackVector<IRL::Pt, 2> intersections;

  // local frame
  const IRL::Pt datum = neighborhood_m->getCenterPolygon().calculateCentroid();
  const IRL::Normal local_normal =
      calculatePolygonNormal(neighborhood_m->getCenterPolygon());
  const IRL::ReferenceFrame local_frame =
      referenceFrameFromNormal(local_normal);

  // slicing about local normal [0,pi)
  std::vector<double> theta_list, k_theta_list;
  for (UnsignedIndex_t s = 0; s < nslices; s++) {
    // rotation angle
    double theta = M_PI * static_cast<double>(s) / static_cast<double>(nslices);

    // rotating local polygon frame
    const IRL::UnitQuaternion rotation(theta, local_frame[2]);
    const IRL::ReferenceFrame rotated_local_frame = rotation * local_frame;

    // slicing plane
    const IRL::Plane slicing_plane(rotated_local_frame[1],
                                   rotated_local_frame[1] * datum);

    // slicing polygons in the neighborhood
    std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_list;
    std::vector<double> weights;

    for (UnsignedIndex_t p = 0; p < npolygons; p++) {
      getIntersectionPts(neighborhood_m->getPolygon(p), slicing_plane,
                         &intersections);
      if (intersections.size() != 2) continue;
      // end points for sampling
      IRL::Pt start_point = intersections[0];
      IRL::Pt end_point = intersections[1];
      end_points_list.push_back({start_point, end_point});
      // weights
      double vf_weight = neighborhood_m->getWeight(p);
      IRL::Normal neighbor_normal =
          calculatePolygonNormal(neighborhood_m->getPolygon(p));
      double n_weight = getNormalWeight(local_frame[2], neighbor_normal);
      double d_weight = getDistanceWeight(
          datum, neighborhood_m->getPolygon(p).calculateCentroid());
      double weight = vf_weight * d_weight * n_weight;
      weights.insert(weights.end(), nsamples_per_segment, weight);
    }

    // sampling points along line segments in rotated frame
    const int num_points = end_points_list.size() * nsamples_per_segment;
    std::vector<Eigen::Vector2d> points_rotated_frame;
    points_rotated_frame.reserve(num_points);
    sampleLocalPoints(rotated_local_frame, datum, end_points_list,
                      nsamples_per_segment, &points_rotated_frame);

    // taubin circle fit
    double k_theta = getSignedTaubinCurvature(points_rotated_frame, weights);
    if (std::isfinite(k_theta)) {
      theta_list.push_back(theta);
      k_theta_list.push_back(k_theta);
    }
  }

  auto [k1, k2, phi, ok] = getPrincipalCurvatures(theta_list, k_theta_list);

  // Darboux frame
  const IRL::UnitQuaternion rotate_phi(phi, local_frame[2]);
  const IRL::ReferenceFrame darboux_frame = rotate_phi * local_frame;

  // new paraboloid
  IRL::Paraboloid paraboloid(datum, darboux_frame, 0.5 * k1, 0.5 * k2);

  return paraboloid;
}

}  // namespace IRL
