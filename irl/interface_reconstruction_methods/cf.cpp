// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/interface_reconstruction_methods/cf.h"

namespace IRL {

Paraboloid CircleFit_3D::solve(const JibbenNeighborhood* a_neighborhood_pointer,
                               const double a_h) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_m = a_neighborhood_pointer;

  if (a_h > 0.0) {
    h_m = a_h;
    delta_m = 2.5 * h_m;
  } else {
    delta_m = neighborhood_m->getDelta();
    h_m = delta_m / 2.5;
  }

  return this->solve();
}

void CircleFit_3D::getIntersectionPts(
    const IRL::Polygon& a_polygon, const IRL::Plane& a_cutting_plane,
    IRL::StackVector<IRL::Pt, 2>* a_intersection_pts) {
  a_intersection_pts->resize(0);

  const int nverts = a_polygon.getNumberOfVertices();
  if (nverts < 2) {
    return;
  }

  double distance = a_cutting_plane.signedDistanceToPoint(a_polygon[0]);

  for (int n = 0; n < nverts; ++n) {
    const int next_id = (n + 1) % nverts;
    const double next_distance =
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

double CircleFit_3D::getVfWeight(const double& a_vfrac) {
  const double limit_vfrac = 0.05;

  if (a_vfrac < limit_vfrac) {
    return 0.5 - 0.5 * std::cos(M_PI * a_vfrac / limit_vfrac);
  } else if (a_vfrac > 1.0 - limit_vfrac) {
    return 0.5 - 0.5 * std::cos(M_PI * (1.0 - a_vfrac) / limit_vfrac);
  } else {
    return 1.0;
  }
}

double CircleFit_3D::getNormalWeight(const IRL::Normal& a_nref,
                                     const IRL::Normal& a_nloc) {
  const double n_dot =
      a_nref[0] * a_nloc[0] + a_nref[1] * a_nloc[1] + a_nref[2] * a_nloc[2];

  return std::max(0.0, n_dot);
}

double CircleFit_3D::getDistanceWeight(const IRL::Pt& a_pref,
                                       const IRL::Pt& a_ploc) {
  const double distance = IRL::magnitude(a_ploc - a_pref);
  const double distance_ndim = distance / delta_m;

  if (distance_ndim >= 1.0) {
    return 0.0;
  }

  return (1.0 + 4.0 * distance_ndim) * std::pow(1.0 - distance_ndim, 4.0);
}

double CircleFit_3D::getTotalWeight(const double& a_vfrac,
                                    const IRL::Normal& a_nref,
                                    const IRL::Normal& a_nloc,
                                    const IRL::Pt& a_pref,
                                    const IRL::Pt& a_ploc) {
  const double vf_weight = getVfWeight(a_vfrac);
  const double normal_weight = getNormalWeight(a_nref, a_nloc);
  const double distance_weight = getDistanceWeight(a_pref, a_ploc);

  // return vf_weight * normal_weight * distance_weight;
  return 1.0;
}

IRL::Normal CircleFit_3D::calculatePolygonNormal(
    const IRL::Polygon& a_polygon) {
  const int nverts = a_polygon.getNumberOfVertices();

  if (nverts < 3) {
    return IRL::Normal(0.0, 0.0, 0.0);
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

  return IRL::Normal(0.0, 0.0, 0.0);
}

IRL::ReferenceFrame CircleFit_3D::referenceFrameFromNormal(
    const IRL::Normal a_normal) {
  IRL::ReferenceFrame frame;

  IRL::Normal normal = a_normal;
  normal.normalize();

  int largest_dir = 0;
  if (std::fabs(normal[largest_dir]) < std::fabs(normal[1])) {
    largest_dir = 1;
  }
  if (std::fabs(normal[largest_dir]) < std::fabs(normal[2])) {
    largest_dir = 2;
  }

  if (largest_dir == 0) {
    frame[0] = IRL::crossProduct(normal, IRL::Normal(0.0, 1.0, 0.0));
  } else if (largest_dir == 1) {
    frame[0] = IRL::crossProduct(normal, IRL::Normal(0.0, 0.0, 1.0));
  } else {
    frame[0] = IRL::crossProduct(normal, IRL::Normal(1.0, 0.0, 0.0));
  }

  frame[0].normalize();
  frame[1] = IRL::crossProduct(normal, frame[0]);
  frame[2] = normal;

  return frame;
}

IRL::Pt CircleFit_3D::toLocal(const IRL::ReferenceFrame& a_local_frame,
                              const IRL::Pt& a_local_datum,
                              const IRL::Pt& a_global_pt) {
  // The local slice plane uses coordinates:
  // x_local = tangent direction,
  // y_local = normal direction,
  // z_local = in-plane unused direction.
  Eigen::Vector3d e1(a_local_frame[0][0], a_local_frame[0][1],
                     a_local_frame[0][2]);
  Eigen::Vector3d e2(a_local_frame[2][0], a_local_frame[2][1],
                     a_local_frame[2][2]);
  Eigen::Vector3d e3(a_local_frame[1][0], a_local_frame[1][1],
                     a_local_frame[1][2]);

  Eigen::Matrix3d R;
  R.col(0) = e1;
  R.col(1) = e2;
  R.col(2) = e3;

  const Eigen::Vector3d origin(a_local_datum[0], a_local_datum[1],
                               a_local_datum[2]);
  const Eigen::Vector3d global_pt(a_global_pt[0], a_global_pt[1],
                                  a_global_pt[2]);

  const Eigen::Vector3d local_pt = R.transpose() * (global_pt - origin);

  return IRL::Pt(local_pt[0], local_pt[1], local_pt[2]);
}

std::vector<double> CircleFit_3D::getTaubinMomentTerms(const double& x0,
                                                       const double& y0,
                                                       const double& dx,
                                                       const double& dy) {
  std::vector<double> terms(10, 0.0);

  // 0 -> Mzz
  terms[0] = (dx * dx * dx * dx) / 5.0 + (2.0 * dx * dx * dy * dy) / 5.0 +
             (dy * dy * dy * dy) / 5.0 + dx * dx * dx * x0 + dx * dy * dy * x0 +
             2.0 * dx * dx * x0 * x0 + (2.0 * dy * dy * x0 * x0) / 3.0 +
             2.0 * dx * x0 * x0 * x0 + x0 * x0 * x0 * x0 + dx * dx * dy * y0 +
             dy * dy * dy * y0 + (8.0 * dx * dy * x0 * y0) / 3.0 +
             2.0 * dy * x0 * x0 * y0 + (2.0 * dx * dx * y0 * y0) / 3.0 +
             2.0 * dy * dy * y0 * y0 + 2.0 * dx * x0 * y0 * y0 +
             2.0 * x0 * x0 * y0 * y0 + 2.0 * dy * y0 * y0 * y0 +
             y0 * y0 * y0 * y0;

  // 1 -> Mxz
  terms[1] = (dx * dx * dx) / 4.0 + (dx * dy * dy) / 4.0 + dx * dx * x0 +
             (dy * dy * x0) / 3.0 + (3.0 * dx * x0 * x0) / 2.0 + x0 * x0 * x0 +
             (2.0 * dx * dy * y0) / 3.0 + dy * x0 * y0 + (dx * y0 * y0) / 2.0 +
             x0 * y0 * y0;

  // 2 -> Myz
  terms[2] = (dx * dx * dy) / 4.0 + (dy * dy * dy) / 4.0 +
             (2.0 * dx * dy * x0) / 3.0 + (dy * x0 * x0) / 2.0 +
             (dx * dx * y0) / 3.0 + dy * dy * y0 + dx * x0 * y0 + x0 * x0 * y0 +
             (3.0 * dy * y0 * y0) / 2.0 + y0 * y0 * y0;

  // 3 -> Mz
  terms[3] =
      (dx * dx) / 3.0 + (dy * dy) / 3.0 + dx * x0 + x0 * x0 + dy * y0 + y0 * y0;

  // 4 -> Mxx
  terms[4] = (dx * dx) / 3.0 + dx * x0 + x0 * x0;

  // 5 -> Mxy
  terms[5] = (dx * dy) / 3.0 + (dy * x0) / 2.0 + (dx * y0) / 2.0 + x0 * y0;

  // 6 -> Mx
  terms[6] = dx / 2.0 + x0;

  // 7 -> Myy
  terms[7] = (dy * dy) / 3.0 + dy * y0 + y0 * y0;

  // 8 -> My
  terms[8] = dy / 2.0 + y0;

  // 9 -> Length normalization term.
  terms[9] = 1.0;

  return terms;
}

void CircleFit_3D::getTaubinMatrices(
    const std::vector<std::pair<IRL::Pt, IRL::Pt>>& a_end_points,
    const std::vector<double>& a_weights, Eigen::Matrix4d* a_M,
    Eigen::Matrix4d* a_C) {
  a_M->setZero();
  a_C->setZero();

  const std::size_t nsegments = a_end_points.size();

  for (std::size_t i = 0; i < nsegments; ++i) {
    const IRL::Pt x0 = a_end_points[i].first;
    const IRL::Pt x1 = a_end_points[i].second;

    const double dx = x1[0] - x0[0];
    const double dy = x1[1] - x0[1];

    const std::vector<double> terms =
        getTaubinMomentTerms(x0[0], x0[1], dx, dy);

    const double w = a_weights[i];

    const double Mzz = terms[0];
    const double Mxz = terms[1];
    const double Myz = terms[2];
    const double Mz = terms[3];
    const double Mxx = terms[4];
    const double Mxy = terms[5];
    const double Mx = terms[6];
    const double Myy = terms[7];
    const double My = terms[8];

    (*a_M)(0, 0) += w * Mzz;
    (*a_M)(0, 1) += w * Mxz;
    (*a_M)(0, 2) += w * Myz;
    (*a_M)(0, 3) += w * Mz;

    (*a_M)(1, 0) += w * Mxz;
    (*a_M)(1, 1) += w * Mxx;
    (*a_M)(1, 2) += w * Mxy;
    (*a_M)(1, 3) += w * Mx;

    (*a_M)(2, 0) += w * Myz;
    (*a_M)(2, 1) += w * Mxy;
    (*a_M)(2, 2) += w * Myy;
    (*a_M)(2, 3) += w * My;

    (*a_M)(3, 0) += w * Mz;
    (*a_M)(3, 1) += w * Mx;
    (*a_M)(3, 2) += w * My;
    (*a_M)(3, 3) += w * terms[9];
  }

  // Taubin constraint matrix.
  (*a_C)(0, 0) = 4.0 * (*a_M)(0, 3);
  (*a_C)(0, 1) = 2.0 * (*a_M)(1, 3);
  (*a_C)(0, 2) = 2.0 * (*a_M)(2, 3);

  (*a_C)(1, 0) = (*a_C)(0, 1);
  (*a_C)(1, 1) = (*a_M)(3, 3);

  (*a_C)(2, 0) = (*a_C)(0, 2);
  (*a_C)(2, 2) = (*a_M)(3, 3);
}

double CircleFit_3D::getSignedTaubinCurvature(
    const std::vector<std::pair<IRL::Pt, IRL::Pt>>& a_end_points,
    const std::vector<double>& a_weights) {
  const std::size_t nsegments = a_end_points.size();

  if (nsegments < 2 || nsegments != a_weights.size()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  Eigen::Matrix4d C = Eigen::Matrix4d::Zero();

  getTaubinMatrices(a_end_points, a_weights, &M, &C);

  Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
  ges.compute(M, C);

  const auto& evals = ges.eigenvalues();
  const auto& evecs = ges.eigenvectors();

  int best_index = -1;
  double best_lambda = std::numeric_limits<double>::infinity();

  for (int i = 0; i < 4; ++i) {
    const auto lam = evals(i);

    if (std::abs(lam.imag()) > 1.0e-9) {
      continue;
    }

    const double lambda_real = lam.real();

    if (lambda_real <= 0.0) {
      continue;
    }

    if (lambda_real < best_lambda) {
      best_lambda = lambda_real;
      best_index = i;
    }
  }

  if (best_index < 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const Eigen::Vector4cd v_c = evecs.col(best_index);
  const Eigen::Vector4d a = v_c.real();

  const double A = a(0);
  const double B = a(1);
  const double Cc = a(2);
  const double D = a(3);

  if (std::abs(A) < 1.0e-14) {
    return 0.0;
  }

  const double num = B * B + Cc * Cc - 4.0 * A * D;

  if (num <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double R = 0.5 * std::sqrt(num) / std::abs(A);

  const double yc_local = -Cc / (2.0 * A);
  const double sign = yc_local >= 0.0 ? -1.0 : 1.0;

  return sign * (1.0 / R);
}

std::tuple<double, double, double, bool> CircleFit_3D::getPrincipalCurvatures(
    const std::vector<double>& a_theta, const std::vector<double>& a_k_theta) {
  const std::size_t m = a_theta.size();

  if (m < 3 || m != a_k_theta.size()) {
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(), false};
  }

  Eigen::MatrixXd X(m, 3);
  Eigen::VectorXd y(m);

  for (std::size_t s = 0; s < m; ++s) {
    const double th = a_theta[s];

    X(s, 0) = 1.0;
    X(s, 1) = std::cos(2.0 * th);
    X(s, 2) = std::sin(2.0 * th);

    y(s) = a_k_theta[s];
  }

  const Eigen::Vector3d a = X.colPivHouseholderQr().solve(y);
  const double res_norm = (X * a - y).norm();

  if (!std::isfinite(res_norm)) {
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(), false};
  }

  const double a0 = a(0);
  const double a1 = a(1);
  const double a2 = a(2);

  const double Delta = std::sqrt(a1 * a1 + a2 * a2);

  const double eps = 1.0e-14;

  if (Delta < eps) {
    const double k1 = a0;
    const double k2 = a0;
    const double phi = 0.0;
    return {k1, k2, phi, true};
  }

  const double k1 = a0 + Delta;
  const double k2 = a0 - Delta;
  const double phi = 0.5 * std::atan2(a2, a1);

  return {k1, k2, phi, true};
}

Paraboloid CircleFit_3D::solve(void) {
  const UnsignedIndex_t npolygons = neighborhood_m->size();

  const UnsignedIndex_t nslices = 18;
  IRL::StackVector<IRL::Pt, 2> intersections;

  const IRL::Polygon& center_polygon = neighborhood_m->getCenterPolygon();
  const IRL::Pt datum = center_polygon.calculateCentroid();

  const IRL::Normal polygon_normal = calculatePolygonNormal(center_polygon);
  const IRL::ReferenceFrame polygon_frame =
      referenceFrameFromNormal(polygon_normal);

  if (npolygons < 2 || center_polygon.getNumberOfVertices() <= 2) {
    return IRL::Paraboloid(datum, polygon_frame, 0.0, 0.0);
  }

  std::vector<double> theta_list;
  std::vector<double> k_theta_list;

  theta_list.reserve(nslices);
  k_theta_list.reserve(nslices);

  for (UnsignedIndex_t s = 0; s < nslices; ++s) {
    const double theta =
        M_PI * static_cast<double>(s) / static_cast<double>(nslices);

    const IRL::UnitQuaternion rotation(theta, polygon_frame[2]);
    const IRL::ReferenceFrame rotated_polygon_frame = rotation * polygon_frame;

    const IRL::Plane slicing_plane(rotated_polygon_frame[1],
                                   rotated_polygon_frame[1] * datum);

    std::vector<std::pair<IRL::Pt, IRL::Pt>> end_points_local_frame;
    std::vector<double> weights;

    end_points_local_frame.reserve(npolygons);
    weights.reserve(npolygons);

    for (UnsignedIndex_t p = 0; p < npolygons; ++p) {
      const IRL::Polygon polygon = neighborhood_m->getPolygon(p);

      if (polygon.getNumberOfVertices() <= 2) {
        continue;
      }

      getIntersectionPts(polygon, slicing_plane, &intersections);

      if (intersections.size() != 2) {
        continue;
      }

      const IRL::Pt p1_local =
          toLocal(rotated_polygon_frame, datum, intersections[0]);
      const IRL::Pt p2_local =
          toLocal(rotated_polygon_frame, datum, intersections[1]);

      end_points_local_frame.push_back({p1_local, p2_local});

      const double vfrac = neighborhood_m->getWeight(p);
      const IRL::Normal neighbor_normal = calculatePolygonNormal(polygon);
      const IRL::Pt neighbor_centroid = polygon.calculateCentroid();

      const double weight =
          getTotalWeight(vfrac, rotated_polygon_frame[2], neighbor_normal,
                         datum, neighbor_centroid);

      weights.push_back(weight);
    }

    const double k_theta =
        getSignedTaubinCurvature(end_points_local_frame, weights);

    if (std::isfinite(k_theta)) {
      theta_list.push_back(theta);
      k_theta_list.push_back(k_theta);
    }
  }

  auto [k1, k2, phi, ok] = getPrincipalCurvatures(theta_list, k_theta_list);

  if (!ok) {
    return IRL::Paraboloid(datum, polygon_frame, 0.0, 0.0);
  }

  const IRL::UnitQuaternion rotate_phi(phi, polygon_frame[2]);
  const IRL::ReferenceFrame darboux_frame = rotate_phi * polygon_frame;

  IRL::Paraboloid paraboloid(datum, darboux_frame, 0.5 * k1, 0.5 * k2);

  return paraboloid;
}

}  // namespace IRL