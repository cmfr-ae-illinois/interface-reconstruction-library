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

Jibben_3D::Jibben_3D(const JibbenNeighborhood* a_neighborhood_pointer,
                     const double a_delta) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_m = a_neighborhood_pointer;
  delta_m = a_delta;
}

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

  // storing paraboloid coefficients
  for (UnsignedIndex_t i = 0; i < 6; i++) {
    coefficients_m[i] = svec(i);
  }

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

const double Jibben_3D::getVolumeError(const double& dx) const {
  const UnsignedIndex_t n_polygons = neighborhood_m->size();

  // looping over neighboring polygons
  double I = 0.0;
  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const UnsignedIndex_t n_vertices = polygon.getNumberOfVertices();
    if (n_vertices == 0) continue;
    const Pt centroid = polygon.calculateCentroid();
    const Normal& normal = polygon.getPlaneOfExistence().normal();
    if (normal[2] < DBL_EPSILON) continue;

    // parameters of the plane
    std::array<double, 3> coefficients_plane{};
    coefficients_plane[0] = centroid * normal / normal[2];
    coefficients_plane[1] = -normal[0] / normal[2];
    coefficients_plane[2] = -normal[1] / normal[2];

    // looping over vertices of a polygon to compute the integrals
    std::array<double, 6> integrals{};
    for (UnsignedIndex_t v = 0; v < n_vertices; v++) {
      const UnsignedIndex_t vn = (v + 1) % n_vertices;
      const double xv = polygon[v][0], yv = polygon[v][1];
      const double xvn = polygon[vn][0], yvn = polygon[vn][1];
      integrals[0] += (xv * yvn - xvn * yv) / 2.0;
      integrals[1] += (xv + xvn) * (xv * yvn - xvn * yv) / 6.0;
      integrals[2] += (yv + yvn) * (xv * yvn - xvn * yv) / 6.0;
      integrals[3] += (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0;
      integrals[4] +=
          (yvn - yv) *
          (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
           2.0 * xv * xvn * yvn + xvn * xvn * yv + 3.0 * xvn * xvn * yvn) /
          24.0;
      integrals[5] += (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;
    }

    double I_parabola = 0.0, I_plane = 0.0;
    for (UnsignedIndex_t c = 0; c < coefficients_m.size(); c++) {
      if (c < 3) {
        I_plane += coefficients_plane[c] * integrals[c];
      }
      I_parabola += coefficients_m[c] * integrals[c];
    }

    // difference in volume between plane and paraboloid
    I += std::abs(I_parabola - I_plane);
  }
  I = std::abs(I);

  // projected area of polygons on local xy plane
  double projected_area = 0.0;
  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const Normal& normal = polygon.getPlaneOfExistence().normal();
    if (normal[2] < DBL_EPSILON) continue;
    const UnsignedIndex_t n_vertices = polygon.getNumberOfVertices();
    Polygon projected_polygon;
    for (UnsignedIndex_t v = 0; v < n_vertices; v++) {
      projected_polygon.addVertex(Pt(polygon[v][0], polygon[v][1], 0.0));
    }
    projected_polygon.calculateAndSetPlaneOfExistence();
    projected_area += std::abs(projected_polygon.calculateVolume());
  }

  // non-dimensionalization
  I *= 1.0 / (projected_area * dx);
  // I /= (dx * dx * dx);

  return I;
}

const double Jibben_3D::getVolumeErrorSquared(const double& dx) const {
  const UnsignedIndex_t n_polygons = neighborhood_m->size();

  // looping over neighboring polygons
  double I = 0.0;
  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const UnsignedIndex_t n_vertices = polygon.getNumberOfVertices();
    if (n_vertices == 0) continue;
    const Pt centroid = polygon.calculateCentroid();
    const Normal& normal = polygon.getPlaneOfExistence().normal();
    if (normal[2] < DBL_EPSILON) continue;

    // coefficients of the plane
    std::array<double, 3> b{};
    b[0] = centroid * normal / normal[2];
    b[1] = -normal[0] / normal[2];
    b[2] = -normal[1] / normal[2];

    // coefficients of the paraboloid
    std::array<double, 6> c{};
    for (UnsignedIndex_t i = 0; i < 6; i++) {
      c[i] = coefficients_m[i];
    }

    // looping over vertices of a polygon to compute the integrals
    std::array<double, 15> integrals{};
    for (UnsignedIndex_t v = 0; v < n_vertices; v++) {
      const UnsignedIndex_t vn = (v + 1) % n_vertices;
      const double xv = polygon[v][0], yv = polygon[v][1];
      const double xvn = polygon[vn][0], yvn = polygon[vn][1];
      const double dxv = xvn - xv, dyv = yvn - yv;

      // ∫dA
      integrals[0] += (xv * yvn - xvn * yv) / 2.0;

      // ∫x dA
      integrals[1] += (xv + xvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫y dA
      integrals[2] += (yv + yvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫x^2 dA
      integrals[3] += (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0;

      // ∫xy dA
      integrals[4] +=
          (yvn - yv) *
          (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
           2.0 * xv * xvn * yvn + xvn * xvn * yv + 3.0 * xvn * xvn * yvn) /
          24.0;

      // ∫y^2 dA
      integrals[5] += (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;

      // ∫x^3 dA
      integrals[6] +=
          (dyv * ((dxv * dxv * dxv * dxv) / 5. + dxv * dxv * dxv * xv +
                  2. * (dxv * dxv) * (xv * xv) + 2. * dxv * (xv * xv * xv) +
                  xv * xv * xv * xv)) /
          4.;

      // ∫x^2 y dA
      integrals[7] +=
          (dyv * (5. * yv * (dxv + 2. * xv) *
                      (dxv * dxv + 2. * dxv * xv + 2. * (xv * xv)) +
                  dyv * (4. * (dxv * dxv * dxv) + 15. * (dxv * dxv) * xv +
                         20. * dxv * (xv * xv) + 10. * (xv * xv * xv)))) /
          60.;

      // ∫x y^2 dA
      integrals[8] += -0.016666666666666666 *
                      (dxv * (10. * (yv * yv * yv) * (dxv + 2. * xv) +
                              10. * dyv * (yv * yv) * (2. * dxv + 3. * xv) +
                              5. * (dyv * dyv) * yv * (3. * dxv + 4. * xv) +
                              dyv * dyv * dyv * (4. * dxv + 5. * xv)));

      // ∫y^3 dA
      integrals[9] +=
          -0.25 * (dxv * ((dyv * dyv * dyv * dyv) / 5. + dyv * dyv * dyv * yv +
                          2. * (dyv * dyv) * (yv * yv) +
                          2. * dyv * (yv * yv * yv) + yv * yv * yv * yv));

      // ∫x^4 dA
      integrals[10] +=
          (dyv *
           ((dxv * dxv * dxv * dxv * dxv) / 6. + dxv * dxv * dxv * dxv * xv +
            (5. * (dxv * dxv * dxv) * (xv * xv)) / 2. +
            (10. * (dxv * dxv) * (xv * xv * xv)) / 3. +
            (5. * dxv * (xv * xv * xv * xv)) / 2. + xv * xv * xv * xv * xv)) /
          5.;

      // ∫x^3 y dA
      integrals[11] +=
          (dyv * (dxv * dxv * dxv * dxv * (5. * dyv + 6. * yv) +
                  6. * (dxv * dxv * dxv) * (4. * dyv + 5. * yv) * xv +
                  15. * (dxv * dxv) * (3. * dyv + 4. * yv) * (xv * xv) +
                  20. * dxv * (2. * dyv + 3. * yv) * (xv * xv * xv) +
                  15. * (dyv + 2. * yv) * (xv * xv * xv * xv))) /
          120.;

      // ∫x^2 y^2 dA
      integrals[12] +=
          -0.005555555555555556 *
          (dxv * (20. * (yv * yv * yv) *
                      (dxv * dxv + 3. * dxv * xv + 3. * (xv * xv)) +
                  15. * dyv * (yv * yv) *
                      (3. * (dxv * dxv) + 8. * dxv * xv + 6. * (xv * xv)) +
                  6. * (dyv * dyv) * yv *
                      (6. * (dxv * dxv) + 15. * dxv * xv + 10. * (xv * xv)) +
                  dyv * dyv * dyv *
                      (10. * (dxv * dxv) + 24. * dxv * xv + 15. * (xv * xv))));

      // ∫x y^3 dA
      integrals[13] +=
          -0.008333333333333333 *
          (dxv * (15. * (yv * yv * yv * yv) * (dxv + 2. * xv) +
                  20. * dyv * (yv * yv * yv) * (2. * dxv + 3. * xv) +
                  15. * (dyv * dyv) * (yv * yv) * (3. * dxv + 4. * xv) +
                  6. * (dyv * dyv * dyv) * yv * (4. * dxv + 5. * xv) +
                  dyv * dyv * dyv * dyv * (5. * dxv + 6. * xv)));

      // ∫y^4 dA
      integrals[14] +=
          -0.2 *
          (dxv *
           ((dyv * dyv * dyv * dyv * dyv) / 6. + dyv * dyv * dyv * dyv * yv +
            (5. * (dyv * dyv * dyv) * (yv * yv)) / 2. +
            (10. * (dyv * dyv) * (yv * yv * yv)) / 3. +
            (5. * dyv * (yv * yv * yv * yv)) / 2. + yv * yv * yv * yv * yv));
    }

    I += integrals[0] * (b[0] * b[0] - 2. * b[0] * c[0] + c[0] * c[0]) +
         integrals[1] * (2. * b[0] * b[1] - 2. * b[1] * c[0] -
                         2. * b[0] * c[1] + 2. * c[0] * c[1]) +
         integrals[2] * (2. * b[0] * b[2] - 2. * b[2] * c[0] -
                         2. * b[0] * c[2] + 2. * c[0] * c[2]) +
         integrals[10] * (c[3] * c[3]) +
         integrals[3] * (b[1] * b[1] - 2. * b[1] * c[1] + c[1] * c[1] -
                         2. * b[0] * c[3] + 2. * c[0] * c[3]) +
         integrals[6] * (-2. * b[1] * c[3] + 2. * c[1] * c[3]) +
         2.0 * integrals[11] * c[3] * c[4] +
         integrals[4] *
             (2. * b[1] * b[2] - 2. * b[2] * c[1] - 2. * b[1] * c[2] +
              2. * c[1] * c[2] - 2. * b[0] * c[4] + 2. * c[0] * c[4]) +
         integrals[7] * (-2. * b[2] * c[3] + 2. * c[2] * c[3] -
                         2. * b[1] * c[4] + 2. * c[1] * c[4]) +
         2.0 * integrals[13] * c[4] * c[5] + integrals[14] * (c[5] * c[5]) +
         integrals[5] * (b[2] * b[2] - 2. * b[2] * c[2] + c[2] * c[2] -
                         2. * b[0] * c[5] + 2. * c[0] * c[5]) +
         integrals[8] * (-2. * b[2] * c[4] + 2. * c[2] * c[4] -
                         2. * b[1] * c[5] + 2. * c[1] * c[5]) +
         integrals[9] * (-2. * b[2] * c[5] + 2. * c[2] * c[5]) +
         integrals[12] * (c[4] * c[4] + 2. * c[3] * c[5]);
  }

  // projected area of polygons on local xy plane
  double projected_area = 0.0;
  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const Normal& normal = polygon.getPlaneOfExistence().normal();
    if (normal[2] < DBL_EPSILON) continue;
    const UnsignedIndex_t n_vertices = polygon.getNumberOfVertices();
    Polygon projected_polygon;
    for (UnsignedIndex_t v = 0; v < n_vertices; v++) {
      projected_polygon.addVertex(Pt(polygon[v][0], polygon[v][1], 0.0));
    }
    projected_polygon.calculateAndSetPlaneOfExistence();
    projected_area += std::abs(projected_polygon.calculateVolume());
  }

  // non-dimensionalization
  I = I / (projected_area * dx * dx);

  return I;
}

const double Jibben_3D::getVolumeErrorSquared_w1(const double& dx) const {
  const UnsignedIndex_t n_polygons = neighborhood_m->size();

  // looping over neighboring polygons
  double I = 0.0;
  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const UnsignedIndex_t n_vertices = polygon.getNumberOfVertices();
    if (n_vertices == 0) continue;
    const Pt centroid = polygon.calculateCentroid();
    const Normal& normal = polygon.getPlaneOfExistence().normal();
    if (normal[2] < DBL_EPSILON) continue;

    const double distance = IRL::magnitude(centroid);
    if (distance > delta_m) {
      continue;
    }
    const double r = distance / delta_m;
    const double distance_weight =
        (1.0 + 4.0 * r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
    const double weight = distance_weight;

    // coefficients of the plane
    std::array<double, 3> b{};
    b[0] = centroid * normal / normal[2];
    b[1] = -normal[0] / normal[2];
    b[2] = -normal[1] / normal[2];

    // coefficients of the paraboloid
    std::array<double, 6> c{};
    for (UnsignedIndex_t i = 0; i < 6; i++) {
      c[i] = coefficients_m[i];
    }

    // looping over vertices of a polygon to compute the integrals
    std::array<double, 15> integrals{};
    for (UnsignedIndex_t v = 0; v < n_vertices; v++) {
      const UnsignedIndex_t vn = (v + 1) % n_vertices;
      const double xv = polygon[v][0], yv = polygon[v][1];
      const double xvn = polygon[vn][0], yvn = polygon[vn][1];
      const double dxv = xvn - xv, dyv = yvn - yv;

      // ∫dA
      integrals[0] += weight * (xv * yvn - xvn * yv) / 2.0;

      // ∫x dA
      integrals[1] += weight * (xv + xvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫y dA
      integrals[2] += weight * (yv + yvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫x^2 dA
      integrals[3] +=
          weight * (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0;

      // ∫xy dA
      integrals[4] +=
          weight * (yvn - yv) *
          (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
           2.0 * xv * xvn * yvn + xvn * xvn * yv + 3.0 * xvn * xvn * yvn) /
          24.0;

      // ∫y^2 dA
      integrals[5] +=
          weight * (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;

      // ∫x^3 dA
      integrals[6] +=
          weight *
          (dyv * ((dxv * dxv * dxv * dxv) / 5. + dxv * dxv * dxv * xv +
                  2. * (dxv * dxv) * (xv * xv) + 2. * dxv * (xv * xv * xv) +
                  xv * xv * xv * xv)) /
          4.;

      // ∫x^2 y dA
      integrals[7] +=
          weight *
          (dyv * (5. * yv * (dxv + 2. * xv) *
                      (dxv * dxv + 2. * dxv * xv + 2. * (xv * xv)) +
                  dyv * (4. * (dxv * dxv * dxv) + 15. * (dxv * dxv) * xv +
                         20. * dxv * (xv * xv) + 10. * (xv * xv * xv)))) /
          60.;

      // ∫x y^2 dA
      integrals[8] += weight * -0.016666666666666666 *
                      (dxv * (10. * (yv * yv * yv) * (dxv + 2. * xv) +
                              10. * dyv * (yv * yv) * (2. * dxv + 3. * xv) +
                              5. * (dyv * dyv) * yv * (3. * dxv + 4. * xv) +
                              dyv * dyv * dyv * (4. * dxv + 5. * xv)));

      // ∫y^3 dA
      integrals[9] +=
          weight * -0.25 *
          (dxv * ((dyv * dyv * dyv * dyv) / 5. + dyv * dyv * dyv * yv +
                  2. * (dyv * dyv) * (yv * yv) + 2. * dyv * (yv * yv * yv) +
                  yv * yv * yv * yv));

      // ∫x^4 dA
      integrals[10] +=
          weight *
          (dyv *
           ((dxv * dxv * dxv * dxv * dxv) / 6. + dxv * dxv * dxv * dxv * xv +
            (5. * (dxv * dxv * dxv) * (xv * xv)) / 2. +
            (10. * (dxv * dxv) * (xv * xv * xv)) / 3. +
            (5. * dxv * (xv * xv * xv * xv)) / 2. + xv * xv * xv * xv * xv)) /
          5.;

      // ∫x^3 y dA
      integrals[11] +=
          weight *
          (dyv * (dxv * dxv * dxv * dxv * (5. * dyv + 6. * yv) +
                  6. * (dxv * dxv * dxv) * (4. * dyv + 5. * yv) * xv +
                  15. * (dxv * dxv) * (3. * dyv + 4. * yv) * (xv * xv) +
                  20. * dxv * (2. * dyv + 3. * yv) * (xv * xv * xv) +
                  15. * (dyv + 2. * yv) * (xv * xv * xv * xv))) /
          120.;

      // ∫x^2 y^2 dA
      integrals[12] +=
          weight * -0.005555555555555556 *
          (dxv * (20. * (yv * yv * yv) *
                      (dxv * dxv + 3. * dxv * xv + 3. * (xv * xv)) +
                  15. * dyv * (yv * yv) *
                      (3. * (dxv * dxv) + 8. * dxv * xv + 6. * (xv * xv)) +
                  6. * (dyv * dyv) * yv *
                      (6. * (dxv * dxv) + 15. * dxv * xv + 10. * (xv * xv)) +
                  dyv * dyv * dyv *
                      (10. * (dxv * dxv) + 24. * dxv * xv + 15. * (xv * xv))));

      // ∫x y^3 dA
      integrals[13] +=
          weight * -0.008333333333333333 *
          (dxv * (15. * (yv * yv * yv * yv) * (dxv + 2. * xv) +
                  20. * dyv * (yv * yv * yv) * (2. * dxv + 3. * xv) +
                  15. * (dyv * dyv) * (yv * yv) * (3. * dxv + 4. * xv) +
                  6. * (dyv * dyv * dyv) * yv * (4. * dxv + 5. * xv) +
                  dyv * dyv * dyv * dyv * (5. * dxv + 6. * xv)));

      // ∫y^4 dA
      integrals[14] +=
          weight * -0.2 *
          (dxv *
           ((dyv * dyv * dyv * dyv * dyv) / 6. + dyv * dyv * dyv * dyv * yv +
            (5. * (dyv * dyv * dyv) * (yv * yv)) / 2. +
            (10. * (dyv * dyv) * (yv * yv * yv)) / 3. +
            (5. * dyv * (yv * yv * yv * yv)) / 2. + yv * yv * yv * yv * yv));
    }

    I += integrals[0] * (b[0] * b[0] - 2. * b[0] * c[0] + c[0] * c[0]) +
         integrals[1] * (2. * b[0] * b[1] - 2. * b[1] * c[0] -
                         2. * b[0] * c[1] + 2. * c[0] * c[1]) +
         integrals[2] * (2. * b[0] * b[2] - 2. * b[2] * c[0] -
                         2. * b[0] * c[2] + 2. * c[0] * c[2]) +
         integrals[10] * (c[3] * c[3]) +
         integrals[3] * (b[1] * b[1] - 2. * b[1] * c[1] + c[1] * c[1] -
                         2. * b[0] * c[3] + 2. * c[0] * c[3]) +
         integrals[6] * (-2. * b[1] * c[3] + 2. * c[1] * c[3]) +
         2.0 * integrals[11] * c[3] * c[4] +
         integrals[4] *
             (2. * b[1] * b[2] - 2. * b[2] * c[1] - 2. * b[1] * c[2] +
              2. * c[1] * c[2] - 2. * b[0] * c[4] + 2. * c[0] * c[4]) +
         integrals[7] * (-2. * b[2] * c[3] + 2. * c[2] * c[3] -
                         2. * b[1] * c[4] + 2. * c[1] * c[4]) +
         2.0 * integrals[13] * c[4] * c[5] + integrals[14] * (c[5] * c[5]) +
         integrals[5] * (b[2] * b[2] - 2. * b[2] * c[2] + c[2] * c[2] -
                         2. * b[0] * c[5] + 2. * c[0] * c[5]) +
         integrals[8] * (-2. * b[2] * c[4] + 2. * c[2] * c[4] -
                         2. * b[1] * c[5] + 2. * c[1] * c[5]) +
         integrals[9] * (-2. * b[2] * c[5] + 2. * c[2] * c[5]) +
         integrals[12] * (c[4] * c[4] + 2. * c[3] * c[5]);
  }

  // projected area of polygons on local xy plane
  double projected_area = 0.0;
  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const Normal& normal = polygon.getPlaneOfExistence().normal();
    if (normal[2] < DBL_EPSILON) continue;

    const Pt centroid = polygon.calculateCentroid();
    const double distance = IRL::magnitude(centroid);
    if (distance > delta_m) {
      continue;
    }
    const double r = distance / delta_m;
    const double distance_weight =
        (1.0 + 4.0 * r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
    const double weight = distance_weight;

    const UnsignedIndex_t n_vertices = polygon.getNumberOfVertices();
    Polygon projected_polygon;
    for (UnsignedIndex_t v = 0; v < n_vertices; v++) {
      projected_polygon.addVertex(Pt(polygon[v][0], polygon[v][1], 0.0));
    }
    projected_polygon.calculateAndSetPlaneOfExistence();
    projected_area += weight * std::abs(projected_polygon.calculateVolume());
  }

  // non-dimensionalization
  I = I / (projected_area * dx * dx);
  return I;
}

const double Jibben_3D::getVolumeErrorSquared_w2(const double& dx) const {
  const UnsignedIndex_t n_polygons = neighborhood_m->size();

  // looping over neighboring polygons
  double I = 0.0;
  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const UnsignedIndex_t n_vertices = polygon.getNumberOfVertices();
    if (n_vertices == 0) continue;
    const Pt centroid = polygon.calculateCentroid();
    const Normal& normal = polygon.getPlaneOfExistence().normal();
    if (normal[2] < DBL_EPSILON) continue;

    const double distance = IRL::magnitude(centroid);
    if (distance > delta_m) {
      continue;
    }
    const double r = distance / delta_m;
    const double distance_weight =
        (1.0 + 4.0 * r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
    const double weight = distance_weight;

    // coefficients of the plane
    std::array<double, 3> b{};
    b[0] = centroid * normal / normal[2];
    b[1] = -normal[0] / normal[2];
    b[2] = -normal[1] / normal[2];

    // coefficients of the paraboloid
    std::array<double, 6> c{};
    for (UnsignedIndex_t i = 0; i < 6; i++) {
      c[i] = coefficients_m[i];
    }

    // looping over vertices of a polygon to compute the integrals
    std::array<double, 15> integrals{};
    for (UnsignedIndex_t v = 0; v < n_vertices; v++) {
      const UnsignedIndex_t vn = (v + 1) % n_vertices;
      const double xv = polygon[v][0], yv = polygon[v][1];
      const double xvn = polygon[vn][0], yvn = polygon[vn][1];
      const double dxv = xvn - xv, dyv = yvn - yv;

      // ∫dA
      integrals[0] += weight * (xv * yvn - xvn * yv) / 2.0;

      // ∫x dA
      integrals[1] += weight * (xv + xvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫y dA
      integrals[2] += weight * (yv + yvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫x^2 dA
      integrals[3] +=
          weight * (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0;

      // ∫xy dA
      integrals[4] +=
          weight * (yvn - yv) *
          (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
           2.0 * xv * xvn * yvn + xvn * xvn * yv + 3.0 * xvn * xvn * yvn) /
          24.0;

      // ∫y^2 dA
      integrals[5] +=
          weight * (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;

      // ∫x^3 dA
      integrals[6] +=
          weight *
          (dyv * ((dxv * dxv * dxv * dxv) / 5. + dxv * dxv * dxv * xv +
                  2. * (dxv * dxv) * (xv * xv) + 2. * dxv * (xv * xv * xv) +
                  xv * xv * xv * xv)) /
          4.;

      // ∫x^2 y dA
      integrals[7] +=
          weight *
          (dyv * (5. * yv * (dxv + 2. * xv) *
                      (dxv * dxv + 2. * dxv * xv + 2. * (xv * xv)) +
                  dyv * (4. * (dxv * dxv * dxv) + 15. * (dxv * dxv) * xv +
                         20. * dxv * (xv * xv) + 10. * (xv * xv * xv)))) /
          60.;

      // ∫x y^2 dA
      integrals[8] += weight * -0.016666666666666666 *
                      (dxv * (10. * (yv * yv * yv) * (dxv + 2. * xv) +
                              10. * dyv * (yv * yv) * (2. * dxv + 3. * xv) +
                              5. * (dyv * dyv) * yv * (3. * dxv + 4. * xv) +
                              dyv * dyv * dyv * (4. * dxv + 5. * xv)));

      // ∫y^3 dA
      integrals[9] +=
          weight * -0.25 *
          (dxv * ((dyv * dyv * dyv * dyv) / 5. + dyv * dyv * dyv * yv +
                  2. * (dyv * dyv) * (yv * yv) + 2. * dyv * (yv * yv * yv) +
                  yv * yv * yv * yv));

      // ∫x^4 dA
      integrals[10] +=
          weight *
          (dyv *
           ((dxv * dxv * dxv * dxv * dxv) / 6. + dxv * dxv * dxv * dxv * xv +
            (5. * (dxv * dxv * dxv) * (xv * xv)) / 2. +
            (10. * (dxv * dxv) * (xv * xv * xv)) / 3. +
            (5. * dxv * (xv * xv * xv * xv)) / 2. + xv * xv * xv * xv * xv)) /
          5.;

      // ∫x^3 y dA
      integrals[11] +=
          weight *
          (dyv * (dxv * dxv * dxv * dxv * (5. * dyv + 6. * yv) +
                  6. * (dxv * dxv * dxv) * (4. * dyv + 5. * yv) * xv +
                  15. * (dxv * dxv) * (3. * dyv + 4. * yv) * (xv * xv) +
                  20. * dxv * (2. * dyv + 3. * yv) * (xv * xv * xv) +
                  15. * (dyv + 2. * yv) * (xv * xv * xv * xv))) /
          120.;

      // ∫x^2 y^2 dA
      integrals[12] +=
          weight * -0.005555555555555556 *
          (dxv * (20. * (yv * yv * yv) *
                      (dxv * dxv + 3. * dxv * xv + 3. * (xv * xv)) +
                  15. * dyv * (yv * yv) *
                      (3. * (dxv * dxv) + 8. * dxv * xv + 6. * (xv * xv)) +
                  6. * (dyv * dyv) * yv *
                      (6. * (dxv * dxv) + 15. * dxv * xv + 10. * (xv * xv)) +
                  dyv * dyv * dyv *
                      (10. * (dxv * dxv) + 24. * dxv * xv + 15. * (xv * xv))));

      // ∫x y^3 dA
      integrals[13] +=
          weight * -0.008333333333333333 *
          (dxv * (15. * (yv * yv * yv * yv) * (dxv + 2. * xv) +
                  20. * dyv * (yv * yv * yv) * (2. * dxv + 3. * xv) +
                  15. * (dyv * dyv) * (yv * yv) * (3. * dxv + 4. * xv) +
                  6. * (dyv * dyv * dyv) * yv * (4. * dxv + 5. * xv) +
                  dyv * dyv * dyv * dyv * (5. * dxv + 6. * xv)));

      // ∫y^4 dA
      integrals[14] +=
          weight * -0.2 *
          (dxv *
           ((dyv * dyv * dyv * dyv * dyv) / 6. + dyv * dyv * dyv * dyv * yv +
            (5. * (dyv * dyv * dyv) * (yv * yv)) / 2. +
            (10. * (dyv * dyv) * (yv * yv * yv)) / 3. +
            (5. * dyv * (yv * yv * yv * yv)) / 2. + yv * yv * yv * yv * yv));
    }

    double volume_term =
        integrals[0] * (b[0] * b[0] - 2. * b[0] * c[0] + c[0] * c[0]) +
        integrals[1] * (2. * b[0] * b[1] - 2. * b[1] * c[0] - 2. * b[0] * c[1] +
                        2. * c[0] * c[1]) +
        integrals[2] * (2. * b[0] * b[2] - 2. * b[2] * c[0] - 2. * b[0] * c[2] +
                        2. * c[0] * c[2]) +
        integrals[10] * (c[3] * c[3]) +
        integrals[3] * (b[1] * b[1] - 2. * b[1] * c[1] + c[1] * c[1] -
                        2. * b[0] * c[3] + 2. * c[0] * c[3]) +
        integrals[6] * (-2. * b[1] * c[3] + 2. * c[1] * c[3]) +
        2.0 * integrals[11] * c[3] * c[4] +
        integrals[4] *
            (2. * b[1] * b[2] - 2. * b[2] * c[1] - 2. * b[1] * c[2] +
             2. * c[1] * c[2] - 2. * b[0] * c[4] + 2. * c[0] * c[4]) +
        integrals[7] * (-2. * b[2] * c[3] + 2. * c[2] * c[3] -
                        2. * b[1] * c[4] + 2. * c[1] * c[4]) +
        2.0 * integrals[13] * c[4] * c[5] + integrals[14] * (c[5] * c[5]) +
        integrals[5] * (b[2] * b[2] - 2. * b[2] * c[2] + c[2] * c[2] -
                        2. * b[0] * c[5] + 2. * c[0] * c[5]) +
        integrals[8] * (-2. * b[2] * c[4] + 2. * c[2] * c[4] -
                        2. * b[1] * c[5] + 2. * c[1] * c[5]) +
        integrals[9] * (-2. * b[2] * c[5] + 2. * c[2] * c[5]) +
        integrals[12] * (c[4] * c[4] + 2. * c[3] * c[5]);

    // projected area
    Polygon projected_polygon;
    for (UnsignedIndex_t v = 0; v < n_vertices; v++) {
      projected_polygon.addVertex(Pt(polygon[v][0], polygon[v][1], 0.0));
    }
    projected_polygon.calculateAndSetPlaneOfExistence();
    double projected_area = std::abs(projected_polygon.calculateVolume());

    I += (volume_term / projected_area);
  }

  // projected area of polygons on local xy plane
  double weight_sum = 0.0;
  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const Normal& normal = polygon.getPlaneOfExistence().normal();
    if (normal[2] < DBL_EPSILON) continue;

    const Pt centroid = polygon.calculateCentroid();
    const double distance = IRL::magnitude(centroid);
    if (distance > delta_m) {
      continue;
    }
    const double r = distance / delta_m;
    const double distance_weight =
        (1.0 + 4.0 * r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
    const double weight = distance_weight;
    weight_sum += weight;
  }

  // non-dimensionalization
  I = I / (weight_sum * dx * dx);
  return I;
}

const double Jibben_3D::getNormalMetric(void) const {
  const UnsignedIndex_t n_polygons = neighborhood_m->size();

  // normal of target plic in local frame
  const IRL::Normal target_normal(0., 0., 1.);

  // looping over neighboring polygons
  double I = 0.0;
  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const Normal& normal = polygon.getPlaneOfExistence().normal();

    I += 0.5 * (1.0 - std::abs(target_normal * normal));
  }

  return I;
}

const double Jibben_3D::getNormalScatterMetric(void) const {
  const UnsignedIndex_t n_polygons = neighborhood_m->size();

  Eigen::Matrix3d C = Eigen::Matrix3d::Zero();

  UnsignedIndex_t count = 0;
  for (UnsignedIndex_t n = 0; n < n_polygons; ++n) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const Normal& normal = polygon.getPlaneOfExistence().normal();

    Eigen::Vector3d ni(normal[0], normal[1], normal[2]);
    const double norm = ni.norm();
    if (norm <= std::numeric_limits<double>::epsilon()) continue;

    ni /= norm;
    C.noalias() += ni * ni.transpose();
    ++count;
  }

  if (count == 0) return 1.0;

  C /= static_cast<double>(count);  //  make trace(C)=1

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
  if (es.info() != Eigen::Success) return 1.0;

  const Eigen::Vector3d evals = es.eigenvalues();
  const double lambda_max = evals[2];

  // Scatter energy: 0 (perfectly coherent) -> 2/3 (isotropic)
  const double E = 1.0 - lambda_max;

  return E;
}

const double Jibben_3D::getAngularVariance(void) const {
  const UnsignedIndex_t n_polygons = neighborhood_m->size();
  double mean = 0.0;
  std::vector<double> thetas(n_polygons);

  for (UnsignedIndex_t n = 0; n < n_polygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    Normal normal = polygon.getPlaneOfExistence().normal();
    normal.normalize();
    // angle between neighboring and target normal
    double c = normal * IRL::Normal(0., 0., 1.);
    c = std::abs(c);
    c = std::clamp(c, 0.0, 1.0);
    const double theta = std::acos(c);
    thetas[n] = theta;
    mean += theta;
  }
  mean /= static_cast<double>(n_polygons);

  // computing standard deviation
  double var_theta = 0.0;
  for (const auto& theta : thetas) {
    var_theta += (theta - mean) * (theta - mean);
  }
  var_theta /= static_cast<double>(thetas.size());
  double std_theta = std::sqrt(var_theta);

  return std_theta;
}

Paraboloid Jibben_3D::solve2(const JibbenNeighborhood* a_neighborhood_pointer,
                             const double a_delta) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_m = a_neighborhood_pointer;

  Eigen::Matrix<double, 6, 6> Amat = Eigen::Matrix<double, 6, 6>::Zero();
  Eigen::Vector<double, 6> dvec = Eigen::Vector<double, 6>::Zero();

  const double delta = neighborhood_m->getDelta();
  if (a_delta < 0.0) {
    delta_m = delta;
  } else {
    delta_m = a_delta;
  }

  const UnsignedIndex_t npolygons = neighborhood_m->size();
  // std::cout << "Number of polygons = " << npolygons << std::endl;

  // looping over polygons
  for (UnsignedIndex_t n = 0; n < npolygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const double input_weight = neighborhood_m->getWeight(n);

    const UnsignedIndex_t nvertices = polygon.getNumberOfVertices();
    if (nvertices == 0) continue;

    const Pt centroid = polygon.calculateCentroid();
    const Normal& normal = polygon.getPlaneOfExistence().normal();

    if (normal[2] < DBL_EPSILON) continue;

    // distance weight
    const double distance = IRL::magnitude(centroid);
    if (distance > delta_m) {
      continue;
    }
    const double r = distance / delta_m;
    const double distance_weight =
        (1.0 + 4.0 * r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
    double weight = input_weight * distance_weight;

    // weight = 1.0;

    // parameters of the plane
    Eigen::Vector<double, 3> b = Eigen::Vector<double, 3>::Zero();
    b(0) = (centroid * normal) / normal[2];
    b(1) = -normal[0] / normal[2];
    b(2) = -normal[1] / normal[2];

    // calculating monomial integrals by looping over vertices
    Eigen::Vector<double, 15> integrals = Eigen::Vector<double, 15>::Zero();

    for (UnsignedIndex_t v = 0; v < nvertices; v++) {
      const UnsignedIndex_t vn = (v + 1) % nvertices;
      const double xv = polygon[v][0], yv = polygon[v][1];
      const double xvn = polygon[vn][0], yvn = polygon[vn][1];
      const double dxv = xvn - xv, dyv = yvn - yv;

      // ∫dA
      integrals(0) += weight * (xv * yvn - xvn * yv) / 2.0;

      // ∫x dA
      integrals(1) += weight * (xv + xvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫y dA
      integrals(2) += weight * (yv + yvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫x^2 dA
      integrals(3) +=
          weight * (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0;

      // ∫xy dA
      integrals(4) +=
          weight * (yvn - yv) *
          (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
           2.0 * xv * xvn * yvn + xvn * xvn * yv + 3.0 * xvn * xvn * yvn) /
          24.0;

      // ∫y^2 dA
      integrals(5) +=
          weight * (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;

      // ∫x^3 dA
      integrals(6) +=
          weight *
          (dyv * ((dxv * dxv * dxv * dxv) / 5. + dxv * dxv * dxv * xv +
                  2. * (dxv * dxv) * (xv * xv) + 2. * dxv * (xv * xv * xv) +
                  xv * xv * xv * xv)) /
          4.;

      // ∫x^2 y dA
      integrals(7) +=
          weight *
          (dyv * (5. * yv * (dxv + 2. * xv) *
                      (dxv * dxv + 2. * dxv * xv + 2. * (xv * xv)) +
                  dyv * (4. * (dxv * dxv * dxv) + 15. * (dxv * dxv) * xv +
                         20. * dxv * (xv * xv) + 10. * (xv * xv * xv)))) /
          60.;

      // ∫x y^2 dA
      integrals(8) += weight * -0.016666666666666666 *
                      (dxv * (10. * (yv * yv * yv) * (dxv + 2. * xv) +
                              10. * dyv * (yv * yv) * (2. * dxv + 3. * xv) +
                              5. * (dyv * dyv) * yv * (3. * dxv + 4. * xv) +
                              dyv * dyv * dyv * (4. * dxv + 5. * xv)));

      // ∫y^3 dA
      integrals(9) +=
          weight * -0.25 *
          (dxv * ((dyv * dyv * dyv * dyv) / 5. + dyv * dyv * dyv * yv +
                  2. * (dyv * dyv) * (yv * yv) + 2. * dyv * (yv * yv * yv) +
                  yv * yv * yv * yv));

      // ∫x^4 dA
      integrals(10) +=
          weight *
          (dyv *
           ((dxv * dxv * dxv * dxv * dxv) / 6. + dxv * dxv * dxv * dxv * xv +
            (5. * (dxv * dxv * dxv) * (xv * xv)) / 2. +
            (10. * (dxv * dxv) * (xv * xv * xv)) / 3. +
            (5. * dxv * (xv * xv * xv * xv)) / 2. + xv * xv * xv * xv * xv)) /
          5.;

      // ∫x^3 y dA
      integrals(11) +=
          weight *
          (dyv * (dxv * dxv * dxv * dxv * (5. * dyv + 6. * yv) +
                  6. * (dxv * dxv * dxv) * (4. * dyv + 5. * yv) * xv +
                  15. * (dxv * dxv) * (3. * dyv + 4. * yv) * (xv * xv) +
                  20. * dxv * (2. * dyv + 3. * yv) * (xv * xv * xv) +
                  15. * (dyv + 2. * yv) * (xv * xv * xv * xv))) /
          120.;

      // ∫x^2 y^2 dA
      integrals(12) +=
          weight * -0.005555555555555556 *
          (dxv * (20. * (yv * yv * yv) *
                      (dxv * dxv + 3. * dxv * xv + 3. * (xv * xv)) +
                  15. * dyv * (yv * yv) *
                      (3. * (dxv * dxv) + 8. * dxv * xv + 6. * (xv * xv)) +
                  6. * (dyv * dyv) * yv *
                      (6. * (dxv * dxv) + 15. * dxv * xv + 10. * (xv * xv)) +
                  dyv * dyv * dyv *
                      (10. * (dxv * dxv) + 24. * dxv * xv + 15. * (xv * xv))));

      // ∫x y^3 dA
      integrals(13) +=
          weight * -0.008333333333333333 *
          (dxv * (15. * (yv * yv * yv * yv) * (dxv + 2. * xv) +
                  20. * dyv * (yv * yv * yv) * (2. * dxv + 3. * xv) +
                  15. * (dyv * dyv) * (yv * yv) * (3. * dxv + 4. * xv) +
                  6. * (dyv * dyv * dyv) * yv * (4. * dxv + 5. * xv) +
                  dyv * dyv * dyv * dyv * (5. * dxv + 6. * xv)));

      // ∫y^4 dA
      integrals(14) +=
          weight * -0.2 *
          (dxv *
           ((dyv * dyv * dyv * dyv * dyv) / 6. + dyv * dyv * dyv * dyv * yv +
            (5. * (dyv * dyv * dyv) * (yv * yv)) / 2. +
            (10. * (dyv * dyv) * (yv * yv * yv)) / 3. +
            (5. * dyv * (yv * yv * yv * yv)) / 2. + yv * yv * yv * yv * yv));
    }

    // updating matrix and vector to compute paraboloid coefficient
    Amat(0, 0) += integrals(0);
    Amat(0, 1) += integrals(1);
    Amat(0, 2) += integrals(2);
    Amat(0, 3) += integrals(3);
    Amat(0, 4) += integrals(4);
    Amat(0, 5) += integrals(5);
    Amat(1, 1) += integrals(3);
    Amat(1, 2) += integrals(4);
    Amat(1, 3) += integrals(6);
    Amat(1, 4) += integrals(7);
    Amat(1, 5) += integrals(8);
    Amat(2, 2) += integrals(5);
    Amat(2, 3) += integrals(7);
    Amat(2, 4) += integrals(8);
    Amat(2, 5) += integrals(9);
    Amat(3, 3) += integrals(10);
    Amat(3, 4) += integrals(11);
    Amat(3, 5) += integrals(12);
    Amat(4, 4) += integrals(12);
    Amat(4, 5) += integrals(13);
    Amat(5, 5) += integrals(14);

    dvec(0) += b(0) * integrals(0) + b(1) * integrals(1) + b(2) * integrals(2);
    dvec(1) += b(0) * integrals(1) + b(1) * integrals(3) + b(2) * integrals(4);
    dvec(2) += b(0) * integrals(2) + b(1) * integrals(4) + b(2) * integrals(5);
    dvec(3) += b(0) * integrals(3) + b(1) * integrals(6) + b(2) * integrals(7);
    dvec(4) += b(0) * integrals(4) + b(1) * integrals(7) + b(2) * integrals(8);
    dvec(5) += b(0) * integrals(5) + b(1) * integrals(8) + b(2) * integrals(9);
  }

  Amat(1, 0) = Amat(0, 1);
  Amat(2, 0) = Amat(0, 2);
  Amat(2, 1) = Amat(1, 2);
  Amat(3, 0) = Amat(0, 3);
  Amat(3, 1) = Amat(1, 3);
  Amat(3, 2) = Amat(2, 3);
  Amat(4, 0) = Amat(0, 4);
  Amat(4, 1) = Amat(1, 4);
  Amat(4, 2) = Amat(2, 4);
  Amat(4, 3) = Amat(3, 4);
  Amat(5, 0) = Amat(0, 5);
  Amat(5, 1) = Amat(1, 5);
  Amat(5, 2) = Amat(2, 5);
  Amat(5, 3) = Amat(3, 5);
  Amat(5, 4) = Amat(4, 5);
  // solving for paraboloid coefficients
  const Eigen::VectorXd cvec =
      Amat.completeOrthogonalDecomposition().pseudoInverse() * dvec;

  // std::cout << cvec.transpose() << std::endl;

  // storing paraboloid coefficients
  for (UnsignedIndex_t i = 0; i < 6; i++) {
    coefficients_m[i] = cvec(i);
  }

  const double a = cvec[0], b = cvec[1], c = cvec[2], d = cvec[3], e = cvec[4],
               f = cvec[5];
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

  // curvature checks
  //   std::cout << "a = " << paraboloid.getAlignedParaboloid().a() <<
  //   std::endl; std::cout << "b = " << paraboloid.getAlignedParaboloid().b()
  //   << std::endl;

  return paraboloid;
}

Paraboloid Jibben_3D::solveCubic(
    const JibbenNeighborhood* a_neighborhood_pointer, const double a_delta) {
  assert(a_neighborhood_pointer != nullptr);
  neighborhood_m = a_neighborhood_pointer;

  Eigen::Matrix<double, 10, 10> Amat = Eigen::Matrix<double, 10, 10>::Zero();
  Eigen::Vector<double, 10> dvec = Eigen::Vector<double, 10>::Zero();

  const double delta = neighborhood_m->getDelta();
  if (delta_m < 0.0) {
    delta_m = delta;
  }

  const UnsignedIndex_t npolygons = neighborhood_m->size();

  // looping over polygons
  int polygon_count = 1;
  for (UnsignedIndex_t n = 0; n < npolygons; n++) {
    const auto& polygon = neighborhood_m->getPolygon(n);
    const double input_weight = neighborhood_m->getWeight(n);

    const UnsignedIndex_t nvertices = polygon.getNumberOfVertices();
    if (nvertices == 0) continue;

    const Pt centroid = polygon.calculateCentroid();
    const Normal& normal = polygon.getPlaneOfExistence().normal();

    if (normal[2] < DBL_EPSILON) continue;

    // distance weight
    const double distance = IRL::magnitude(centroid);
    // if (distance > delta_m) {
    //   continue;
    // }
    const double r = distance / delta_m;
    const double distance_weight =
        (1.0 + 4.0 * r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
    double weight = input_weight * distance_weight;

    weight = 1.0;

    // parameters of the plane
    Eigen::Vector<double, 3> b = Eigen::Vector<double, 3>::Zero();
    b(0) = (centroid * normal) / normal[2];
    b(1) = -normal[0] / normal[2];
    b(2) = -normal[1] / normal[2];

    // calculating monomial integrals by looping over vertices
    Eigen::Vector<double, 28> integrals = Eigen::Vector<double, 28>::Zero();
    Eigen::Matrix<double, 7, 7> Mmat = Eigen::Matrix<double, 7, 7>::Zero();

    for (UnsignedIndex_t v = 0; v < nvertices; v++) {
      const UnsignedIndex_t vn = (v + 1) % nvertices;
      const double xv = polygon[v][0], yv = polygon[v][1];
      const double xvn = polygon[vn][0], yvn = polygon[vn][1];
      const double dxv = xvn - xv, dyv = yvn - yv;

      // ∫dA
      Mmat(0, 0) += weight * (xv * yvn - xvn * yv) / 2.0;

      // ∫x dA
      Mmat(1, 0) += weight * (xv + xvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫y dA
      Mmat(0, 1) += weight * (yv + yvn) * (xv * yvn - xvn * yv) / 6.0;

      // ∫x^2 dA
      Mmat(2, 0) +=
          weight * (xv + xvn) * (xv * xv + xvn * xvn) * (yvn - yv) / 12.0;

      // ∫xy dA
      Mmat(1, 1) +=
          weight * (yvn - yv) *
          (3.0 * xv * xv * yv + xv * xv * yvn + 2.0 * xv * xvn * yv +
           2.0 * xv * xvn * yvn + xvn * xvn * yv + 3.0 * xvn * xvn * yvn) /
          24.0;

      // ∫y^2 dA
      Mmat(0, 2) +=
          weight * (xv - xvn) * (yv + yvn) * (yv * yv + yvn * yvn) / 12.0;

      // ∫x^3 dA
      Mmat(3, 0) +=
          weight *
          (dyv * ((dxv * dxv * dxv * dxv) / 5. + dxv * dxv * dxv * xv +
                  2. * (dxv * dxv) * (xv * xv) + 2. * dxv * (xv * xv * xv) +
                  xv * xv * xv * xv)) /
          4.;

      // ∫x^2 y dA
      Mmat(2, 1) +=
          weight *
          (dyv * (5. * yv * (dxv + 2. * xv) *
                      (dxv * dxv + 2. * dxv * xv + 2. * (xv * xv)) +
                  dyv * (4. * (dxv * dxv * dxv) + 15. * (dxv * dxv) * xv +
                         20. * dxv * (xv * xv) + 10. * (xv * xv * xv)))) /
          60.;

      // ∫x y^2 dA
      Mmat(1, 2) += weight * -0.016666666666666666 *
                    (dxv * (10. * (yv * yv * yv) * (dxv + 2. * xv) +
                            10. * dyv * (yv * yv) * (2. * dxv + 3. * xv) +
                            5. * (dyv * dyv) * yv * (3. * dxv + 4. * xv) +
                            dyv * dyv * dyv * (4. * dxv + 5. * xv)));

      // ∫y^3 dA
      Mmat(0, 3) +=
          weight * -0.25 *
          (dxv * ((dyv * dyv * dyv * dyv) / 5. + dyv * dyv * dyv * yv +
                  2. * (dyv * dyv) * (yv * yv) + 2. * dyv * (yv * yv * yv) +
                  yv * yv * yv * yv));

      // ∫x^4 dA
      Mmat(4, 0) +=
          weight *
          (dyv *
           ((dxv * dxv * dxv * dxv * dxv) / 6. + dxv * dxv * dxv * dxv * xv +
            (5. * (dxv * dxv * dxv) * (xv * xv)) / 2. +
            (10. * (dxv * dxv) * (xv * xv * xv)) / 3. +
            (5. * dxv * (xv * xv * xv * xv)) / 2. + xv * xv * xv * xv * xv)) /
          5.;

      // ∫x^3 y dA
      Mmat(3, 1) +=
          weight *
          (dyv * (dxv * dxv * dxv * dxv * (5. * dyv + 6. * yv) +
                  6. * (dxv * dxv * dxv) * (4. * dyv + 5. * yv) * xv +
                  15. * (dxv * dxv) * (3. * dyv + 4. * yv) * (xv * xv) +
                  20. * dxv * (2. * dyv + 3. * yv) * (xv * xv * xv) +
                  15. * (dyv + 2. * yv) * (xv * xv * xv * xv))) /
          120.;

      // ∫x^2 y^2 dA
      Mmat(2, 2) +=
          weight * -0.005555555555555556 *
          (dxv * (20. * (yv * yv * yv) *
                      (dxv * dxv + 3. * dxv * xv + 3. * (xv * xv)) +
                  15. * dyv * (yv * yv) *
                      (3. * (dxv * dxv) + 8. * dxv * xv + 6. * (xv * xv)) +
                  6. * (dyv * dyv) * yv *
                      (6. * (dxv * dxv) + 15. * dxv * xv + 10. * (xv * xv)) +
                  dyv * dyv * dyv *
                      (10. * (dxv * dxv) + 24. * dxv * xv + 15. * (xv * xv))));

      // ∫x y^3 dA
      Mmat(1, 3) +=
          weight * -0.008333333333333333 *
          (dxv * (15. * (yv * yv * yv * yv) * (dxv + 2. * xv) +
                  20. * dyv * (yv * yv * yv) * (2. * dxv + 3. * xv) +
                  15. * (dyv * dyv) * (yv * yv) * (3. * dxv + 4. * xv) +
                  6. * (dyv * dyv * dyv) * yv * (4. * dxv + 5. * xv) +
                  dyv * dyv * dyv * dyv * (5. * dxv + 6. * xv)));

      // ∫y^4 dA
      Mmat(0, 4) +=
          weight * -0.2 *
          (dxv *
           ((dyv * dyv * dyv * dyv * dyv) / 6. + dyv * dyv * dyv * dyv * yv +
            (5. * (dyv * dyv * dyv) * (yv * yv)) / 2. +
            (10. * (dyv * dyv) * (yv * yv * yv)) / 3. +
            (5. * dyv * (yv * yv * yv * yv)) / 2. + yv * yv * yv * yv * yv));

      // ∫x^5 dA
      Mmat(5, 0) += (dyv * ((dxv * dxv * dxv * dxv * dxv * dxv) / 7. +
                            dxv * dxv * dxv * dxv * dxv * xv +
                            3. * (dxv * dxv * dxv * dxv) * (xv * xv) +
                            5. * (dxv * dxv * dxv) * (xv * xv * xv) +
                            5. * (dxv * dxv) * (xv * xv * xv * xv) +
                            3. * dxv * (xv * xv * xv * xv * xv) +
                            xv * xv * xv * xv * xv * xv)) /
                    6.;

      // ∫x^4 y dA
      Mmat(4, 1) +=
          (dyv * (7. * yv * (dxv + 2. * xv) * (dxv * dxv + dxv * xv + xv * xv) *
                      (dxv * dxv + 3. * dxv * xv + 3. * (xv * xv)) +
                  dyv * (6. * (dxv * dxv * dxv * dxv * dxv) +
                         35. * (dxv * dxv * dxv * dxv) * xv +
                         84. * (dxv * dxv * dxv) * (xv * xv) +
                         105. * (dxv * dxv) * (xv * xv * xv) +
                         70. * dxv * (xv * xv * xv * xv) +
                         21. * (xv * xv * xv * xv * xv)))) /
          210.;

      // ∫x^3 y^2 dA
      Mmat(3, 2) +=
          (dyv *
           (dxv * dxv * dxv * dxv *
                (15. * (dyv * dyv) + 35. * dyv * yv + 21. * (yv * yv)) +
            7. * (dxv * dxv * dxv) *
                (10. * (dyv * dyv) + 24. * dyv * yv + 15. * (yv * yv)) * xv +
            21. * (dxv * dxv) *
                (6. * (dyv * dyv) + 15. * dyv * yv + 10. * (yv * yv)) *
                (xv * xv) +
            35. * dxv * (3. * (dyv * dyv) + 8. * dyv * yv + 6. * (yv * yv)) *
                (xv * xv * xv) +
            35. * (dyv * dyv + 3. * dyv * yv + 3. * (yv * yv)) *
                (xv * xv * xv * xv))) /
          420.;

      // ∫x^2 y^3 dA
      Mmat(2, 3) +=
          -0.002380952380952381 *
          (dxv * (35. * (yv * yv * yv * yv) *
                      (dxv * dxv + 3. * dxv * xv + 3. * (xv * xv)) +
                  35. * dyv * (yv * yv * yv) *
                      (3. * (dxv * dxv) + 8. * dxv * xv + 6. * (xv * xv)) +
                  21. * (dyv * dyv) * (yv * yv) *
                      (6. * (dxv * dxv) + 15. * dxv * xv + 10. * (xv * xv)) +
                  7. * (dyv * dyv * dyv) * yv *
                      (10. * (dxv * dxv) + 24. * dxv * xv + 15. * (xv * xv)) +
                  dyv * dyv * dyv * dyv *
                      (15. * (dxv * dxv) + 35. * dxv * xv + 21. * (xv * xv))));

      // ∫x y^4 dA
      Mmat(1, 4) +=
          (dxv * (-(dxv * (6. * (dyv * dyv * dyv * dyv * dyv) +
                           35. * (dyv * dyv * dyv * dyv) * yv +
                           84. * (dyv * dyv * dyv) * (yv * yv) +
                           105. * (dyv * dyv) * (yv * yv * yv) +
                           70. * dyv * (yv * yv * yv * yv) +
                           21. * (yv * yv * yv * yv * yv))) -
                  7. * (dyv + 2. * yv) * (dyv * dyv + dyv * yv + yv * yv) *
                      (dyv * dyv + 3. * dyv * yv + 3. * (yv * yv)) * xv)) /
          210.;

      // ∫y^5 dA
      Mmat(0, 5) += -0.16666666666666666 *
                    (dxv * ((dyv * dyv * dyv * dyv * dyv * dyv) / 7. +
                            dyv * dyv * dyv * dyv * dyv * yv +
                            3. * (dyv * dyv * dyv * dyv) * (yv * yv) +
                            5. * (dyv * dyv * dyv) * (yv * yv * yv) +
                            5. * (dyv * dyv) * (yv * yv * yv * yv) +
                            3. * dyv * (yv * yv * yv * yv * yv) +
                            yv * yv * yv * yv * yv * yv));

      // ∫x^6 dA
      Mmat(6, 0) +=
          (dyv * ((dxv * dxv * dxv * dxv * dxv * dxv * dxv) / 8. +
                  dxv * dxv * dxv * dxv * dxv * dxv * xv +
                  (7. * (dxv * dxv * dxv * dxv * dxv) * (xv * xv)) / 2. +
                  7. * (dxv * dxv * dxv * dxv) * (xv * xv * xv) +
                  (35. * (dxv * dxv * dxv) * (xv * xv * xv * xv)) / 4. +
                  7. * (dxv * dxv) * (xv * xv * xv * xv * xv) +
                  (7. * dxv * (xv * xv * xv * xv * xv * xv)) / 2. +
                  xv * xv * xv * xv * xv * xv * xv)) /
          7.;

      // ∫x^5 y dA
      Mmat(5, 1) +=
          (dyv *
           (dxv * dxv * dxv * dxv * dxv * dxv * (7. * dyv + 8. * yv) +
            8. * (dxv * dxv * dxv * dxv * dxv) * (6. * dyv + 7. * yv) * xv +
            28. * (dxv * dxv * dxv * dxv) * (5. * dyv + 6. * yv) * (xv * xv) +
            56. * (dxv * dxv * dxv) * (4. * dyv + 5. * yv) * (xv * xv * xv) +
            70. * (dxv * dxv) * (3. * dyv + 4. * yv) * (xv * xv * xv * xv) +
            56. * dxv * (2. * dyv + 3. * yv) * (xv * xv * xv * xv * xv) +
            28. * (dyv + 2. * yv) * (xv * xv * xv * xv * xv * xv))) /
          336.;

      // ∫x^4 y^2 dA
      Mmat(4, 2) += (dyv * (28. * (yv * yv) * (dxv + 2. * xv) *
                                (dxv * dxv + dxv * xv + xv * xv) *
                                (dxv * dxv + 3. * dxv * xv + 3. * (xv * xv)) +
                            8. * dyv * yv *
                                (6. * (dxv * dxv * dxv * dxv * dxv) +
                                 35. * (dxv * dxv * dxv * dxv) * xv +
                                 84. * (dxv * dxv * dxv) * (xv * xv) +
                                 105. * (dxv * dxv) * (xv * xv * xv) +
                                 70. * dxv * (xv * xv * xv * xv) +
                                 21. * (xv * xv * xv * xv * xv)) +
                            dyv * dyv *
                                (21. * (dxv * dxv * dxv * dxv * dxv) +
                                 120. * (dxv * dxv * dxv * dxv) * xv +
                                 280. * (dxv * dxv * dxv) * (xv * xv) +
                                 336. * (dxv * dxv) * (xv * xv * xv) +
                                 210. * dxv * (xv * xv * xv * xv) +
                                 56. * (xv * xv * xv * xv * xv)))) /
                    840.;

      // ∫x^3 y^3 dA
      Mmat(3, 3) +=
          -0.0008928571428571428 *
          (dxv * (70. * (yv * yv * yv * yv) * (dxv + 2. * xv) *
                      (dxv * dxv + 2. * dxv * xv + 2. * (xv * xv)) +
                  56. * dyv * (yv * yv * yv) *
                      (4. * (dxv * dxv * dxv) + 15. * (dxv * dxv) * xv +
                       20. * dxv * (xv * xv) + 10. * (xv * xv * xv)) +
                  28. * (dyv * dyv) * (yv * yv) *
                      (10. * (dxv * dxv * dxv) + 36. * (dxv * dxv) * xv +
                       45. * dxv * (xv * xv) + 20. * (xv * xv * xv)) +
                  8. * (dyv * dyv * dyv) * yv *
                      (20. * (dxv * dxv * dxv) + 70. * (dxv * dxv) * xv +
                       84. * dxv * (xv * xv) + 35. * (xv * xv * xv)) +
                  dyv * dyv * dyv * dyv *
                      (35. * (dxv * dxv * dxv) + 120. * (dxv * dxv) * xv +
                       140. * dxv * (xv * xv) + 56. * (xv * xv * xv))));

      // ∫x^2 y^4 dA
      Mmat(2, 4) +=
          (dxv *
           (-(dxv * dxv *
              (21. * (dyv * dyv * dyv * dyv * dyv) +
               120. * (dyv * dyv * dyv * dyv) * yv +
               280. * (dyv * dyv * dyv) * (yv * yv) +
               336. * (dyv * dyv) * (yv * yv * yv) +
               210. * dyv * (yv * yv * yv * yv) +
               56. * (yv * yv * yv * yv * yv))) -
            8. * dxv *
                (6. * (dyv * dyv * dyv * dyv * dyv) +
                 35. * (dyv * dyv * dyv * dyv) * yv +
                 84. * (dyv * dyv * dyv) * (yv * yv) +
                 105. * (dyv * dyv) * (yv * yv * yv) +
                 70. * dyv * (yv * yv * yv * yv) +
                 21. * (yv * yv * yv * yv * yv)) *
                xv -
            28. * (dyv + 2. * yv) * (dyv * dyv + dyv * yv + yv * yv) *
                (dyv * dyv + 3. * dyv * yv + 3. * (yv * yv)) * (xv * xv))) /
          840.;

      // ∫x y^5 dA
      Mmat(1, 5) +=
          -0.002976190476190476 *
          (dxv *
           (28. * (yv * yv * yv * yv * yv * yv) * (dxv + 2. * xv) +
            56. * dyv * (yv * yv * yv * yv * yv) * (2. * dxv + 3. * xv) +
            70. * (dyv * dyv) * (yv * yv * yv * yv) * (3. * dxv + 4. * xv) +
            56. * (dyv * dyv * dyv) * (yv * yv * yv) * (4. * dxv + 5. * xv) +
            28. * (dyv * dyv * dyv * dyv) * (yv * yv) * (5. * dxv + 6. * xv) +
            8. * (dyv * dyv * dyv * dyv * dyv) * yv * (6. * dxv + 7. * xv) +
            dyv * dyv * dyv * dyv * dyv * dyv * (7. * dxv + 8. * xv)));

      // ∫y^6 dA
      Mmat(0, 6) +=
          -0.14285714285714285 *
          (dxv * ((dyv * dyv * dyv * dyv * dyv * dyv * dyv) / 8. +
                  dyv * dyv * dyv * dyv * dyv * dyv * yv +
                  (7. * (dyv * dyv * dyv * dyv * dyv) * (yv * yv)) / 2. +
                  7. * (dyv * dyv * dyv * dyv) * (yv * yv * yv) +
                  (35. * (dyv * dyv * dyv) * (yv * yv * yv * yv)) / 4. +
                  7. * (dyv * dyv) * (yv * yv * yv * yv * yv) +
                  (7. * dyv * (yv * yv * yv * yv * yv * yv)) / 2. +
                  yv * yv * yv * yv * yv * yv * yv));
    }

    // basis functions for cubic surface
    const std::pair<int, int> basis[10] = {
        {0, 0},  // 1
        {1, 0},  // x
        {0, 1},  // y
        {2, 0},  // x^2
        {1, 1},  // xy
        {0, 2},  // y^2
        {3, 0},  // x^3
        {2, 1},  // x^2 y
        {1, 2},  // x y^2
        {0, 3}   // y^3
    };

    // updating A matrix
    for (int i = 0; i < 10; i++) {
      auto [i1, j1] = basis[i];
      for (int j = 0; j < 10; j++) {
        auto [i2, j2] = basis[j];
        Amat(i, j) += Mmat(i1 + i2, j1 + j2);
      }
    }

    // updating d vector
    for (int i = 0; i < 10; i++) {
      auto [i1, j1] = basis[i];
      dvec(i) += b(0) * Mmat(i1, j1) + b(1) * Mmat(i1 + 1, j1) +
                 b(2) * Mmat(i1, j1 + 1);
    }
  }

  // symmetry update
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < i; j++) {
      Amat(i, j) = Amat(j, i);
    }
  }

  // solving for paraboloid coefficients
  const Eigen::VectorXd cvec =
      Amat.completeOrthogonalDecomposition().pseudoInverse() * dvec;

  //   std::cout << cvec.transpose() << std::endl;

  // implicit definition of cubic surface
  auto cubic_surface = [cvec](const Eigen::Vector3d& p) -> double {
    const double x = p(0), y = p(1), z = p(2);
    return z -
           (cvec(0) + cvec(1) * x + cvec(2) * y + cvec(3) * x * x +
            cvec(4) * x * y + cvec(5) * y * y + cvec(6) * x * x * x +
            cvec(7) * x * x * y + cvec(8) * x * y * y + cvec(9) * y * y * y);
  };

  // gradient of cubic surface
  auto cubic_gradient = [cvec](const Eigen::Vector3d& p) -> Eigen::Vector3d {
    const double x = p(0), y = p(1), z = p(2);
    return Eigen::Vector3d(
        -(cvec(1) + 2.0 * cvec(3) * x + cvec(4) * y + 3.0 * cvec(6) * x * x +
          2.0 * cvec(7) * x * y + cvec(8) * y * y),
        -(cvec(2) + cvec(4) * x + 2.0 * cvec(5) * y + cvec(7) * x * x +
          2.0 * cvec(8) * x * y + 3.0 * cvec(9) * y * y),
        1.0);
  };

  // hessia of cubic surface
  auto cubic_hessian = [cvec](const Eigen::Vector3d& p) -> Eigen::Matrix3d {
    const double x = p(0), y = p(1), z = p(2);
    Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
    hessian(0, 0) = -(2.0 * cvec(3) + 6.0 * cvec(6) * x + 2.0 * cvec(7) * y);
    hessian(0, 1) = hessian(1, 0) =
        -(cvec(4) + 2.0 * cvec(7) * x + 2.0 * cvec(8) * y);
    hessian(0, 2) = hessian(2, 0) = 0.0;
    hessian(1, 1) = -(2.0 * cvec(5) + 2.0 * cvec(8) * x + 6.0 * cvec(9) * y);
    hessian(1, 2) = hessian(2, 1) = 0.0;
    hessian(2, 2) = 0.0;
    return hessian;
  };

  // projecting plic centroid onto zero level set of cubic surface
  const Pt target_centroid =
      neighborhood_m->getCenterPolygon().calculateCentroid();
  Eigen::Vector3d x_proj(target_centroid[0], target_centroid[1],
                         target_centroid[2]);

  const int max_iterations = 100;
  const double tolerance = 1e-10;

  for (int iter = 0; iter < max_iterations; iter++) {
    const double f_val = cubic_surface(x_proj);
    if (std::abs(f_val) < tolerance) break;
    const Eigen::Vector3d grad = cubic_gradient(x_proj);
    const double grad_norm_sq = grad.squaredNorm();
    if (grad_norm_sq < DBL_EPSILON) break;  // Avoid division by zero
    x_proj -= (f_val / grad_norm_sq) * grad;
    if (iter == max_iterations - 1) {
      std::cerr << "Warning: Cubic projection did not converge after "
                << max_iterations << " iterations." << std::endl;
    }
  }

  // building paraboloid from cubic surface at projected point
  Pt x_proj_pt(x_proj(0), x_proj(1), x_proj(2));
  Eigen::Vector3d grad_at_proj = cubic_gradient(x_proj);
  Eigen::Matrix3d hessian_at_proj = cubic_hessian(x_proj);
  Paraboloid paraboloid_local = IRL::Paraboloid::fromDerivatives(
      x_proj_pt, grad_at_proj, hessian_at_proj);

  // aligning paraboloid with global frame
  const ReferenceFrame paraboloid_local_frame =
      paraboloid_local.getReferenceFrame();
  const Pt paraboloid_local_datum = paraboloid_local.getDatum();
  const Pt local_origin = neighborhood_m->getDatum();
  const ReferenceFrame local_frame = neighborhood_m->getReferenceFrame();

  Pt paraboloid_datum;
  ReferenceFrame paraboloid_frame;

  for (int e = 0; e < 3; e++) {
    for (int d = 0; d < 3; d++) {
      paraboloid_datum[e] += local_frame[d][e] * paraboloid_local_datum[d];
      paraboloid_frame[0][e] +=
          local_frame[d][e] * paraboloid_local_frame[0][d];
      paraboloid_frame[1][e] +=
          local_frame[d][e] * paraboloid_local_frame[1][d];
      paraboloid_frame[2][e] +=
          local_frame[d][e] * paraboloid_local_frame[2][d];
    }
  }
  paraboloid_datum += local_origin;

  IRL::Paraboloid paraboloid(paraboloid_datum, paraboloid_frame,
                             paraboloid_local.getAlignedParaboloid().a(),
                             paraboloid_local.getAlignedParaboloid().b());

  return paraboloid;
}

}  // namespace IRL
