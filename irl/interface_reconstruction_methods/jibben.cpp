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

      // ∫x^3 dA CHANGE
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
  I = I / (projected_area * dx *
           dx);  // CHANGE EITHER SQUARE ROOT NUM OR dx^2 IN DENOM

  // I = std::sqrt(I) / (projected_area * dx);
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

      // ∫x^3 dA CHANGE
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
  I = I / (projected_area * dx *
           dx);  // CHANGE EITHER SQUARE ROOT NUM OR dx^2 IN DENOM

  // I = std::sqrt(I) / (projected_area * dx);
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

      // ∫x^3 dA CHANGE
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
  I = I /
      (weight_sum * dx * dx);  // CHANGE EITHER SQUARE ROOT NUM OR dx^2 IN DENOM

  // I = std::sqrt(I) / (projected_area * dx);
  return I;
}

}  // namespace IRL
