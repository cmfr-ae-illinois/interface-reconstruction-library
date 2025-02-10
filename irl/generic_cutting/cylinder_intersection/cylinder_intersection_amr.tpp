// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_CYLINDER_INTERSECTION_AMR_TPP_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_CYLINDER_INTERSECTION_AMR_TPP_

#include <float.h>
#include <cassert>
#include <cmath>
#include <fstream>
#include "external/NumericalIntegration/NumericalIntegration.h"

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"

namespace IRL {

class ClippedTriangleCylinderVolume_Functor {
 public:
  ClippedTriangleCylinderVolume_Functor(const Pt& p0, const Pt& p1,
                                        const Pt& p2,
                                        const AlignedCylinder& cylinder)
      : x0(p0[0]),
        y0(p0[1]),
        x1(p1[0]),
        y1(p1[1]),
        x2(p2[0]),
        y2(p2[1]),
        b(cylinder.b()),
        r(cylinder.r()) {}

  double operator()(double t) const {
    const double z01 = maximum(
        0.0, r - ((1.0 - t) * y0 + t * y1) * ((1.0 - t) * y0 + t * y1) * b);
    const double z12 = maximum(
        0.0, r - ((1.0 - t) * y1 + t * y2) * ((1.0 - t) * y1 + t * y2) * b);
    const double z20 = maximum(
        0.0, r - ((1.0 - t) * y2 + t * y0) * ((1.0 - t) * y2 + t * y0) * b);
    return (((1.0 - t) * x0 + t * x1) * (y1 - y0) * sqrt(z01) +
            ((1.0 - t) * x1 + t * x2) * (y2 - y1) * sqrt(z12) +
            ((1.0 - t) * x2 + t * x0) * (y0 - y2) * sqrt(z20));
  }

 private:
  const double x0, y0, x1, y1, x2, y2, b, r;
};

template <>
inline Volume computeMomentContributionClippedTriangle<Volume>(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  if (a_print) {
    amr_triangles_clipped.insert(
        amr_triangles_clipped.end(),
        {a_pt_0[0], a_pt_0[1], a_pt_0[2], a_pt_1[0], a_pt_1[1], a_pt_1[2],
         a_pt_2[0], a_pt_2[1], a_pt_2[2]});
  }

  // Using Gauss-Konrod adaptive quadrature
  const double epsabs = 10.0 * DBL_EPSILON;
  const double epsrel = 0.0;
  const int limit = 64;
  // Define the functor
  ClippedTriangleCylinderVolume_Functor functor(a_pt_0, a_pt_1, a_pt_2,
                                                a_aligned_cylinder);
  // Define the integrator.
  Eigen::Integrator<double> integrator(limit);
  // Integrate.
  const double A =
      integrator.quadratureAdaptive(functor, 0.0, 1.0, epsabs, epsrel,
                                    Eigen::Integrator<double>::GaussKronrod21);

  return A;
}

class ClippedTriangleCylinderXCentroid_Functor {
 public:
  ClippedTriangleCylinderXCentroid_Functor(const Pt& p0, const Pt& p1,
                                        const Pt& p2,
                                        const AlignedCylinder& cylinder)
      : x0(p0[0]),
        y0(p0[1]),
        x1(p1[0]),
        y1(p1[1]),
        x2(p2[0]),
        y2(p2[1]),
        b(cylinder.b()),
        r(cylinder.r()) {}

  double operator()(double t) const {
    const double HALF = 1.0 / 2.0;
    const double z01 = maximum(
        0.0, r - ((1.0 - t) * y0 + t * y1) * ((1.0 - t) * y0 + t * y1) * b);
    const double z12 = maximum(
        0.0, r - ((1.0 - t) * y1 + t * y2) * ((1.0 - t) * y1 + t * y2) * b);
    const double z20 = maximum(
        0.0, r - ((1.0 - t) * y2 + t * y0) * ((1.0 - t) * y2 + t * y0) * b);
    const double lx01 = (1.0 - t) * x0 + t * x1;
    const double lx12 = (1.0 - t) * x1 + t * x2;
    const double lx20 = (1.0 - t) * x2 + t * x0;
    return (HALF * lx01 * lx01 * (y1 - y0) * sqrt(z01) +
            HALF * lx12 * lx12 * (y2 - y1) * sqrt(z12) +
            HALF * lx20 * lx20 * (y0 - y2) * sqrt(z20));
  }

 private:
  const double x0, y0, x1, y1, x2, y2, b, r;
};

class ClippedTriangleCylinderYCentroid_Functor {
 public:
  ClippedTriangleCylinderYCentroid_Functor(const Pt& p0, const Pt& p1,
                                        const Pt& p2,
                                        const AlignedCylinder& cylinder)
      : x0(p0[0]),
        y0(p0[1]),
        x1(p1[0]),
        y1(p1[1]),
        x2(p2[0]),
        y2(p2[1]),
        b(cylinder.b()),
        r(cylinder.r()) {}

  double operator()(double t) const {
    const double z01 = maximum(
        0.0, r - ((1.0 - t) * y0 + t * y1) * ((1.0 - t) * y0 + t * y1) * b);
    const double z12 = maximum(
        0.0, r - ((1.0 - t) * y1 + t * y2) * ((1.0 - t) * y1 + t * y2) * b);
    const double z20 = maximum(
        0.0, r - ((1.0 - t) * y2 + t * y0) * ((1.0 - t) * y2 + t * y0) * b);
    const double lx01 = (1.0 - t) * x0 + t * x1;
    const double lx12 = (1.0 - t) * x1 + t * x2;
    const double lx20 = (1.0 - t) * x2 + t * x0;
    const double ly01 = (1.0 - t) * y0 + t * y1;
    const double ly12 = (1.0 - t) * y1 + t * y2;
    const double ly20 = (1.0 - t) * y2 + t * y0;
    return (lx01 * ly01 * (y1 - y0) * sqrt(z01) +
            lx12 * ly12 * (y2 - y1) * sqrt(z12) +
            lx20 * ly20 * (y0 - y2) * sqrt(z20));
  }

 private:
  const double x0, y0, x1, y1, x2, y2, b, r;
};

template <>
inline VolumeMoments computeMomentContributionClippedTriangle<VolumeMoments>(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  if (a_print) {
    amr_triangles_clipped.insert(
        amr_triangles_clipped.end(),
        {a_pt_0[0], a_pt_0[1], a_pt_0[2], a_pt_1[0], a_pt_1[1], a_pt_1[2],
         a_pt_2[0], a_pt_2[1], a_pt_2[2]});
  }

  // Using Gauss-Konrod adaptive quadrature
  const double epsabs = 10.0 * DBL_EPSILON;
  const double epsrel = 0.0;
  const int limit = 64;

  auto moments = VolumeMoments::fromScalarConstant(0.0);

  // Volume

  // Define the functor
  ClippedTriangleCylinderVolume_Functor functor(a_pt_0, a_pt_1, a_pt_2,
                                                a_aligned_cylinder);
  // Define the integrator.
  Eigen::Integrator<double> integrator(limit);
  // Integrate.
  moments.volume() =
      integrator.quadratureAdaptive(functor, 0.0, 1.0, epsabs, epsrel,
                                    Eigen::Integrator<double>::GaussKronrod21);

  // Moments

  // Define the functor
  ClippedTriangleCylinderXCentroid_Functor functorX(a_pt_0, a_pt_1, a_pt_2,
                                                a_aligned_cylinder);
  // Define the integrator.
  Eigen::Integrator<double> integratorX(limit);

  ClippedTriangleCylinderYCentroid_Functor functorY(a_pt_0, a_pt_1, a_pt_2,
                                                a_aligned_cylinder);

  // Define the integrator.
  Eigen::Integrator<double> integratorY(limit);

  // Integrate.
  moments.centroid()[0] =
      integratorX.quadratureAdaptive(functorX, 0.0, 1.0, epsabs, epsrel,
                                    Eigen::Integrator<double>::GaussKronrod21);
  moments.centroid()[1] =
      integratorY.quadratureAdaptive(functorY, 0.0, 1.0, epsabs, epsrel,
                                    Eigen::Integrator<double>::GaussKronrod21);

  // Z component has an simple analitic expression
  moments.centroid()[2] = a_signed_area *
        (6.0 * a_aligned_cylinder.r() -
          (a_pt_0[1] * a_pt_0[1] + a_pt_1[1] * a_pt_1[1] + a_pt_2[1] * a_pt_2[1] +
           a_pt_0[1] * a_pt_1[1] + a_pt_1[1] * a_pt_2[1] + a_pt_2[1] * a_pt_0[1]) *
          a_aligned_cylinder.b()) / 12.0;

  return moments;
}

template <>
inline Volume computeMomentContributionUnclippedTriangle<Volume>(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  if (a_print) {
    amr_triangles_unclipped.insert(
        amr_triangles_unclipped.end(),
        {a_pt_0[0], a_pt_0[1], a_pt_0[2], a_pt_1[0], a_pt_1[1], a_pt_1[2],
         a_pt_2[0], a_pt_2[1], a_pt_2[2]});
  }
  return (a_pt_0[2] + a_pt_1[2] + a_pt_2[2]) * a_signed_area / 3.0;
}

template <>
inline VolumeMoments
computeMomentContributionUnclippedTriangle<VolumeMoments>(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  if (a_print) {
    amr_triangles_unclipped.insert(
        amr_triangles_unclipped.end(),
        {a_pt_0[0], a_pt_0[1], a_pt_0[2], a_pt_1[0], a_pt_1[1], a_pt_1[2],
         a_pt_2[0], a_pt_2[1], a_pt_2[2]});
  }
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  moments.volume() = (a_pt_0[2] + a_pt_1[2] + a_pt_2[2]) * a_signed_area / 3.0;
  moments.centroid()[0] = (2. * a_pt_0[0] * a_pt_0[2] + a_pt_1[0] * a_pt_0[2] +
                           a_pt_2[0] * a_pt_0[2] + a_pt_0[0] * a_pt_1[2] +
                           2. * a_pt_1[0] * a_pt_1[2] + a_pt_2[0] * a_pt_1[2] +
                           a_pt_0[0] * a_pt_2[2] + a_pt_1[0] * a_pt_2[2] +
                           2. * a_pt_2[0] * a_pt_2[2]) *
                          a_signed_area / 12.0;
  moments.centroid()[1] = (2. * a_pt_0[1] * a_pt_0[2] + a_pt_1[1] * a_pt_0[2] +
                           a_pt_2[1] * a_pt_0[2] + a_pt_0[1] * a_pt_1[2] +
                           2. * a_pt_1[1] * a_pt_1[2] + a_pt_2[1] * a_pt_1[2] +
                           a_pt_0[1] * a_pt_2[2] + a_pt_1[1] * a_pt_2[2] +
                           2. * a_pt_2[1] * a_pt_2[2]) *
                          a_signed_area / 12.0;
  moments.centroid()[2] =
      (a_pt_0[2] * a_pt_0[2] + a_pt_0[2] * a_pt_1[2] + a_pt_1[2] * a_pt_1[2] +
       a_pt_0[2] * a_pt_2[2] + a_pt_1[2] * a_pt_2[2] + a_pt_2[2] * a_pt_2[2]) *
      a_signed_area / 12.0;
  return moments;
}

bool vertexBelowCylinderAMR(const Pt& a_pt, const AlignedCylinder& a_cylinder) {
  double z = a_cylinder.r() - a_cylinder.b() * a_pt[1] * a_pt[1];
  if (z < 0.0) return false;
  return z - a_pt[2] * a_pt[2] > 0.0;
}

template <class ReturnType>
void computeMomentContributionMixedTriangle(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const Normal& a_normal,
    const double a_signed_area, std::array<ReturnType, 1>& a_moments_to_add,
    const bool a_print) {
  auto below = StackVector<bool, 3>(
      {vertexBelowCylinderAMR(a_pt_0, a_aligned_cylinder),
       vertexBelowCylinderAMR(a_pt_1, a_aligned_cylinder),
       vertexBelowCylinderAMR(a_pt_2, a_aligned_cylinder)});

  if (below[0] && below[1] && below[2]) {
    a_moments_to_add[0] =
        computeMomentContributionUnclippedTriangle<ReturnType>(
            a_aligned_cylinder, a_pt_0, a_pt_1, a_pt_2, a_signed_area, a_print);
  } else if (!below[0] && !below[1] && !below[2]) {
    a_moments_to_add[0] = computeMomentContributionClippedTriangle<ReturnType>(
        a_aligned_cylinder, a_pt_0, a_pt_1, a_pt_2, a_signed_area, a_print);
  } else {
    const double d = a_normal * a_pt_0;
    auto tris =
        std::array<Pt, 6>({a_pt_0, a_pt_1, a_pt_2, a_pt_0, a_pt_1, a_pt_2});
    for (UnsignedIndex_t v = 0; v < 3; ++v) {
      tris[v][2] = std::sqrt(
          std::max(0.0, a_aligned_cylinder.r() -
                            a_aligned_cylinder.b() * tris[v][1] * tris[v][1]));
    }

    static constexpr std::array<std::array<UnsignedIndex_t, 3>, 3> id_lists{
        {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}}};

    for (UnsignedIndex_t v = 0; v < 3; ++v) {
      const auto& ids = id_lists[v];
      if ((below[ids[0]] && !below[ids[1]] && !below[ids[2]]) ||
          (!below[ids[0]] && below[ids[1]] && below[ids[2]])) {
        const double nx0 = a_normal * tris[ids[0]];
        const double nx1 = a_normal * tris[ids[1]];
        const double nx2 = a_normal * tris[ids[2]];
        if (std::fabs(nx0 - nx1) < DBL_EPSILON ||
            std::fabs(nx0 - nx2) < DBL_EPSILON) {
          Pt center = (tris[ids[0]] + tris[ids[1]] + tris[ids[2]]);
          center /= 3.0;
          if (vertexBelowCylinderAMR(center, a_aligned_cylinder)) {
            a_moments_to_add[0] =
                computeMomentContributionUnclippedTriangle<ReturnType>(
                    a_aligned_cylinder, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
                    a_print);
          } else {
            a_moments_to_add[0] =
                computeMomentContributionClippedTriangle<ReturnType>(
                    a_aligned_cylinder, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
                    a_print);
          }
          break;
        } else {
          double t1 = (nx0 - d) / (nx0 - nx1);
          double t2 = (nx0 - d) / (nx0 - nx2);
          t1 = std::min(1.0, std::max(0.0, t1));
          t2 = std::min(1.0, std::max(0.0, t2));
          Pt new_pt_1 = tris[ids[0]] * (1.0 - t1) + t1 * tris[ids[1]];
          Pt new_pt_2 = tris[ids[0]] * (1.0 - t2) + t2 * tris[ids[2]];
          double f1 = a_aligned_cylinder.b() * new_pt_1[1] * new_pt_1[1] +
                      new_pt_1[2] * new_pt_1[2] - a_aligned_cylinder.r();
          double f2 = a_aligned_cylinder.b() * new_pt_2[1] * new_pt_2[1] +
                      new_pt_2[2] * new_pt_2[2] - a_aligned_cylinder.r();
          // Apply Newton steps
          // if (std::abs(f1) > DBL_EPSILON) {
          //   for (int i = 0; i < 5; i++) {
          //     f1 = a_aligned_cylinder.b() * new_pt_1[1] * new_pt_1[1] +
          //          new_pt_1[2] * new_pt_1[2] - a_aligned_cylinder.r();
          //     double fder1 =
          //         2.0 * a_aligned_cylinder.b() *
          //             (tris[ids[1]][1] - tris[ids[0]][1]) * new_pt_1[1] +
          //         2.0 * (tris[ids[1]][2] - tris[ids[0]][2]) * new_pt_1[2];
          //     t1 -= f1 / safelyTiny(fder1);
          //     t1 = std::min(1.0, std::max(0.0, t1));
          //     new_pt_1 = tris[ids[0]] * (1.0 - t1) + t1 * tris[ids[1]];
          //   }
          // }
          // if (std::abs(f2) > DBL_EPSILON) {
          //   for (int i = 0; i < 5; i++) {
          //     f2 = a_aligned_cylinder.b() * new_pt_2[1] * new_pt_2[1] +
          //          new_pt_2[2] * new_pt_2[2] - a_aligned_cylinder.r();
          //     double fder2 =
          //         2.0 * a_aligned_cylinder.b() *
          //             (tris[ids[2]][1] - tris[ids[0]][1]) * new_pt_2[1] +
          //         2.0 * (tris[ids[2]][2] - tris[ids[0]][2]) * new_pt_2[2];
          //     t2 -= f2 / safelyTiny(fder2);
          //     t2 = std::min(1.0, std::max(0.0, t2));
          //     new_pt_2 = tris[ids[0]] * (1.0 - t2) + t2 * tris[ids[2]];
          //   }
          // }
          const double signed_area_1 =
              triangleSignedArea(tris[ids[0]], new_pt_1, new_pt_2);
          const double signed_area_2 =
              triangleSignedArea(new_pt_1, tris[ids[1]], tris[ids[2]]);
          const double signed_area_3 =
              triangleSignedArea(tris[ids[2]], new_pt_2, new_pt_1);
          if (below[ids[0]]) {
            a_moments_to_add[0] =
                computeMomentContributionUnclippedTriangle<ReturnType>(
                    a_aligned_cylinder, tris[3 + ids[0]], new_pt_1, new_pt_2,
                    signed_area_1, a_print) +
                computeMomentContributionClippedTriangle<ReturnType>(
                    a_aligned_cylinder, new_pt_1, tris[3 + ids[1]],
                    tris[3 + ids[2]], signed_area_2, a_print) +
                computeMomentContributionClippedTriangle<ReturnType>(
                    a_aligned_cylinder, tris[3 + ids[2]], new_pt_2, new_pt_1,
                    signed_area_3, a_print);
          } else {
            a_moments_to_add[0] =
                computeMomentContributionClippedTriangle<ReturnType>(
                    a_aligned_cylinder, tris[3 + ids[0]], new_pt_1, new_pt_2,
                    signed_area_1, a_print) +
                computeMomentContributionUnclippedTriangle<ReturnType>(
                    a_aligned_cylinder, new_pt_1, tris[3 + ids[1]],
                    tris[3 + ids[2]], signed_area_2, a_print) +
                computeMomentContributionUnclippedTriangle<ReturnType>(
                    a_aligned_cylinder, tris[3 + ids[2]], new_pt_2, new_pt_1,
                    signed_area_3, a_print);
          }
          break;
        }
      }
    }
  }
}

std::pair<bool, bool> computeZBounds(const AlignedCylinder& a_aligned_cylinder,
                                     const Pt& a_pt_0, const Pt& a_pt_1,
                                     const Pt& a_pt_2) {
  // Compute bounding box of triangle
  std::array<double, 6> tri_bounds;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    tri_bounds[2 * d] = std::min({a_pt_0[d], a_pt_1[d], a_pt_2[d]});
    tri_bounds[2 * d + 1] = std::max({a_pt_0[d], a_pt_1[d], a_pt_2[d]});
  }
  // Compute bouding box of cylinder (over-estimate)
  double z0 = a_aligned_cylinder.r() -
              a_aligned_cylinder.b() * tri_bounds[2] * tri_bounds[2];
  double z1 = a_aligned_cylinder.r() -
              a_aligned_cylinder.b() * tri_bounds[3] * tri_bounds[3];
  z0 = std::sqrt(std::max(0.0, z0));
  z1 = std::sqrt(std::max(0.0, z1));
  if (a_aligned_cylinder.b() >= 0) {
    const double cyl_max = tri_bounds[2] * tri_bounds[3] <= 0.0
                              ? std::sqrt(a_aligned_cylinder.r())
                              : std::max({z0, z1});
    const double cyl_min = std::min({z0, z1});
    // Return min/max of z height of triangle and cylinder
    return std::pair<bool, bool>({cyl_min > tri_bounds[5] + AMR_DBL_EPSILON,
                                  tri_bounds[4] > cyl_max + AMR_DBL_EPSILON});
  } else {
    const double cyl_max = std::max({z0, z1});
    const double cyl_min = tri_bounds[2] * tri_bounds[3] <= 0.0
                              ? std::sqrt(a_aligned_cylinder.r())
                              : std::min({z0, z1});
    // Return min/max of z height of triangle and cylinder
    return std::pair<bool, bool>({cyl_min > tri_bounds[5] + AMR_DBL_EPSILON,
                                  tri_bounds[4] > cyl_max + AMR_DBL_EPSILON});
  }
}

template <class ReturnType>
void computeMomentContributionAMR(
    const AlignedCylinder& a_aligned_cylinder, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const Normal& a_normal,
    const double a_signed_area, const UnsignedIndex_t a_amr_level,
    const UnsignedIndex_t a_max_amr_level,
    std::array<std::pair<ReturnType, ReturnType>, 1>& a_full_moments_ref,
    std::array<std::pair<ReturnType, ReturnType>, 1>& a_full_moments,
    const bool a_print) {
  std::array<ReturnType, 1> moments_to_add;

  // Compute z-bounds of triangle and cylinder
  auto z_limits = computeZBounds(a_aligned_cylinder, a_pt_0, a_pt_1, a_pt_2);

  if (z_limits.first) {
    // Max of triangle is smaller than min of paraboloid
    moments_to_add[0] = computeMomentContributionUnclippedTriangle<ReturnType>(
        a_aligned_cylinder, a_pt_0, a_pt_1, a_pt_2, a_signed_area, a_print);
    kahanSummationMoments<ReturnType>(a_full_moments, a_full_moments_ref,
                                      moments_to_add);
    return;
  } else if (z_limits.second) {
    // Max of triangle is smaller than min of paraboloid
    moments_to_add[0] = computeMomentContributionClippedTriangle<ReturnType>(
        a_aligned_cylinder, a_pt_0, a_pt_1, a_pt_2, a_signed_area, a_print);
    kahanSummationMoments<ReturnType>(a_full_moments, a_full_moments_ref,
                                      moments_to_add);
    return;
  } else if (a_amr_level == a_max_amr_level) {
    computeMomentContributionMixedTriangle<ReturnType>(
        a_aligned_cylinder, a_pt_0, a_pt_1, a_pt_2, a_normal, a_signed_area,
        moments_to_add, a_print);
    kahanSummationMoments<ReturnType>(a_full_moments, a_full_moments_ref,
                                      moments_to_add);
    return;
  }
  // Refine triangle and call function recursively
  const Pt c0 = 0.5 * (a_pt_0 + a_pt_1);
  const Pt c1 = 0.5 * (a_pt_1 + a_pt_2);
  const Pt c2 = 0.5 * (a_pt_2 + a_pt_0);
  computeMomentContributionAMR<ReturnType>(
      a_aligned_cylinder, a_pt_0, c0, c2, a_normal, 0.25 * a_signed_area,
      a_amr_level + 1, a_max_amr_level, a_full_moments_ref, a_full_moments,
      a_print);
  computeMomentContributionAMR<ReturnType>(
      a_aligned_cylinder, a_pt_1, c1, c0, a_normal, 0.25 * a_signed_area,
      a_amr_level + 1, a_max_amr_level, a_full_moments_ref, a_full_moments,
      a_print);
  computeMomentContributionAMR<ReturnType>(
      a_aligned_cylinder, a_pt_2, c2, c1, a_normal, 0.25 * a_signed_area,
      a_amr_level + 1, a_max_amr_level, a_full_moments_ref, a_full_moments,
      a_print);
  computeMomentContributionAMR<ReturnType>(
      a_aligned_cylinder, c0, c1, c2, a_normal, 0.25 * a_signed_area,
      a_amr_level + 1, a_max_amr_level, a_full_moments_ref, a_full_moments,
      a_print);
}

template <class SegmentedHalfEdgePolyhedronType>
void printClippedPolytope(SegmentedHalfEdgePolyhedronType* a_polytope) {
  const auto number_of_faces = a_polytope->getNumberOfFaces();
  for (UnsignedIndex_t f = 0; f < number_of_faces; ++f) {
    const auto& face = *(*a_polytope)[f];
    const auto starting_half_edge = face.getStartingHalfEdge();
    const auto& ref_pt = starting_half_edge->getVertex()->getLocation().getPt();
    auto current_half_edge =
        starting_half_edge->getNextHalfEdge()->getNextHalfEdge();
    auto prev_pt =
        current_half_edge->getPreviousVertex()->getLocation().getPt();
    do {
      const auto& curr_pt =
          current_half_edge->getVertex()->getLocation().getPt();
      amr_triangles_clipped.insert(
          amr_triangles_clipped.end(),
          {ref_pt[0], ref_pt[1], ref_pt[2], prev_pt[0], prev_pt[1], prev_pt[2],
           curr_pt[0], curr_pt[1], curr_pt[2]});
      prev_pt = curr_pt;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != starting_half_edge);
  }
}

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithCylinderAMR(SegmentedHalfEdgePolyhedronType* a_polytope,
                                   HalfEdgePolytopeType* a_complete_polytope,
                                   const AlignedCylinder& a_cylinder,
                                   const UnsignedIndex_t a_max_amr_level,
                                   const std::string& a_filename) {
  using ReturnScalarType = typename ReturnType::value_type;
  const bool print = !a_filename.empty();

  // We have 1 strategies for computing the moments (and will later on choose
  // the best of those)
  std::array<std::pair<ReturnType, ReturnType>, 1> full_moments;
  std::array<std::pair<ReturnType, ReturnType>, 1> full_moments_ref;
  full_moments[0] =
      std::pair<ReturnType, ReturnType>({ReturnType::fromScalarConstant(0.0),
                                         ReturnType::fromScalarConstant(0.0)});
  full_moments_ref[0] =
      std::pair<ReturnType, ReturnType>({ReturnType::fromScalarConstant(0.0),
                                         ReturnType::fromScalarConstant(0.0)});

  bool elliptic = a_cylinder.b() >= 0.0;

  int nb_cut = elliptic ? 4 : 2;

  // If elliptic cylinder, clip region outside cylinder
  if (elliptic) {
    const double semi_major_axis = std::sqrt(a_cylinder.r());
    const double semi_minor_axis = std::sqrt(a_cylinder.r() / a_cylinder.b());
    SegmentedHalfEdgePolyhedronType dummy_clipped_polytope;
    splitHalfEdgePolytope(a_polytope, &dummy_clipped_polytope,
                          a_complete_polytope,
                          Plane(Normal(0.0, 1.0, 0.0), semi_minor_axis));
    if (print) printClippedPolytope(&dummy_clipped_polytope);
    splitHalfEdgePolytope(a_polytope, &dummy_clipped_polytope,
                          a_complete_polytope,
                          Plane(Normal(0.0, -1.0, 0.0), semi_minor_axis));
    if (print) printClippedPolytope(&dummy_clipped_polytope);
    splitHalfEdgePolytope(a_polytope, &dummy_clipped_polytope,
                          a_complete_polytope,
                          Plane(Normal(0.0, 0.0, 1.0), semi_major_axis));
    if (print) printClippedPolytope(&dummy_clipped_polytope);
    splitHalfEdgePolytope(a_polytope, &dummy_clipped_polytope,
                          a_complete_polytope,
                          Plane(Normal(0.0, 0.0, -1.0), semi_major_axis));
    if (print) printClippedPolytope(&dummy_clipped_polytope);
  }

  // Split polytope into Z+ and Z- contributions
  // SegmentedHalfEdgePolyhedronType bottom_polytope;
  // splitHalfEdgePolytope(a_polytope, &bottom_polytope, a_complete_polytope,
  //                       Plane(Normal(0.0, 0.0, -1.0), 0.0));

  std::vector<SegmentedHalfEdgePolyhedronType*> polytope_list;
  std::vector<AlignedCylinder> cylinder_list;
  std::vector<double> rotation_list;
  SegmentedHalfEdgePolyhedronType p1, p2, p3;

  if (elliptic) {
    const auto rotated_cylinder =
        AlignedCylinder({1.0 / a_cylinder.b(), a_cylinder.r() / a_cylinder.b()});
    splitHalfEdgePolytope(
        a_polytope, &p1, a_complete_polytope,
        Plane(Normal(0.0, 1.0 / std::sqrt(2.0), -1.0 / std::sqrt(2.0)), 0.0));
    splitHalfEdgePolytope(
        a_polytope, &p2, a_complete_polytope,
        Plane(Normal(0.0, -1.0 / std::sqrt(2.0), -1.0 / std::sqrt(2.0)), 0.0));
    polytope_list.push_back(a_polytope);
    cylinder_list.push_back(a_cylinder);
    rotation_list.push_back(0.0);
    // SegmentedHalfEdgePolyhedronType* N_polytope = a_polytope;
    polytope_list.push_back(&p2);
    cylinder_list.push_back(rotated_cylinder);
    rotation_list.push_back(3.0 * M_PI / 2.0);
    // SegmentedHalfEdgePolyhedronType* W_polytope = &p2;
    splitHalfEdgePolytope(
        &p1, &p3, a_complete_polytope,
        Plane(Normal(0.0, -1.0 / std::sqrt(2.0), -1.0 / std::sqrt(2.0)), 0.0));
    polytope_list.push_back(&p1);
    cylinder_list.push_back(rotated_cylinder);
    rotation_list.push_back(M_PI / 2.0);
    // SegmentedHalfEdgePolyhedronType* E_polytope = &p1;
    polytope_list.push_back(&p3);
    cylinder_list.push_back(a_cylinder);
    rotation_list.push_back(M_PI);
    // SegmentedHalfEdgePolyhedronType* S_polytope = &p3;

    // std::array<SegmentedHalfEdgePolyhedronType*, 4> polytope_list = {
    //     N_polytope, E_polytope, S_polytope, W_polytope};
    // std::array<AlignedCylinder, 4> cylinder_list = {a_cylinder, rotated_cylinder,
    //                                                 a_cylinder, rotated_cylinder};
    // std::array<double, 4> rotation_list = {0.0, M_PI / 2.0, M_PI,
    //                                       3.0 * M_PI / 2.0};
  } else {
    splitHalfEdgePolytope(
        a_polytope, &p1, a_complete_polytope,
        Plane(Normal(0.0, 0.0, -1.0), 0.0));
    polytope_list.push_back(a_polytope);
    cylinder_list.push_back(a_cylinder);
    rotation_list.push_back(0.0);
    polytope_list.push_back(&p1);
    cylinder_list.push_back(a_cylinder);
    rotation_list.push_back(M_PI);
  }

  for (int i = 0; i < nb_cut; i++) {
    auto& polytope = polytope_list[i];
    auto& cylinder = cylinder_list[i];
    if (polytope != nullptr) {
      // Rotate polytope
      ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
      UnitQuaternion x_rotation(-rotation_list[i], frame[0]);
      x_rotation.normalize();
      frame = x_rotation * frame;
      // also need to rotate the centroid
      if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
        // std::cout << "hmm " << full_moments[0].first.centroid() << std::endl;
        Pt centroid01(full_moments[0].first.centroid()[0], 
          full_moments[0].first.centroid()[1], full_moments[0].first.centroid()[2]);
        Pt centroid02(full_moments[0].second.centroid()[0],
          full_moments[0].second.centroid()[1], full_moments[0].second.centroid()[2]);
        Pt centroid11(full_moments[1].first.centroid()[0],
          full_moments[1].first.centroid()[1], full_moments[1].first.centroid()[2]);
        Pt centroid12(full_moments[1].second.centroid()[0],
          full_moments[1].second.centroid()[1], full_moments[1].second.centroid()[2]);
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          full_moments[0].first.centroid()[n] = frame[n] * centroid01;
          full_moments[0].second.centroid()[n] = frame[n] * centroid02;
          full_moments[1].first.centroid()[n] = frame[n] * centroid11;
          full_moments[1].second.centroid()[n] = frame[n] * centroid12;
        }
      }
      for (UnsignedIndex_t v = 0; v < polytope->getNumberOfVertices(); ++v) {
        const Pt original_pt = polytope->getVertex(v)->getLocation().getPt();
        Pt pt(0, 0, 0);
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          pt[n] = frame[n] * original_pt;
        }
        polytope->getVertex(v)->setLocation(pt);
      }

      // The AMR moment calculation needs plane information
      const auto number_of_faces = polytope->getNumberOfFaces();
      for (UnsignedIndex_t f = 0; f < number_of_faces; ++f) {
        auto& face = *(*polytope)[f];
        auto normal = Normal(0, 0, 0);
        const auto starting_half_edge = face.getStartingHalfEdge();
        auto current_half_edge = starting_half_edge;
        auto next_half_edge = starting_half_edge->getNextHalfEdge();
        const auto& start_location =
            starting_half_edge->getPreviousVertex()->getLocation();
        do {
          normal += crossProduct(
              current_half_edge->getVertex()->getLocation() - start_location,
              next_half_edge->getVertex()->getLocation() - start_location);
          current_half_edge = next_half_edge;
          next_half_edge = next_half_edge->getNextHalfEdge();
        } while (next_half_edge != starting_half_edge);
        normal.normalize();
        face.setPlane(Plane(normal, normal * start_location));
      }

      int tri_clipped_size = amr_triangles_clipped.size();
      int tri_unclipped_size = amr_triangles_unclipped.size();

      // We loop over each face, triangulate, and split until there are no
      // intersections OR we have reached the max refinement level
      for (UnsignedIndex_t f = 0; f < number_of_faces; ++f) {
        const auto& face = *(*polytope)[f];
        const auto& face_normal = face.getPlane().normal();
        const auto starting_half_edge = face.getStartingHalfEdge();
        const auto& ref_pt =
            starting_half_edge->getVertex()->getLocation().getPt();
        auto current_half_edge =
            starting_half_edge->getNextHalfEdge()->getNextHalfEdge();
        auto prev_pt =
            current_half_edge->getPreviousVertex()->getLocation().getPt();
        do {
          const auto& curr_pt =
              current_half_edge->getVertex()->getLocation().getPt();
          const double signed_area =
              triangleSignedArea(ref_pt, prev_pt, curr_pt);
          computeMomentContributionAMR<ReturnType>(
              cylinder, ref_pt, prev_pt, curr_pt, face_normal, signed_area, 0,
              a_max_amr_level, full_moments_ref, full_moments, print);
          prev_pt = curr_pt;
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
      }

      // Revert rotate polytope
      frame = ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
      x_rotation = UnitQuaternion(rotation_list[i], frame[0]);
      x_rotation.normalize();
      frame = x_rotation * frame;
      if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
        Pt centroid01(full_moments[0].first.centroid()[0], 
          full_moments[0].first.centroid()[1], full_moments[0].first.centroid()[2]);
        Pt centroid02(full_moments[0].second.centroid()[0],
          full_moments[0].second.centroid()[1], full_moments[0].second.centroid()[2]);
        Pt centroid11(full_moments[1].first.centroid()[0],
          full_moments[1].first.centroid()[1], full_moments[1].first.centroid()[2]);
        Pt centroid12(full_moments[1].second.centroid()[0],
          full_moments[1].second.centroid()[1], full_moments[1].second.centroid()[2]);
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          full_moments[0].first.centroid()[n] = frame[n] * centroid01;
          full_moments[0].second.centroid()[n] = frame[n] * centroid02;
          full_moments[1].first.centroid()[n] = frame[n] * centroid11;
          full_moments[1].second.centroid()[n] = frame[n] * centroid12;
        }
      }
      for (UnsignedIndex_t v = 0; v < polytope->getNumberOfVertices(); ++v) {
        const Pt original_pt = polytope->getVertex(v)->getLocation().getPt();
        Pt pt(0, 0, 0);
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          pt[n] = frame[n] * original_pt;
        }
        polytope->getVertex(v)->setLocation(pt);
      }
      if (print) {
        for (std::size_t i = tri_clipped_size; i < amr_triangles_clipped.size();
             i += 3) {
          Pt pt(amr_triangles_clipped[i], amr_triangles_clipped[i + 1],
                amr_triangles_clipped[i + 2]);
          Pt new_pt(0, 0, 0);
          for (UnsignedIndex_t n = 0; n < 3; ++n) {
            new_pt[n] = frame[n] * pt;
          }
          amr_triangles_clipped[i] = new_pt[0];
          amr_triangles_clipped[i + 1] = new_pt[1];
          amr_triangles_clipped[i + 2] = new_pt[2];
        }
        for (std::size_t i = tri_unclipped_size;
             i < amr_triangles_unclipped.size(); i += 3) {
          Pt pt(amr_triangles_unclipped[i], amr_triangles_unclipped[i + 1],
                amr_triangles_unclipped[i + 2]);
          Pt new_pt(0, 0, 0);
          for (UnsignedIndex_t n = 0; n < 3; ++n) {
            new_pt[n] = frame[n] * pt;
          }
          amr_triangles_unclipped[i] = new_pt[0];
          amr_triangles_unclipped[i + 1] = new_pt[1];
          amr_triangles_unclipped[i + 2] = new_pt[2];
        }
      }
    }
  }

  // Print triangles
  if (print) {
    auto write_triangles = [a_filename](
                               const std::string& a_header,
                               const std::string& a_file_prefix,
                               const std::vector<double>& a_triangle_list) {
      // binary file
      char head[80];
      std::strncpy(head, a_header.c_str(), a_header.size() - 1);
      char attribute[2] = "0";
      char dummy[4] = "0";
      const auto nTriLong = amr_triangles_clipped.size() / 3;

      std::ofstream myfile;

      myfile.open(a_file_prefix + a_filename + ".stl",
                  std::ios::out | std::ios::binary);
      myfile.write(head, sizeof(head));
      myfile.write(reinterpret_cast<const char*>(&nTriLong), 4);

      // write down every triangle
      for (std::size_t i = 0; i < a_triangle_list.size(); i += 9) {
        // normal vector coordinates
        myfile.write(dummy, 4);
        myfile.write(dummy, 4);
        myfile.write(dummy, 4);
        for (std::size_t t = 0; t < 9; ++t) {
          const float val = static_cast<float>(a_triangle_list[i + t]);
          myfile.write(reinterpret_cast<const char*>(&val), 4);
        }
        myfile.write(attribute, 2);
      }
    };
    write_triangles("clipped_triangles", "clipped", amr_triangles_clipped);
    write_triangles("unclipped_triangles", "unclipped",
                    amr_triangles_unclipped);
  }

  amr_triangles_unclipped.clear();
  amr_triangles_clipped.clear();

  return full_moments[0].first;
}
}  // namespace IRL
#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_CYLINDER_INTERSECTION_AMR_TPP_
