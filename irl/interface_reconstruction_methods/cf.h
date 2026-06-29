// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_CIRCLE_FIT_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_CIRCLE_FIT_H_

#include <Eigen/Dense>

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/interface_reconstruction_methods/jibben_neighborhood.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

class CircleFit_3D {
 public:
  /// \brief Default constructor.
  CircleFit_3D(void) = default;

  /// \brief Solve the system for the reconstruction.
  ///
  /// a_h is the mesh spacing. If a_h < 0, this assumes
  /// neighborhood.getDelta() = 2.5*h.
  Paraboloid solve(const JibbenNeighborhood* a_neighborhood_pointer,
                   const double a_h = -1.0);

  /// \brief Default destructor.
  ~CircleFit_3D(void) = default;

 private:
  /// \brief Solve the system for the reconstruction.
  Paraboloid solve(void);

  /// \brief Get intersection points from a plane and polygon.
  void getIntersectionPts(const IRL::Polygon& a_polygon,
                          const IRL::Plane& a_cutting_plane,
                          IRL::StackVector<IRL::Pt, 2>* a_intersection_pts);

  /// \brief Volume-fraction weight.
  double getVfWeight(const double& a_vfrac);

  /// \brief Normal alignment weight.
  double getNormalWeight(const IRL::Normal& a_nref, const IRL::Normal& a_nloc);

  /// \brief Distance weight.
  double getDistanceWeight(const IRL::Pt& a_pref, const IRL::Pt& a_ploc);

  /// \brief Total weight.
  double getTotalWeight(const double& a_vfrac, const IRL::Normal& a_nref,
                        const IRL::Normal& a_nloc, const IRL::Pt& a_pref,
                        const IRL::Pt& a_ploc);

  /// \brief Calculate polygon normal.
  IRL::Normal calculatePolygonNormal(const IRL::Polygon& a_polygon);

  /// \brief Construct reference frame from normal.
  IRL::ReferenceFrame referenceFrameFromNormal(const IRL::Normal a_normal);

  /// \brief Convert a global point to a rotated local 2D slice frame.
  IRL::Pt toLocal(const IRL::ReferenceFrame& a_local_frame,
                  const IRL::Pt& a_local_datum, const IRL::Pt& a_global_pt);

  /// \brief Moment terms for one line segment in the local slice frame.
  std::vector<double> getTaubinMomentTerms(const double& x0, const double& y0,
                                           const double& dx, const double& dy);

  /// \brief Build Taubin moment and constraint matrices.
  void getTaubinMatrices(
      const std::vector<std::pair<IRL::Pt, IRL::Pt>>& a_end_points,
      const std::vector<double>& a_weights, Eigen::Matrix4d* a_M,
      Eigen::Matrix4d* a_C);

  /// \brief Signed curvature from line-segment Taubin circle fit.
  double getSignedTaubinCurvature(
      const std::vector<std::pair<IRL::Pt, IRL::Pt>>& a_end_points,
      const std::vector<double>& a_weights);

  /// \brief Calculate principal curvatures and principal direction angle.
  std::tuple<double, double, double, bool> getPrincipalCurvatures(
      const std::vector<double>& a_theta, const std::vector<double>& a_k_theta);

  /// \brief Storage of the stencil information.
  const JibbenNeighborhood* neighborhood_m = nullptr;

  /// \brief Mesh spacing.
  double h_m = -1.0;

  /// \brief Weighting radius, usually 2.5*h.
  double delta_m = -1.0;
};

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_CIRCLE_FIT_H_