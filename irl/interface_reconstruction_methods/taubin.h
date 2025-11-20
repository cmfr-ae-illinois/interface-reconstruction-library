// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_TAUBIN_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_TAUBIN_H_

#include <Eigen/Dense>
#include <cassert>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/interface_reconstruction_methods/jibben_neighborhood.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

class Taubin_3D {
 public:
  /// \brief Default constructor.
  Taubin_3D(void) = default;

  /// \brief Solve the system for the reconstruction
  Paraboloid solve(const JibbenNeighborhood* a_neighborhood_pointer,
                   const double a_delta = -1.0);

  /// \brief Default destructor.
  ~Taubin_3D(void) = default;

 private:
  /// \brief Solve the system for the reconstruction
  Paraboloid solve(void);

  /// \brief Get intersection points from a plane and polygon
  void getIntersectionPts(const IRL::Polygon& a_polygon,
                          const IRL::Plane& a_cutting_plane,
                          IRL::StackVector<IRL::Pt, 2>* a_intersection_pts);

  /// \brief normal weights
  double getNormalWeight(const IRL::Normal& a_nref, const IRL::Normal& a_nloc);

  /// \brief distance weights
  double getDistanceWeight(const IRL::Pt& a_pref, const IRL::Pt& a_ploc);

  /// \brief Calculating polygon normal
  IRL::Normal calculatePolygonNormal(const IRL::Polygon& a_polygon);

  /// \brief reference frame from normal
  IRL::ReferenceFrame referenceFrameFromNormal(const IRL::Normal a_normal);

  /// \brief calculating signed taubin curvature
  double getSignedTaubinCurvature(const std::vector<Eigen::Vector2d>& a_points,
                                  const std::vector<double>& a_weights);

  /// \brief calculating principal directions and curvatures
  std::tuple<double, double, double, bool> getPrincipalCurvatures(
      const std::vector<double>& a_theta, const std::vector<double>& a_k_theta);

  /// \brief sampling points along a slice in rotated frame
  void sampleLocalPoints(
      const IRL::ReferenceFrame& local_frame, const IRL::Pt& local_datum,
      const std::vector<std::pair<IRL::Pt, IRL::Pt>>& end_points_list,
      const int& num_samples_per_segment,
      std::vector<Eigen::Vector2d>* points_local_frame);

  /// \brief Storage of the stencil information
  const JibbenNeighborhood* neighborhood_m;
  /// \brief Weighting function radius
  double delta_m;
};

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_TAUBIN_H_
