// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_VARIANT_ADVECTOR_RECONSTRUCTION_TYPES_H_
#define EXAMPLES_VARIANT_ADVECTOR_RECONSTRUCTION_TYPES_H_

#include <string>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/interface_reconstruction_methods/plicnet.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/moments/volume_moments.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include "examples/variant_advector/data.h"
#include "examples/variant_advector/vtk.h"

#include "Eigen/Dense"

void getReconstruction(const std::string& a_reconstruction_method,
                       const Data<IRL::VolumeMoments>& a_liq_moments,
                       const Data<IRL::VolumeMoments>& a_gas_moments,
                       const double a_dt, const Data<double>& a_U,
                       const Data<double>& a_V, const Data<double>& a_W,
                       Data<IRL::SeparatorVariant>* a_interface,
                       std::vector<InterfaceScalarField>* a_scalar_fields);

struct ELVIRA {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr);
};

struct LVIRA {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct PLICnet {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct PU {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct Jibben {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false,
      Data<IRL::Pt>* a_centroids = nullptr, Data<double>* a_areas = nullptr,
      Data<double>* a_errors = nullptr);
};

struct iJibben {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false,
      Data<IRL::Pt>* a_centroids = nullptr, Data<double>* a_areas = nullptr,
      Data<double>* a_errors = nullptr);
};

struct Jibben2 {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false,
      Data<IRL::Pt>* a_centroids = nullptr, Data<double>* a_areas = nullptr,
      Data<double>* a_errors = nullptr);
};

struct JibbenCubic {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false,
      Data<IRL::Pt>* a_centroids = nullptr, Data<double>* a_areas = nullptr,
      Data<double>* a_errors = nullptr);
};

struct JibbenM {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr);
};

struct MixedPLICJibben {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr);
};

struct SlicesParabola {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct SlicesTaubin {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct SlicesTaubin2 {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct SlicesTaubin3 {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct SlicesTaubinLS {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct SlicesTaubinS {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct PLICalign {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct PLICalign2 {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct MossoSwartz {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct JibbenPU {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct Testing {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct Hybrid2 {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct JibbenSqPU {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

struct SlicesParticle {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      std::vector<InterfaceScalarField>* a_scalar_fields = nullptr,
      const bool a_plic_already_built = false);
};

const IRL::ReferenceFrame referenceFrameFromNormal(const IRL::Normal a_normal);

// temporary particle method function declarations
std::vector<Eigen::Vector2d> ComputeParticlePositions(const int& N,
                                                      const Eigen::Vector2d& p,
                                                      const double& phi,
                                                      const double& theta,
                                                      const double& hp);

Eigen::Vector2d ComputeParticleForce(
    const Eigen::Vector2d& x,
    const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& segs,
    const double& eta);

std::vector<Eigen::Vector2d> InitializeParticlePositions(
    const std::pair<Eigen::Vector2d, Eigen::Vector2d>& target_endpoints,
    const double& hp, const int& N);

double ComputeParticleForceProjection(
    const int& N, const double& phi, const double& theta, const double& hp,
    const bool& iswrtPhi, const std::vector<Eigen::Vector2d> particle_forces);

double getParticleMethodCurvature(
    const std::pair<Eigen::Vector2d, Eigen::Vector2d>& target_endpoints,
    const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& endpoints,
    const int& N, const double& Hp, const double& h, const double& eta);

void recenterMoments(IRL::VolumeMoments* moments, const IRL::Pt& center);

void correctInterfaceBorders(Data<IRL::SeparatorVariant>* a_interface);

namespace details {
inline IRL::Paraboloid fromSphere(const IRL::Pt& a_center,
                                  const double a_radius,
                                  const IRL::Normal& a_normal) {
  const double curvature = 1.0 / a_radius;
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
  frame[1] = crossProduct(a_normal, frame[0]);
  frame[2] = a_normal;

  return IRL::Paraboloid(a_center + a_radius * a_normal, frame, 0.5 * curvature,
                         0.5 * curvature);
}

}  // namespace details

#endif  // EXAMPLES_VARIANT_ADVECTOR_RECONSTRUCTION_TYPES_H_
