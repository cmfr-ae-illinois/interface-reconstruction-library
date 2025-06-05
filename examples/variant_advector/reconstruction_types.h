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
#include "irl/moments/volume_moments.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include "examples/variant_advector/data.h"

void getReconstruction(
    const std::string& a_reconstruction_method,
    const Data<IRL::VolumeMoments>& a_liq_moments,
    const Data<IRL::VolumeMoments>& a_gas_moments,
    Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface,
    const double a_dt, const Data<double>& a_U, const Data<double>& a_V,
    const Data<double>& a_W, Data<IRL::SeparatorVariant>* a_interface);

struct PLIC {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface);
};

struct Jibben {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface,
      const bool a_plic_already_built = false);
};

struct MixedPLICJibben {
  static void getReconstruction(
      const Data<IRL::VolumeMoments>& a_liq_moments,
      const Data<IRL::VolumeMoments>& a_gas_moments, const double a_dt,
      const Data<double>& a_U, const Data<double>& a_V, const Data<double>& a_W,
      Data<IRL::SeparatorVariant>* a_interface,
      Data<IRL::LocalizedSeparatorVariantLink>* a_link_localized_interface);
};

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
