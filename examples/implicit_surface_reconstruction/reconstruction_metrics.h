// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_RECONSTRUCTION_METRICS_H_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_RECONSTRUCTION_METRICS_H_

#include "examples/implicit_surface_reconstruction/basic_mesh.h"
#include "examples/implicit_surface_reconstruction/data.h"
#include "examples/implicit_surface_reconstruction/initialization.h"

// Indices for raw moment layout (ORDER>=2 assumed if you use M2)
enum : int {
  M0 = 0,
  M1x = 1,
  M1y = 2,
  M1z = 3,
  Mxx = 4,
  Mxy = 5,
  Mxz = 6,
  Myy = 7,
  Myz = 8,
  Mzz = 9
};

struct MomentDiffNorms {
  // volume
  double vol_M0_Linf, vol_M0_L2;
  double vol_M1_Linf, vol_M1_L2;
  double vol_M2_Linf, vol_M2_L2;
  // surface
  double surf_M0_Linf, surf_M0_L2;
  double surf_M1_Linf, surf_M1_L2;
  double surf_M2_Linf, surf_M2_L2;
};

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
MomentDiffNorms compute_moment_diff_norms(
    const Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                         IRL::GeneralSurfaceMoments3D<SM_ORDER>>>& A,
    const Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                         IRL::GeneralSurfaceMoments3D<SM_ORDER>>>& B);

#include "examples/implicit_surface_reconstruction/reconstruction_metrics.tpp"

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_RECONSTRUCTION_METRICS_H_
