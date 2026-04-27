// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/interface_reconstruction_methods/c_reconstruction_metrics_interface.h"

#include <cassert>

extern "C" {

double c_reconstructionMetricWithJibben3D(
    const c_JibbenNeigh* a_jibben_neighborhood) {
  return IRL::reconstructionMetricWithJibben3D(
      *(a_jibben_neighborhood->obj_ptr));
}

double c_angularVarianceMetricWithJibben3D(
    const c_JibbenNeigh* a_jibben_neighborhood) {
  return IRL::angularVarianceMetricWithJibben3D(
      *(a_jibben_neighborhood->obj_ptr));
}

double c_volumeErrorSquaredWithJibben3D(
    const c_JibbenNeigh* a_jibben_neighborhood, const double* a_dx) {
  return IRL::volumeErrorSquaredWithJibben3D(*(a_jibben_neighborhood->obj_ptr),
                                             *a_dx);
}

}  // end extern C
