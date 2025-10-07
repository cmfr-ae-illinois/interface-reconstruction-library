// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_RADIUS_SOLVER_CYLINDER_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_RADIUS_SOLVER_CYLINDER_TPP_

#include <algorithm>

#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_initializer.tpp"

namespace IRL {

template <class CellType>
ProgressiveRadiusSolverCylinder<CellType>::ProgressiveRadiusSolverCylinder(
    const CellType& a_cell, const double a_volume_fraction,
    const double a_volume_fraction_tolerance, const Cylinder& a_reconstruction)
    : target_volume_fraction_m(a_volume_fraction),
      volume_fraction_tolerance_m(a_volume_fraction_tolerance),
      reconstruction_m(a_reconstruction) {
  assert(a_volume_fraction >= 0.0);
  assert(a_volume_fraction <= 1.0);
  this->solveForRadius(a_cell);
}

template <class CellType>
void ProgressiveRadiusSolverCylinder<CellType>::solve(
    const CellType& a_cell, const double a_volume_fraction,
    const double a_volume_fraction_tolerance,
    const Cylinder& a_reconstruction) {
  assert(a_volume_fraction >= 0.0);
  assert(a_volume_fraction <= 1.0);
  target_volume_fraction_m = a_volume_fraction;
  volume_fraction_tolerance_m = a_volume_fraction_tolerance;
  reconstruction_m = a_reconstruction;
  this->solveForRadius(a_cell);
}

template <class CellType>
Cylinder ProgressiveRadiusSolverCylinder<CellType>::getCylinder(void) {
  reconstruction_m.setAlignedCylinder(AlignedCylinder({b_m, r_m}));
  return reconstruction_m;
}

template <class CellType>
void ProgressiveRadiusSolverCylinder<CellType>::solveForRadius(
    const CellType& a_cell) {
  auto copy_cell = CellType(a_cell);

  // Calculate volume of cell
  auto cell_volume = copy_cell.calculateVolume();
  double length_scale = std::cbrt(cell_volume);

  // Move cell to local frame of reference of the paraboloid
  const double b = reconstruction_m.getAlignedCylinder().b();
  const double radius = reconstruction_m.getAlignedCylinder().r();

  double interval_min = 0.0;
  double interval_max = radius;

  double vfrac_max =
      getVolumeMoments<Volume>(copy_cell, reconstruction_m) / cell_volume;
  UnsignedIndex_t iter = 0;
  while (iter < 40 && vfrac_max < 1.0) {
    interval_max *= 2.0;
    reconstruction_m.setAlignedCylinder(AlignedCylinder(
        {b < 0.0 ? b * interval_max / radius : b, interval_max}));
    vfrac_max =
        getVolumeMoments<Volume>(copy_cell, reconstruction_m) / cell_volume;
    iter++;
  }

  // Perform bisection since secant failed to find answer within tolerance.
  // Move cell back to its initial position
  std::array<double, 3> bounding_values{
      {interval_min, 0.5 * (interval_min + interval_max), interval_max}};
  for (UnsignedIndex_t iter = 0; iter < max_bisection_iter; ++iter) {
    bounding_values[1] = 0.5 * (bounding_values[0] + bounding_values[2]);
    reconstruction_m.setAlignedCylinder(AlignedCylinder(
        {b < 0.0 ? b * bounding_values[1] / radius : b, bounding_values[1]}));
    const double vfrac_cut =
        getVolumeMoments<Volume>(copy_cell, reconstruction_m) / cell_volume;
    if (vfrac_cut > target_volume_fraction_m + volume_fraction_tolerance_m) {
      bounding_values[2] = bounding_values[1];
    } else if (vfrac_cut <
               target_volume_fraction_m - volume_fraction_tolerance_m) {
      bounding_values[0] = bounding_values[1];
    } else {
      r_m = bounding_values[1];
      b_m = b < 0.0 ? b * bounding_values[1] / radius : b;
      return;
    }
  }

  std::cout << "WARNING: Bisection max iter reached inside volume fraction "
               "routine for cylinder"
            << std::endl;
  r_m = bounding_values[1];
  b_m = b < 0.0 ? b * bounding_values[1] / radius : b;

  return;
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_RADIUS_SOLVER_CYLINDER_TPP_
