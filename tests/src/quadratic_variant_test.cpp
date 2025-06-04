// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// #define DEBUG_CYL_IRL
// #define NUDGE_REGION

#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/rotations.h"

#include <sys/time.h>
#include <cmath>
#include <random>
#include <variant>

#include "gtest/gtest.h"

#include "irl/generic_cutting/generic_cutting.h"

#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection.h"
#include "irl/moments/general_moments.h"

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/interface_reconstruction_methods/progressive_radius_solver_cylinder.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

#include "irl/variant_reconstruction/separator_variant.h"

namespace {
using namespace IRL;

TEST(QuadraticVariant, GetIf) {
  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto cell =
      RectangularCuboid::fromBoundingPts(Pt(-1, -1, -1), Pt(1, 1, 1));

  SeparatorVariant interface;
  interface = Paraboloid(datum, frame, 0.5, 0.5);

  Volume volume;
  if (const auto plane = std::get_if<PlanarSeparator>(&interface))
    volume = getVolumeMoments<Volume>(cell, *plane);
  else if (const auto paraboloid = std::get_if<Paraboloid>(&interface))
    volume = getVolumeMoments<Volume>(cell, *paraboloid);
  else if (const auto cylinder = std::get_if<Cylinder>(&interface))
    volume = getVolumeMoments<Volume>(cell, *cylinder);
  EXPECT_NEAR(volume, 4.0 - 4.0 / 3.0, 1.0e-14);

  interface = Cylinder(datum, frame, 1.0, 0.25);
  if (const auto plane = std::get_if<PlanarSeparator>(&interface))
    volume = getVolumeMoments<Volume>(cell, *plane);
  else if (const auto paraboloid = std::get_if<Paraboloid>(&interface))
    volume = getVolumeMoments<Volume>(cell, *paraboloid);
  else if (const auto cylinder = std::get_if<Cylinder>(&interface))
    volume = getVolumeMoments<Volume>(cell, *cylinder);
  EXPECT_NEAR(volume, M_PI * 0.25 * 2.0, 1.0e-14);
}

TEST(QuadraticVariant, getVolumeMoments) {
  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto cell =
      RectangularCuboid::fromBoundingPts(Pt(-1, -1, -1), Pt(1, 1, 1));

  SeparatorVariant interface;
  interface = Paraboloid(datum, frame, 0.5, 0.5);

  auto volume = getVolumeMoments<Volume>(cell, interface);
  EXPECT_NEAR(volume, 4.0 - 4.0 / 3.0, 1.0e-14);

  interface = Cylinder(datum, frame, 1.0, 0.25);
  volume = getVolumeMoments<Volume>(cell, interface);
  EXPECT_NEAR(volume, M_PI * 0.25 * 2.0, 1.0e-14);
}

TEST(QuadraticVariant, Localizer) {
  const auto ijk_frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto jki_frame =
      ReferenceFrame(Normal(0, 1, 0), Normal(0, 0, 1), Normal(1, 0, 0));

  // Create two adjacent cells
  std::array<RectangularCuboid, 3> cells;
  cells[0] = RectangularCuboid::fromBoundingPts(Pt(0, 0, 0), Pt(1, 1, 1));
  cells[1] = RectangularCuboid::fromBoundingPts(Pt(0, 1, 0), Pt(1, 2, 1));
  cells[2] = RectangularCuboid::fromBoundingPts(Pt(0, 2, 0), Pt(1, 3, 1));

  // Get localizers corresponding to cells
  std::array<PlanarLocalizer, 3> cell_localizers;
  for (int i = 0; i < 3; i++) cell_localizers[i] = cells[i].getLocalizer();

  // Create 2 interfaces: 1 paraboloid and 1 cylinder
  std::array<SeparatorVariant, 3> interfaces;
  interfaces[0] = Paraboloid(Pt(0.5, 0.5, 0.5), ijk_frame, 0.25, 0.25);
  interfaces[1] = Cylinder(Pt(0.0, 1.0, 0.0), ijk_frame, 1.0, 0.25);
  interfaces[2] = Cylinder(Pt(0.0, 2.0, 0.0), jki_frame, 1.0, 0.25);

  auto volume0 = getVolumeMoments<Volume>(cells[0], interfaces[0]);
  EXPECT_NEAR(volume0, 0.5 - 1.0 / 24.0, 1.0e-14);
  auto volume1 = getVolumeMoments<Volume>(cells[1], interfaces[1]);
  EXPECT_NEAR(volume1, M_PI / 16.0, 1.0e-14);
  auto volume2 = getVolumeMoments<Volume>(cells[2], interfaces[2]);
  EXPECT_NEAR(volume2, M_PI / 16.0, 1.0e-14);

  // Link interfaces to cell localizers
  std::array<LocalizedSeparatorVariantLink, 3> linked_localized_interfaces;
  for (int i = 0; i < 3; i++)
    linked_localized_interfaces[i] =
        LocalizedSeparatorVariantLink(&(cell_localizers[i]), &(interfaces[i]));

  // Connect linked localized interfaces
  for (int i = 0; i < 3; i++) {
    linked_localized_interfaces[i].setId(i);
    for (int j = 0; j < 6; j++)
      linked_localized_interfaces[i].setEdgeConnectivity(j, nullptr);
  }
  linked_localized_interfaces[0].setEdgeConnectivity(
      3, &(linked_localized_interfaces[1]));
  linked_localized_interfaces[1].setEdgeConnectivity(
      2, &(linked_localized_interfaces[0]));
  linked_localized_interfaces[1].setEdgeConnectivity(
      3, &(linked_localized_interfaces[2]));
  linked_localized_interfaces[2].setEdgeConnectivity(
      2, &(linked_localized_interfaces[1]));

  // Create large cell that overlaps all cells
  const auto large_cell =
      RectangularCuboid::fromBoundingPts(Pt(0, 0.5, 0), Pt(1, 2.5, 1));

  // Compute volume of intersection between large cell and both localized
  // interfaces
  auto volumefrom0 =
      getVolumeMoments<Volume>(large_cell, linked_localized_interfaces[0]);
  EXPECT_NEAR(volumefrom0, (0.5 - 1.0 / 24.0) / 2.0 + 1.5 * M_PI / 16.0,
              1.0e-14);
  auto volumefrom1 =
      getVolumeMoments<Volume>(large_cell, linked_localized_interfaces[1]);
  EXPECT_NEAR(volumefrom1, (0.5 - 1.0 / 24.0) / 2.0 + 1.5 * M_PI / 16.0,
              1.0e-14);
  auto volumefrom2 =
      getVolumeMoments<Volume>(large_cell, linked_localized_interfaces[2]);
  EXPECT_NEAR(volumefrom2, (0.5 - 1.0 / 24.0) / 2.0 + 1.5 * M_PI / 16.0,
              1.0e-14);
}
}  // namespace
