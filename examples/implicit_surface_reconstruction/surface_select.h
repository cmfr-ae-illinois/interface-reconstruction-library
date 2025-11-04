// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_IMPLICIT_SURFACE_SURFACE_SELECT_H_
#define EXAMPLES_IMPLICIT_SURFACE_SURFACE_SELECT_H_

#include <iostream>
#include <string>
#include <variant>

#include "examples/implicit_surface_reconstruction/basic_mesh.h"
#include "irl/geometry/implicit_surfaces/general_implicit_surface.h"
#include "irl/geometry/implicit_surfaces/implicit_surfaces.h"

// refine levels for initializing moments (for Nx = 256)
constexpr std::size_t SPHERE_MAX_REFINE = 5;
constexpr std::size_t ELLIPSOID_MAX_REFINE = 5;
constexpr std::size_t GENUS_MAX_REFINE = 5;
constexpr std::size_t ORTHOCIRCLE_MAX_REFINE = 5;

using SphereVariant = IRL::Sphere<double, SPHERE_MAX_REFINE>;
using EllipsoidVariant = IRL::Ellipsoid<double, ELLIPSOID_MAX_REFINE>;
using GenusVariant = IRL::Genus<double, GENUS_MAX_REFINE>;
using OrthocircleVariant = IRL::Orthocircle<double, ORTHOCIRCLE_MAX_REFINE>;

using SurfaceVariant = std::variant<SphereVariant, EllipsoidVariant,
                                    GenusVariant, OrthocircleVariant>;

inline SurfaceVariant makeSurface(const std::string& name, BasicMesh& mesh) {
  if (name == "sphere") {
    mesh.setCellBoundaries(IRL::Pt(-0.18, -0.18, -0.18),
                           IRL::Pt(0.18, 0.18, 0.18));
    return SphereVariant(0., 0., 0., 0.15);
  } else if (name == "ellipsoid") {
    mesh.setCellBoundaries(IRL::Pt(-0.32, -0.32, -0.32),
                           IRL::Pt(0.32, 0.32, 0.32));
    return EllipsoidVariant(0., 0., 0., 0.30, 0.15, 0.10);
  } else if (name == "genus") {
    mesh.setCellBoundaries(IRL::Pt(-2, -2, -2), IRL::Pt(2, 2, 2));
    return GenusVariant{};
  } else if (name == "orthocircle") {
    mesh.setCellBoundaries(IRL::Pt(-1.5, -1.5, -1.5), IRL::Pt(1.5, 1.5, 1.5));
    return OrthocircleVariant{};
  }
  throw std::runtime_error("Unknown surface: " + name +
                           " (expected: sphere|ellipsoid|genus|orthocircle)");
}

inline std::string binary_filename(const std::string& outdir,
                                   const std::string& shape, int Nx) {
  return outdir + "/" + shape + "_moments_Nx" + std::to_string(Nx) + ".bin";
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_SURFACE_SELECT_H_
