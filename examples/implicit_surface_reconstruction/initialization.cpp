// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/implicit_surface_reconstruction/initialization.h"

// finding mixed cells
std::vector<std::tuple<int, int, int>> getCellStatus(Data<int>* cell_status) {
  std::vector<std::tuple<int, int, int>> mixed_cells_list;

  const BasicMesh& mesh = cell_status->getMesh();
  IRL::Sphere<double, 0> sphere(0., 0., 0., 0.15);

  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
        IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
        IRL::RectangularCuboid cell =
            IRL::RectangularCuboid::fromBoundingPts(x0, x1);
        IRL::ImplicitSurfaceCutter<IRL::Sphere<double, 0>, IRL::Volume> cutter(
            sphere, cell);
        (*cell_status)(i, j, k) = cutter.getBaseCellStatus();
        if ((*cell_status)(i, j, k) == 0)
          mixed_cells_list.emplace_back(i, j, k);
      }
    }
  }

  return mixed_cells_list;
}