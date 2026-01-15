// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_IMPLICIT_SURFACE_CUTTING_CUT_IMPLICIT_SURFACE_H_
#define IRL_GENERIC_CUTTING_IMPLICIT_SURFACE_CUTTING_CUT_IMPLICIT_SURFACE_H_

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "examples/variant_advector/basic_mesh.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace IRL {

/// \brief Cutting cell by implicit surface
template <class SurfaceType, class ReturnType,
          class CellType = RectangularCuboid>
class ImplicitSurfaceCutter {
 public:
  using Vec3 = Eigen::Vector3d;
  using Mat3 = Eigen::Matrix3d;

  static_assert(std::is_same<typename SurfaceType::Scalar, double>::value,
                "SurfaceType must use double as ScalarType");

  static_assert(std::is_same<ReturnType, Volume>::value ||
                    std::is_same<ReturnType, VolumeMoments>::value ||
                    std::is_same<ReturnType, GeneralMoments3D<2>>::value,
                "ReturnType only supports Volume, "
                "VolumeMoments, or GeneralMoments3D<2>");

  static_assert(std::is_same<CellType, RectangularCuboid>::value,
                "ImplicitSurfaceCutter only supports "
                "RectangularCuboid as CellType");

  ImplicitSurfaceCutter(const SurfaceType& surface, const CellType& base_cell);

  ReturnType getMoments() const { return m_moments; }

 private:
  enum class CellStatus { Above, Below, Mixed };

  static constexpr size_t maxRefineLevel = SurfaceType::getMaxRefineLevel();

  struct Node {
    CellType cell;
    int level;
    CellStatus status;
    std::array<std::unique_ptr<Node>, 8> children;

    Node(const CellType& c, int lvl);
    bool isLeaf() const { return !children[0]; }
  };

  void divideCell(Node* node);

  void collectMixedLeaves(Node* node);

  Vec3 projectPointOnSurface(const Vec3& p) const;

  std::vector<Vec3> getSamplePoints(const CellType& cell,
                                    const bool& use_stencil) const;

  CellStatus getCellStatusFor(const CellType& cell) const;

  ReturnType computeMoments() const;

  Paraboloid makeParaboloidFor(const CellType& cell) const;

  SurfaceType m_surface;
  CellType m_cell;
  std::vector<Node*> mixed_leaves_;
  std::unique_ptr<Node> root_;
  ReturnType m_moments;
};

}  // namespace IRL

#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.tpp"

#endif  // IRL_GENERIC_CUTTING_IMPLICIT_SURFACE_CUTTING_CUT_IMPLICIT_SURFACE_H_