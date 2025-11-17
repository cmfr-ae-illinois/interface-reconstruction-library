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

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/moments/general_surface_moments.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

namespace IRL {

/// \brief Cutting cell by implicit surface
template <class ImplicitSurfaceType, class ReturnType,
          class CellType = RectangularCuboid>
class ImplicitSurfaceCutter {
 public:
  using Vec3 = Eigen::Vector3d;
  using Mat3 = Eigen::Matrix3d;

  static_assert(
      std::is_same<typename ImplicitSurfaceType::Scalar, double>::value,
      "ImplicitSurfaceType must use double as ScalarType");

  static_assert(std::is_same<ReturnType, Volume>::value ||
                    std::is_same<ReturnType, VolumeMoments>::value ||
                    std::is_same<ReturnType, GeneralMoments3D<2>>::value,
                "ReturnType only supports Volume, "
                "VolumeMoments, or GeneralMoments3D<2>");

  static_assert(std::is_same<CellType, RectangularCuboid>::value,
                "ImplicitSurfaceCutter only supports "
                "RectangularCuboid as CellType");

  ImplicitSurfaceCutter(const ImplicitSurfaceType& surface,
                        const CellType& base_cell);

  ReturnType computeVolumeMoments() const;

  template <std::size_t ORDER>
  GeneralSurfaceMoments3D<ORDER> computeSurfaceMoments(
      const bool useAdaptive = true,
      const Eigen::Integrator<double, 2>::QuadratureRule quadratureRule =
          Eigen::Integrator<double, 2>::GaussKronrod15,
      const int npts = 50) const;

  int getBaseCellStatus() const;

 private:
  enum class CellStatus { Above, Below, Mixed };

  static constexpr size_t maxRefineLevel =
      ImplicitSurfaceType::getMaxRefineLevel();

  struct Node {
    CellType cell;
    int level;
    CellStatus status;
    std::array<std::unique_ptr<Node>, 8> children;

    Node(const CellType& c, int lvl);
    bool isLeaf() const { return !children[0]; }
  };

  void divideCell(Node* node);

  void collectLeaves(Node* node);

  Vec3 projectPointOnSurface(const Vec3& p) const;

  std::vector<Vec3> getSamplePoints(const CellType& cell,
                                    const bool& use_stencil) const;

  CellStatus getCellStatusFor(const CellType& cell) const;

  Paraboloid makeParaboloidFor(const CellType& cell) const;

  ImplicitSurfaceType m_surface;
  CellType m_cell;
  std::vector<Node*> mixed_leaves_;
  std::vector<Node*> below_leaves_;
  std::vector<Node*> above_leaves_;
  std::vector<Paraboloid> m_paraboloids;
  std::unique_ptr<Node> root_;
};

}  // namespace IRL

#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.tpp"

#endif  // IRL_GENERIC_CUTTING_IMPLICIT_SURFACE_CUTTING_CUT_IMPLICIT_SURFACE_H_