// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_IMPLICIT_SURFACE_CUTTING_CUT_IMPLICIT_SURFACE_TPP_
#define IRL_GENERIC_CUTTING_IMPLICIT_SURFACE_CUTTING_CUT_IMPLICIT_SURFACE_TPP_

#include "irl/generic_cutting/implicit_surface_cutting/cut_implicit_surface.h"

namespace IRL {

// constructor definition
template <class ImplicitSurfaceType, class ReturnType, class CellType>
ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType, CellType>::
    ImplicitSurfaceCutter(const ImplicitSurfaceType& surface,
                          const CellType& base_cell)
    : m_surface(surface), m_cell(base_cell) {
  root_ = std::make_unique<Node>(base_cell, 0);
  divideCell(root_.get());
  collectLeaves(root_.get());
}

// constructor definition for node of octree
template <class ImplicitSurfaceType, class ReturnType, class CellType>
ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType, CellType>::Node::Node(
    const CellType& c, int lvl)
    : cell(c), level(lvl), status(CellStatus::Mixed) {}

// projecting point onto implicit surface
template <class ImplicitSurfaceType, class ReturnType, class CellType>
typename ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType, CellType>::Vec3
ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                      CellType>::projectPointOnSurface(const Vec3& p) const {
  const int max_iter = 200;
  const double tol = 1e-10;
  Vec3 x_proj = p;
  for (int i = 0; i < max_iter; i++) {
    double f = m_surface.F(x_proj[0], x_proj[1], x_proj[2]);
    Vec3 g = m_surface.gradF(x_proj[0], x_proj[1], x_proj[2]);
    double g_norm2 = g.squaredNorm();
    if (g_norm2 < 1e-14) break;
    x_proj -= (f / g_norm2) * g;
    if (std::abs(f) < tol) break;
    if (i == (max_iter - 1)) {
      std::cout << "Max iterations reached. Projection incomplete. "
                << "f = " << std::abs(f) << std::endl;
    }
  }
  return x_proj;
}

// sample pooints on/around surface for getting cell status
template <class ImplicitSurfaceType, class ReturnType, class CellType>
std::vector<typename ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                                           CellType>::Vec3>
ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                      CellType>::getSamplePoints(const CellType& cell,
                                                 const bool& use_stencil)
    const {
  Pt x0 = cell.getLowerLimits();
  Pt x1 = cell.getUpperLimits();
  double dx = std::abs((x1 - x0)[0]);
  double dy = std::abs((x1 - x0)[1]);
  double dz = std::abs((x1 - x0)[2]);

  double sx = x0[0] - (use_stencil ? dx : 0.0);
  double sy = x0[1] - (use_stencil ? dy : 0.0);
  double sz = x0[2] - (use_stencil ? dz : 0.0);
  double hx = (use_stencil ? 3.0 * dx : dx);
  double hy = (use_stencil ? 3.0 * dy : dy);
  double hz = (use_stencil ? 3.0 * dz : dz);

  std::vector<Vec3> pts;
  // corners
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        pts.emplace_back(sx + static_cast<double>(i) * hx,
                         sy + static_cast<double>(j) * hy,
                         sz + static_cast<double>(k) * hz);
      }
    }
  }
  // face centers
  pts.emplace_back(sx + hx / 2.0, sy + hy / 2.0, sz);
  pts.emplace_back(sx + hx / 2.0, sy, sz + hz / 2.0);
  pts.emplace_back(sx, sy + hy / 2.0, sz + hz / 2.0);
  pts.emplace_back(sx + hx, sy + hy / 2.0, sz + hz / 2.0);
  pts.emplace_back(sx + hx / 2.0, sy + hy, sz + hz / 2.0);
  pts.emplace_back(sx + hx / 2.0, sy + hy / 2.0, sz + hz);

  return pts;
}

// cell status
template <class ImplicitSurfaceType, class ReturnType, class CellType>
typename ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                               CellType>::CellStatus
ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                      CellType>::getCellStatusFor(const CellType& cell) const {
  // check cell points
  std::vector<Vec3> pts = getSamplePoints(cell, false);
  bool all_pos = false, all_neg = false;
  for (const auto& pt : pts) {
    double Fval = m_surface.F(pt[0], pt[1], pt[2]);
    if (Fval > 0.0)
      all_pos = true;
    else if (Fval < 0.0)
      all_neg = true;
    if (all_pos && all_neg) return CellStatus::Mixed;
  }

  // check stencil poitns
  pts = getSamplePoints(cell, true);
  for (const auto& pt : pts) {
    double Fval = m_surface.F(pt[0], pt[1], pt[2]);
    if ((Fval > 0.0 && !all_pos) || (Fval < 0.0 && !all_neg)) {
      Paraboloid paraboloid = makeParaboloidFor(cell);
      double liq_vf =
          getVolumeMoments<Volume>(cell, paraboloid) / cell.calculateVolume();
      bool not_mixed = true;
      if (liq_vf >= IRL::global_constants::VF_LOW &&
          liq_vf <= IRL::global_constants::VF_HIGH)
        not_mixed = false;
      return not_mixed ? (all_pos ? CellStatus::Above : CellStatus::Below)
                       : CellStatus::Mixed;
    }
  }
  return all_pos ? CellStatus::Above : CellStatus::Below;
}

// subdivide node if mixed and below max level
template <class ImplicitSurfaceType, class ReturnType, class CellType>
void ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                           CellType>::divideCell(Node* node) {
  node->status = getCellStatusFor(node->cell);
  if (node->status == CellStatus::Mixed &&
      node->level < static_cast<int>(maxRefineLevel)) {
    Pt n_x0 = node->cell.getLowerLimits();
    Pt n_x1 = node->cell.getUpperLimits();
    Pt n_xm = 0.5 * (n_x0 + n_x1);
    for (int i = 0; i < 8; i++) {
      Pt child_x0, child_x1;
      child_x0[0] = (i & 1) ? n_xm[0] : n_x0[0];
      child_x0[1] = (i & 2) ? n_xm[1] : n_x0[1];
      child_x0[2] = (i & 4) ? n_xm[2] : n_x0[2];
      child_x1[0] = (i & 1) ? n_x1[0] : n_xm[0];
      child_x1[1] = (i & 2) ? n_x1[1] : n_xm[1];
      child_x1[2] = (i & 4) ? n_x1[2] : n_xm[2];

      CellType subcell = CellType::fromBoundingPts(child_x0, child_x1);
      node->children[i] = std::make_unique<Node>(subcell, node->level + 1);
      divideCell(node->children[i].get());
    }
  }
}

// pointers to above, below and mixed leaves
template <class ImplicitSurfaceType, class ReturnType, class CellType>
void ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                           CellType>::collectLeaves(Node* node) {
  if (node->isLeaf()) {
    if (node->status == CellStatus::Mixed) {
      mixed_leaves_.push_back(node);
      m_paraboloids.push_back(makeParaboloidFor(node->cell));
    } else if (node->status == CellStatus::Below) {
      below_leaves_.push_back(node);
    } else {
      above_leaves_.push_back(node);
    }
  } else {
    for (auto& c : node->children) collectLeaves(c.get());
  }
}

// project and construct paraboloid
template <class ImplicitSurfaceType, class ReturnType, class CellType>
Paraboloid
ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                      CellType>::makeParaboloidFor(const CellType& cell) const {
  Pt x0 = cell.getLowerLimits();
  Pt x1 = cell.getUpperLimits();
  Pt xm = (x0 + x1) * 0.5;
  Vec3 center(xm[0], xm[1], xm[2]);
  Vec3 x_proj = projectPointOnSurface(center);
  Pt x_proj_pt(x_proj(0), x_proj(1), x_proj(2));
  Vec3 g = m_surface.gradF(x_proj(0), x_proj(1), x_proj(2));
  Mat3 H = m_surface.hessF(x_proj(0), x_proj(1), x_proj(2));

  return Paraboloid::fromDerivatives(x_proj_pt, g, H);
}

// accumulate volume moments
template <class ImplicitSurfaceType, class ReturnType, class CellType>
ReturnType ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                                 CellType>::computeVolumeMoments() const {
  ReturnType sum{};  // accumulator
  ReturnType c{};    // compensation

  auto kahan_add = [&](const ReturnType& x) {
    const ReturnType y = x - c;
    const ReturnType t = sum + y;
    c = (t - sum) - y;
    sum = t;
  };

  // mixed contribution
  for (std::size_t i = 0; i < mixed_leaves_.size(); i++) {
    const ReturnType moment_contribution =
        getVolumeMoments<ReturnType>(mixed_leaves_[i]->cell, m_paraboloids[i]);
    kahan_add(moment_contribution);
  }
  // below contribution
  for (std::size_t i = 0; i < below_leaves_.size(); i++) {
    const ReturnType moment_contribution =
        getVolumeMoments<ReturnType>(below_leaves_[i]->cell);
    kahan_add(moment_contribution);
  }

  return sum;
}

// accumulate surface moments
template <class ImplicitSurfaceType, class ReturnType, class CellType>
template <std::size_t ORDER>
GeneralSurfaceMoments3D<ORDER>
ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType, CellType>::
    computeSurfaceMoments(
        const bool useAdaptive,
        const Eigen::Integrator<double, 2>::QuadratureRule quadratureRule,
        const int npts) const {
  using VolumeMomentsAndSurface =
      AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

  GeneralSurfaceMoments3D<ORDER> sum{};  // accumulator
  GeneralSurfaceMoments3D<ORDER> c{};    // compensation

  auto kahan_add = [&](const GeneralSurfaceMoments3D<ORDER>& x) {
    const auto y = x - c;
    const auto t = sum + y;
    c = (t - sum) - y;
    sum = t;
  };

  for (std::size_t i = 0; i < mixed_leaves_.size(); i++) {
    auto surface = getVolumeMoments<VolumeMomentsAndSurface>(
                       mixed_leaves_[i]->cell, m_paraboloids[i])
                       .getSurface();
    const auto term = surface.template getSurfaceMoments<ORDER>(
        useAdaptive, quadratureRule, npts);
    kahan_add(term);
  }

  return sum;
}

// -1: below, 0: mixed, 1: above
template <class ImplicitSurfaceType, class ReturnType, class CellType>
int ImplicitSurfaceCutter<ImplicitSurfaceType, ReturnType,
                          CellType>::getBaseCellStatus() const {
  // auto cell = root_.get();
  root_->status = getCellStatusFor(root_->cell);
  if (root_->status == CellStatus::Above) {
    return 1;
  } else if (root_->status == CellStatus::Below) {
    return -1;
  } else {
    return 0;
  }
}

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_IMPLICIT_SURFACE_CUTTING_CUT_IMPLICIT_SURFACE_TPP_