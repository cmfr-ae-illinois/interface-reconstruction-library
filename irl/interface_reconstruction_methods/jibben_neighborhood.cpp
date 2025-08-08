// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/interface_reconstruction_methods/jibben_neighborhood.h"

namespace IRL {

JibbenNeighborhood::JibbenNeighborhood(void)
    : center_cell_index_m(static_cast<UnsignedIndex_t>(-1)),
      is_localized_m(false),
      delta_m(0.0) {}

void JibbenNeighborhood::addMember(const Polygon& a_polygon,
                                   const double a_weight) {
  polygons_m.push_back(a_polygon);
  weights_m.push_back(a_weight);
}

void JibbenNeighborhood::emptyNeighborhood(void) {
  polygons_m.resize(0);
  weights_m.resize(0);
  center_cell_index_m = static_cast<UnsignedIndex_t>(-1);
  is_localized_m = false;
}

void JibbenNeighborhood::localize(void) {
  if ((not is_localized_m) && polygons_m.size() > 0) {
    this->checkCenterStencilSet();
    this->checkIndex(center_cell_index_m);
    const Polygon& central_polygon = polygons_m[center_cell_index_m];
    datum_m = central_polygon.calculateCentroid();
    frame_m = ReferenceFrame::fromNormal(
        central_polygon.getPlaneOfExistence().normal());
    for (auto& polygon : polygons_m) {
      // Move vertices
      for (auto& pt : polygon) {
        const Pt tmp_pt = pt - datum_m;
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          pt[d] = frame_m[d] * tmp_pt;
        }
      }
      if (polygon.getNumberOfVertices() > 0) {
        // Reorient plane of polygon
        const auto& plane = polygon.getPlaneOfExistence();
        Normal normal = plane.normal();
        const Normal tmp_normal = plane.normal();
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          normal[d] = frame_m[d] * tmp_normal;
        }
        polygon.setPlaneOfExistence(Plane(normal, normal * polygon[0]));
      }
    }
    is_localized_m = true;
  }
}

void JibbenNeighborhood::setCenterOfStencil(const UnsignedIndex_t a_index) {
  this->checkIndex(a_index);
  center_cell_index_m = a_index;
}

void JibbenNeighborhood::setDelta(const double a_delta) { delta_m = a_delta; }

const Polygon& JibbenNeighborhood::getCenterPolygon(void) const {
  this->checkCenterStencilSet();
  this->checkIndex(center_cell_index_m);
  return polygons_m[center_cell_index_m];
}

const double& JibbenNeighborhood::getWeight(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return weights_m[a_index];
}

UnsignedIndex_t JibbenNeighborhood::getCenterOfStencilIndex(void) const {
  this->checkCenterStencilSet();
  this->checkIndex(center_cell_index_m);
  return center_cell_index_m;
}

const Polygon& JibbenNeighborhood::getPolygon(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return polygons_m[a_index];
}

const double& JibbenNeighborhood::getDelta(void) const { return delta_m; }

const Pt& JibbenNeighborhood::getDatum(void) const { return datum_m; }

const ReferenceFrame& JibbenNeighborhood::getReferenceFrame(void) const {
  return frame_m;
}

void JibbenNeighborhood::resize(const UnsignedIndex_t a_size) {
  polygons_m.resize(a_size);
  weights_m.resize(a_size);
}

void JibbenNeighborhood::reserve(const UnsignedIndex_t a_size) {
  polygons_m.reserve(a_size);
  weights_m.reserve(a_size);
}

UnsignedIndex_t JibbenNeighborhood::size(void) const {
  assert(polygons_m.size() == weights_m.size());
  return static_cast<UnsignedIndex_t>(polygons_m.size());
}

typename JibbenNeighborhood::iterator JibbenNeighborhood::begin(void) noexcept {
  return polygons_m.begin();
}

typename JibbenNeighborhood::const_iterator JibbenNeighborhood::begin(
    void) const noexcept {
  return this->cbegin();
}

typename JibbenNeighborhood::const_iterator JibbenNeighborhood::end(
    void) const noexcept {
  return this->cend();
}

typename JibbenNeighborhood::const_iterator JibbenNeighborhood::cbegin(
    void) const noexcept {
  return polygons_m.cbegin();
}

typename JibbenNeighborhood::iterator JibbenNeighborhood::end(void) noexcept {
  return polygons_m.end();
}

typename JibbenNeighborhood::const_iterator JibbenNeighborhood::cend(
    void) const noexcept {
  return polygons_m.cend();
}

void JibbenNeighborhood::checkIndex(UnsignedIndex_t a_index) const {
  assert(a_index < polygons_m.size());
}

void JibbenNeighborhood::checkCenterStencilSet(void) const {
  assert(center_cell_index_m != static_cast<UnsignedIndex_t>(-1));
}

}  // namespace IRL
