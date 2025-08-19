// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_JIBBEN_NEIGHBORHOOD_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_JIBBEN_NEIGHBORHOOD_H_

#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/polygons/polygon.h"

namespace IRL {

class JibbenNeighborhood {
  using iterator = typename std::vector<Polygon>::iterator;
  using const_iterator = typename std::vector<Polygon>::const_iterator;

 public:
  JibbenNeighborhood(void);

  void addMember(const Polygon& a_polygon, const double a_weight = 1.0);

  void emptyNeighborhood(void);

  void localize(void);

  void resize(const UnsignedIndex_t a_size);

  void reserve(const UnsignedIndex_t a_size);

  void setCenterOfStencil(const UnsignedIndex_t a_index);

  void setDelta(const double a_delta);

  UnsignedIndex_t getCenterOfStencilIndex(void) const;

  const Polygon& getCenterPolygon(void) const;

  const double& getWeight(const UnsignedIndex_t a_index) const;

  const Polygon& getPolygon(const UnsignedIndex_t a_index) const;

  const double& getDelta(void) const;

  const Pt& getDatum(void) const;

  const ReferenceFrame& getReferenceFrame(void) const;

  UnsignedIndex_t size(void) const;
  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator cend(void) const noexcept;

  ~JibbenNeighborhood(void) = default;

 private:
  void checkIndex(UnsignedIndex_t a_index) const;
  void checkCenterStencilSet(void) const;

  bool is_localized_m;
  double delta_m;
  Pt datum_m;
  ReferenceFrame frame_m;
  std::vector<double> weights_m;
  UnsignedIndex_t center_cell_index_m;
  std::vector<Polygon> polygons_m;
};

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_JIBBEN_NEIGHBORHOOD_H_
