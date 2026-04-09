// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_NEIGHBORHOOD_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_NEIGHBORHOOD_H_

#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/variant_reconstruction/separator_variant.h"

namespace IRL {

class PUNeighborhood {
 public:
  PUNeighborhood(void);

  void addMember(const SeparatorVariant& a_separator, const Pt& a_centroid,
                 const double a_weight = 1.0);

  void setMember(const UnsignedIndex_t a_index,
                 const SeparatorVariant& a_separator, const Pt& a_centroid,
                 const double a_weight = 1.0);

  void emptyNeighborhood(void);

  void resize(const UnsignedIndex_t a_size);

  void reserve(const UnsignedIndex_t a_size);

  void setCenterOfStencil(const UnsignedIndex_t a_index);

  // functions that do not have c and fortran interface
  UnsignedIndex_t size(void) const;

  const SeparatorVariant& getSeparator(const UnsignedIndex_t a_index) const;

  const double& getWeight(const UnsignedIndex_t a_index) const;

  const Pt& getCentroid(const UnsignedIndex_t a_index) const;

  const UnsignedIndex_t& getCenterOfStencil(void) const;

  ~PUNeighborhood(void) = default;

 private:
  std::vector<double> weights_m;
  std::vector<SeparatorVariant> separators_m;
  std::vector<Pt> centroids_m;
  UnsignedIndex_t center_cell_index_m;
};

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PU_NEIGHBORHOOD_H_
