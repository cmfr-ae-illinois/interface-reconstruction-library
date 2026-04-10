// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Andrew Cahaly <ajc428@cornell.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PLICNET_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PLICNET_H_

namespace IRL {

class PLICNet {
 public:
  /// \brief Default constructor.
  PLICNet(void);

  /// \brief Clear stencil information.
  void clear(void);

  /// \brief Set member information given (ijk) position in stencil
  /// in the discrete space {-1,0,1}^3
  void setMember(const Pt& a_lower_cell_pt, const Pt& a_upper_cell_pt,
                 const double a_volume_fraction, const Pt& a_centroid_0,
                 const Pt& a_centroid_1, const int a_i, const int a_j,
                 const int a_k);

  /// \brief Reflects moments, executes neural net, and returns planar separator
  // in original orientation
  PlanarSeparator getPlanarSeparator(void);

  /// \brief Default destructor.
  ~PLICNet(void) = default;

 private:
  void calculateNormal(void);
  void flip(void);
  void computeStencilMoments(void);
  void reflectMomentsX(void);
  void reflectMomentsY(void);
  void reflectMomentsZ(void);
  void reflectMomentsXY(void);
  void reflectMomentsYZ(void);
  void reflectMomentsXZ(void);
  void reflectMoments(void);

  /// \brief Storage of the stencil information
  std::array<double, 189> data_m;
  bool flipped_m, central_cell_defined_m;
  RectangularCuboid cell_m;
  UnsignedIndex_t dir_m, dir2_m;
  double stencil_m0_m;
  Pt stencil_m1_m;
  Normal normal_m;
};

}  // namespace IRL

#include "irl/interface_reconstruction_methods/plicnet.tpp"

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PLICNET_H_
