// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_CYLINDER_RECONSTRUCTIONS_CYLINDER_H_
#define IRL_CYLINDER_RECONSTRUCTIONS_CYLINDER_H_

#include <math.h>
// extern "C" {
// #include <quadmath.h>
// }
#include <cassert>
#include <ostream>

#include "irl/data_structures/stack_vector.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/graphs/un_directed_graph_node.h"
#include "irl/quadratic_reconstruction/quadratic_helper.h"
#include "irl/cylinder_reconstruction/aligned_cylinder.h"
// #include "irl/paraboloid_reconstruction/gradient_paraboloid.h"
#include "irl/planar_reconstruction/joined_reconstructions.h"
#include "irl/planar_reconstruction/planar_localizer.h"
#include "irl/planar_reconstruction/reconstruction_link.h"

namespace IRL {

// Cylinder represented by a datum,
// reference frame, and an AlignedCylinder defined as
// z^2 + b*y^2 = r
template <class ScalarType>
class CylinderBase {
 public:
  using value_type = ScalarType;
  CylinderBase(void);

  /// a_reference_frame should have the normal vector in a_reference_frame[2],
  /// and the tangent vectors corresponding to the axe of the cylinder and a_coef_b in
  /// element 0 and 1, respectively. (a_coef_r correspond to the radius squared, like for AlignedCylinderBase)
  CylinderBase(const PtBase<ScalarType>& a_datum,
                 const ReferenceFrameBase<ScalarType>& a_reference_frame,
                 const ScalarType a_coef_b, const ScalarType a_coef_r);

  static CylinderBase createAlwaysAbove(void);

  static CylinderBase createAlwaysBelow(void);

  void setDatum(const PtBase<ScalarType>& a_datum);
  void setReferenceFrame(
      const ReferenceFrameBase<ScalarType>& a_reference_frame);
  void setAlignedCylinder(
      const AlignedCylinderBase<ScalarType>& a_aligned_cylinder);

  const PtBase<ScalarType>& getDatum(void) const;
  const ReferenceFrameBase<ScalarType>& getReferenceFrame(void) const;
  const AlignedCylinderBase<ScalarType>& getAlignedCylinder(void) const;

  /// Indicates that the intersection should actually be performed.
  void markAsRealReconstruction(void);

  /// Marks cylinder as being above any polyhedron (so any polyhedron will
  /// be unclipped).
  void markAsAlwaysAbove(void);

  /// Marks cylinder as being below any polyhedron (so any polyhedron will
  /// be clipped).
  void markAsAlwaysBelow(void);

  /// Whether the cylinder has been set to be above any polyhedron.
  bool isAlwaysAbove(void) const;

  /// Whether the cylinder has been set to be below any polyhedron.
  bool isAlwaysBelow(void) const;

  /// Cylinder cannot be a flipped reconstruction. Add this for ease of use
  /// with other routines that usually take planar reconstructions.
  static constexpr bool isFlipped(void) { return false; }

  /// \brief Since localizers are always convex, never flip.
  static constexpr ScalarType flip(void) { return static_cast<ScalarType>(1); }

  /// \brief Return if cutting for gas phase is needed.
  static constexpr bool isNotFlipped(void) { return true; }

  void serialize(ByteBuffer* a_buffer) const;
  void unpackSerialized(ByteBuffer* a_buffer);

  ~CylinderBase(void) = default;

 private:
  PtBase<ScalarType> datum_m;
  ReferenceFrameBase<ScalarType> frame_m;
  AlignedCylinderBase<ScalarType> cylinder_m;
  std::array<bool, 2> place_infinite_shortcut_m;
};

using Cylinder = CylinderBase<double>;

template <class ScalarType>
using LocalizedCylinder =
    JoinedReconstructions<PlanarLocalizer, CylinderBase<ScalarType>>;
// using LocalizedParaboloid = LocalizedParaboloidBase<double>;

template <class ScalarType>
using LocalizedCylinderLink =
    ReconstructionLink<LocalizedCylinder<ScalarType>, UnDirectedGraphNode>;
// using LocalizedParaboloidLink = LocalizedParaboloidLinkBase<double>;

template <class ScalarType>
inline PtBase<ScalarType> conicCenter(
    const PlaneBase<ScalarType>& a_plane,
    const AlignedCylinderBase<ScalarType>& a_cylinder);

template <class ScalarType>
inline NormalBase<ScalarType> getCylinderSurfaceNormal(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const PtBase<ScalarType>& a_pt);

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongLineOntoCylinder(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt);

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongHalfLineOntoCylinder(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt);

template <class ScalarType>
inline std::ostream& operator<<(std::ostream& out,
                                const CylinderBase<ScalarType>& a_cylinder);

}  // namespace IRL

#include "irl/cylinder_reconstruction/cylinder.tpp"

#endif  // IRL_CYLINDER_RECONSTRUCTIONS_CYLINDER_H_
