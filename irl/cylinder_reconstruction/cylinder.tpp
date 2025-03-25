// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_CYLINDER_RECONSTRUCTIONS_CYLINDER_TPP_
#define IRL_CYLINDER_RECONSTRUCTIONS_CYLINDER_TPP_

#include "irl/cylinder_reconstruction/cylinder.h"

namespace IRL {

template <class ScalarType>
inline CylinderBase<ScalarType>::CylinderBase(void)
    : datum_m(),
      frame_m(),
      cylinder_m(),
      place_infinite_shortcut_m({false, false}) {}

template <class ScalarType>
inline CylinderBase<ScalarType>::CylinderBase(
    const PtBase<ScalarType>& a_datum,
    const ReferenceFrameBase<ScalarType>& a_reference_frame,
    const ScalarType a_coef_b, const ScalarType a_coef_r)
    : datum_m(a_datum),
      frame_m(a_reference_frame),
      cylinder_m(std::array<ScalarType, 2>({a_coef_b, a_coef_r})),
      place_infinite_shortcut_m({false, false}) {
  // assert(frame_m.isOrthonormalBasis());
}

template <class ScalarType>
inline CylinderBase<ScalarType> CylinderBase<ScalarType>::createAlwaysAbove(
    void) {
  CylinderBase<ScalarType> par;
  par.markAsAlwaysAbove();
  return par;
}

template <class ScalarType>
inline CylinderBase<ScalarType> CylinderBase<ScalarType>::createAlwaysBelow(
    void) {
  CylinderBase<ScalarType> par;
  par.markAsAlwaysBelow();
  return par;
}

template <class ScalarType>
inline void CylinderBase<ScalarType>::setDatum(
    const PtBase<ScalarType>& a_datum) {
  datum_m = a_datum;
}

template <class ScalarType>
inline void CylinderBase<ScalarType>::setReferenceFrame(
    const ReferenceFrameBase<ScalarType>& a_reference_frame) {
  assert(a_reference_frame.isOrthonormalBasis());
  frame_m = a_reference_frame;
}

template <class ScalarType>
inline void CylinderBase<ScalarType>::setAlignedCylinder(
    const AlignedCylinderBase<ScalarType>& a_aligned_cylinder) {
  cylinder_m = a_aligned_cylinder;
}

template <class ScalarType>
inline const PtBase<ScalarType>& CylinderBase<ScalarType>::getDatum(
    void) const {
  return datum_m;
}

template <class ScalarType>
inline const ReferenceFrameBase<ScalarType>&
CylinderBase<ScalarType>::getReferenceFrame(void) const {
  return frame_m;
}

template <class ScalarType>
inline const AlignedCylinderBase<ScalarType>&
CylinderBase<ScalarType>::getAlignedCylinder(void) const {
  return cylinder_m;
}

template <class ScalarType>
inline void CylinderBase<ScalarType>::markAsRealReconstruction(void) {
  place_infinite_shortcut_m[0] = false;
  place_infinite_shortcut_m[1] = false;
}

template <class ScalarType>
inline void CylinderBase<ScalarType>::markAsAlwaysAbove(void) {
  place_infinite_shortcut_m[0] = true;
  place_infinite_shortcut_m[1] = false;
}

template <class ScalarType>
inline void CylinderBase<ScalarType>::markAsAlwaysBelow(void) {
  place_infinite_shortcut_m[0] = false;
  place_infinite_shortcut_m[1] = true;
}

template <class ScalarType>
inline bool CylinderBase<ScalarType>::isAlwaysAbove(void) const {
  return place_infinite_shortcut_m[0];
}

template <class ScalarType>
inline bool CylinderBase<ScalarType>::isAlwaysBelow(void) const {
  return place_infinite_shortcut_m[1];
}

template <class ScalarType>
inline void CylinderBase<ScalarType>::serialize(ByteBuffer* a_buffer) const {
  datum_m.serialize(a_buffer);
  for (std::size_t d = 0; d < 3; ++d) {
    frame_m[d].serialize(a_buffer);
  }
  cylinder_m.serialize(a_buffer);
  const UnsignedIndex_t bool_to_int =
      (place_infinite_shortcut_m[0] ? 1 : 0) +
      2 * (place_infinite_shortcut_m[1] ? 1 : 0);
  a_buffer->pack(&bool_to_int, 1);
}

template <class ScalarType>
inline void CylinderBase<ScalarType>::unpackSerialized(ByteBuffer* a_buffer) {
  datum_m.unpackSerialized(a_buffer);
  NormalBase<ScalarType> normal;
  for (std::size_t d = 0; d < 3; ++d) {
    frame_m[d].unpackSerialized(a_buffer);
  }
  cylinder_m.unpackSerialized(a_buffer);
  UnsignedIndex_t int_to_bool = 0;
  a_buffer->unpack(&int_to_bool, 1);
  place_infinite_shortcut_m[0] = (int_to_bool % 2 == 1) ? true : false;
  place_infinite_shortcut_m[1] = (int_to_bool / 2 == 1) ? true : false;
}

template <class ScalarType>
inline PtBase<ScalarType> conicCenter(
    const PlaneBase<ScalarType>& a_plane,
    const AlignedCylinderBase<ScalarType>& a_cylinder) {
  const auto& face_normal = a_plane.normal();
  const auto& face_distance = a_plane.distance();
  const ScalarType MACHINE_EPSILON = machine_epsilon<ScalarType>();
    // If the intersection is not 2 parallel lines,
    // the intersection will be an ellipse with the center
    // on the Ox axis
  return PtBase<ScalarType>(
      -face_distance / safelyTiny(face_normal[0]),
      static_cast<ScalarType>(0),
      static_cast<ScalarType>(0));
}

template <class ScalarType>
inline NormalBase<ScalarType> getCylinderSurfaceNormal(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const PtBase<ScalarType>& a_pt) {
  return NormalBase<ScalarType>(
      static_cast<ScalarType>(0),
      static_cast<ScalarType>(2) * a_cylinder.b() * a_pt[1],
      static_cast<ScalarType>(2) * a_pt[2]);
};

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongLineOntoCylinder(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  const ScalarType a = (a_cylinder.b() * a_line[1] * a_line[1] +
                        a_line[2] * a_line[2]);
  const ScalarType b = static_cast<ScalarType>(2) * 
                       (a_cylinder.b() * a_line[1] * a_starting_pt[1] +
                        a_line[2] * a_starting_pt[2]);
  const ScalarType c = (a_cylinder.b() * a_starting_pt[1] * a_starting_pt[1] +
                        a_starting_pt[2] * a_starting_pt[2] -
                        a_cylinder.r());
  const StackVector<ScalarType, 2> solutions =
      solveQuadratic<ScalarType>(a, b, c);
  if (solutions.size() == 0) {
    std::cout << "No solution found for projection on cylinder" << a_line
              << a_starting_pt << std::endl;
  }
  assert(solutions.size() > 0);
  if (solutions.size() == 1) {
    return a_starting_pt + a_line * solutions[0];
  } else {
    if (abs(solutions[0]) < abs(solutions[1])) {
      return a_starting_pt + a_line * solutions[0];
    } else {
      return a_starting_pt + a_line * solutions[1];
    }
  }
}

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongHalfLineOntoCylinder(
    const AlignedCylinderBase<ScalarType>& a_cylinder,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  const ScalarType a = (a_cylinder.b() * a_line[1] * a_line[1] +
                        a_line[2] * a_line[2]);
  const ScalarType b = static_cast<ScalarType>(2) * 
                       (a_cylinder.b() * a_line[1] * a_starting_pt[1] +
                        a_line[2] * a_starting_pt[2]);
  const ScalarType c = (a_cylinder.b() * a_starting_pt[1] * a_starting_pt[1] +
                        a_starting_pt[2] * a_starting_pt[2] -
                        a_cylinder.r());
  const StackVector<ScalarType, 2> solutions =
      solveQuadratic<ScalarType>(a, b, c);
  if (solutions.size() == 0) {
    return PtBase<ScalarType>(static_cast<ScalarType>(DBL_MAX),
                              static_cast<ScalarType>(DBL_MAX),
                              static_cast<ScalarType>(DBL_MAX));
  }
  if (solutions.size() == 1) {
    if (solutions[0] < machine_epsilon<ScalarType>()) {
      return PtBase<ScalarType>(static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX));
    }
    return a_starting_pt + a_line * solutions[0];
  } else {
    if ((solutions[1] < static_cast<ScalarType>(0))) {
      return PtBase<ScalarType>(static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX),
                                static_cast<ScalarType>(DBL_MAX));

    } else {
      const ScalarType distance_along_line =
          solutions[0] > static_cast<ScalarType>(0)
              ? minimum(solutions[0], solutions[1])
              : maximum(solutions[0], solutions[1]);
      return a_starting_pt + a_line * distance_along_line;
    }
  }
}

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out, const CylinderBase<ScalarType>& a_cylinder) {
  const auto& datum = a_cylinder.getDatum();
  const auto& frame = a_cylinder.getReferenceFrame();
  const auto& aligned_cylinder = a_cylinder.getAlignedCylinder();

  out << "Datum: " << datum << '\n';
  out << "Frame: \n"
      << frame[0] << '\n'
      << frame[1] << '\n'
      << frame[2] << '\n';
  out << "Aligned cylinder: " << aligned_cylinder << '\n';
  out << "is always above? " << a_cylinder.isAlwaysAbove() << '\n';
  out << "is always below? " << a_cylinder.isAlwaysBelow() << '\n';
  return out;
}

}  // namespace IRL

#endif  // IRL_CYLINDER_RECONSTRUCTIONS_CYLINDER_TPP_
