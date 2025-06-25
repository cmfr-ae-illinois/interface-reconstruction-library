// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_VOLUME_FRACTION_MATCHING_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_VOLUME_FRACTION_MATCHING_TPP_

namespace IRL {

template <class CellType, class ReconstructionType>
inline void setDistanceToMatchVolumeFraction(
    const CellType& a_cell, const double a_volume_fraction,
    ReconstructionType* a_reconstruction,
    const double a_volume_fraction_tolerance) {
  if (wantPurelyInternal(a_volume_fraction) ||
      wantPurelyExternal(a_volume_fraction)) {
    setToPurePhaseReconstruction(a_volume_fraction, a_reconstruction);
  } else {
    setDistanceToMatchVolumeFractionPartialFill(a_cell, a_volume_fraction,
                                                a_reconstruction,
                                                a_volume_fraction_tolerance);
  }
}

template <class CellType, class VolumeFractionArrayType>
inline void setGroupDistanceToMatchVolumeFraction(
    const CellType& a_cell, const VolumeFractionArrayType& a_volume_fraction,
    PlanarSeparatorPathGroup* a_reconstruction,
    const double a_volume_fraction_tolerance) {
  setGroupDistanceToMatchVolumeFractionPartialFill(
      a_cell, a_volume_fraction, a_reconstruction, a_volume_fraction_tolerance);
}

template <class CellType, class PlanarType>
inline void setDistanceToMatchVolumeFractionPartialFill(
    const CellType& a_cell, const double a_volume_fraction,
    PlanarType* a_reconstruction, const double a_volume_fraction_tolerance) {
  assert(a_reconstruction != nullptr);
  runIterativeSolverForDistance(a_cell, a_volume_fraction, a_reconstruction,
                                a_volume_fraction_tolerance);
}

template <class CellType>
inline void setDistanceToMatchVolumeFractionPartialFill(
    const CellType& a_cell, const double a_volume_fraction,
    SeparatorVariant* a_reconstruction,
    const double a_volume_fraction_tolerance) {
  assert(a_reconstruction != nullptr);
  if (PlanarSeparator* separator =
          std::get_if<PlanarSeparator>(a_reconstruction)) {
    setDistanceToMatchVolumeFractionPartialFill(
        a_cell, a_volume_fraction, separator, a_volume_fraction_tolerance);
  } else if (Paraboloid* paraboloid =
                 std::get_if<Paraboloid>(a_reconstruction)) {
    const Pt& datum = paraboloid->getDatum();
    const ReferenceFrame& frame = paraboloid->getReferenceFrame();
    ProgressiveDistanceSolverParaboloid<CellType> solver(
        a_cell, a_volume_fraction, a_volume_fraction_tolerance, *paraboloid);
    paraboloid->setDatum(Pt(datum + solver.getDistance() * frame[2]));
  } else if (Cylinder* cylinder = std::get_if<Cylinder>(a_reconstruction)) {
    const Pt& datum = cylinder->getDatum();
    const ReferenceFrame& frame = cylinder->getReferenceFrame();
    ProgressiveRadiusSolverCylinder<CellType> solver(
        a_cell, a_volume_fraction, a_volume_fraction_tolerance, *cylinder);
    *cylinder = solver.getCylinder();
  } else {
    throw std::runtime_error(
        "Unknown SeparatorVariant type in volume matching routine");
  }
}

template <class ReconstructionType>
inline enable_if_t<not std::is_same_v<ReconstructionType, SeparatorVariant>,
                   void>
setDistanceToMatchVolumeFractionPartialFill(
    const RectangularCuboid& a_cell, const double a_volume_fraction,
    ReconstructionType* a_reconstruction,
    const double a_volume_fraction_tolerance) {
  assert(a_reconstruction != nullptr);
  if (a_reconstruction->getNumberOfPlanes() == 1) {
    (*a_reconstruction)[0].distance() = findDistanceOnePlane(
        a_cell, a_volume_fraction, (*a_reconstruction)[0].normal());
  } else {
    runIterativeSolverForDistance(a_cell, a_volume_fraction, a_reconstruction,
                                  a_volume_fraction_tolerance);
  }
}

template <class ReconstructionType>
inline enable_if_t<not std::is_same_v<ReconstructionType, SeparatorVariant>,
                   void>
setDistanceToMatchVolumeFractionPartialFill(
    const Tet& a_cell, const double a_volume_fraction,
    ReconstructionType* a_reconstruction,
    const double a_volume_fraction_tolerance) {
  assert(a_reconstruction != nullptr);
  if (a_reconstruction->getNumberOfPlanes() == 1) {
    (*a_reconstruction)[0].distance() = findDistanceOnePlane(
        a_cell, a_volume_fraction, (*a_reconstruction)[0].normal());
  } else {
    runIterativeSolverForDistance(a_cell, a_volume_fraction, a_reconstruction,
                                  a_volume_fraction_tolerance);
  }
}

template <class CellType, class VolumeFractionArrayType>
inline void setGroupDistanceToMatchVolumeFractionPartialFill(
    const CellType& a_cell, const VolumeFractionArrayType& a_volume_fraction,
    PlanarSeparatorPathGroup* a_reconstruction,
    const double a_volume_fraction_tolerance) {
  runIterativeSolverForDistance(a_cell, a_volume_fraction, a_reconstruction,
                                a_volume_fraction_tolerance);
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_VOLUME_FRACTION_MATCHING_TPP_