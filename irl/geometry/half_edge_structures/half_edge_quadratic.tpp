// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_QUADRATIC_TPP_
#define IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_QUADRATIC_TPP_

namespace IRL {

template <class VertexType>
inline void doubleLinkHalfEdges(
    HalfEdgeQuadratic<VertexType>* a_starting_half_edge,
    HalfEdgeQuadratic<VertexType>* a_ending_half_edge) {
  a_starting_half_edge->setNextHalfEdge(a_ending_half_edge);
  a_ending_half_edge->setPreviousHalfEdge(a_starting_half_edge);
}

template <class VertexType>
inline void setMutualOpposites(HalfEdgeQuadratic<VertexType>* a_half_edge_0,
                               HalfEdgeQuadratic<VertexType>* a_half_edge_1) {
  a_half_edge_0->setOppositeHalfEdge(a_half_edge_1);
  a_half_edge_1->setOppositeHalfEdge(a_half_edge_0);
}

template <class VertexType>
inline HalfEdgeQuadratic<VertexType>::HalfEdgeQuadratic(void)
    : previous_m(nullptr),
      next_m(nullptr),
      opposite_m(nullptr),
      end_point_m(nullptr),
      corresponding_face_m(nullptr) {}

template <class VertexType>
inline HalfEdgeQuadratic<VertexType>::HalfEdgeQuadratic(
    VertexType* a_vertex, HalfEdgeQuadratic* a_previous,
    HalfEdgeQuadratic* a_next, FaceQuadratic<HalfEdgeQuadratic>* a_face)
    : previous_m(a_previous),
      next_m(a_next),
      opposite_m(nullptr),
      end_point_m(a_vertex),
      corresponding_face_m(a_face) {}

template <class VertexType>
inline void HalfEdgeQuadratic<VertexType>::setPreviousHalfEdge(
    HalfEdgeQuadratic* a_previous) {
  previous_m = a_previous;
}
template <class VertexType>
inline HalfEdgeQuadratic<VertexType>*
HalfEdgeQuadratic<VertexType>::getPreviousHalfEdge(void) {
  return previous_m;
}
template <class VertexType>
inline const HalfEdgeQuadratic<VertexType>*
HalfEdgeQuadratic<VertexType>::getPreviousHalfEdge(void) const {
  return previous_m;
}

template <class VertexType>
inline void HalfEdgeQuadratic<VertexType>::setNextHalfEdge(
    HalfEdgeQuadratic* a_next) {
  next_m = a_next;
}
template <class VertexType>
inline HalfEdgeQuadratic<VertexType>*
HalfEdgeQuadratic<VertexType>::getNextHalfEdge(void) {
  return next_m;
}
template <class VertexType>
inline const HalfEdgeQuadratic<VertexType>*
HalfEdgeQuadratic<VertexType>::getNextHalfEdge(void) const {
  return next_m;
}

template <class VertexType>
inline void HalfEdgeQuadratic<VertexType>::setOppositeHalfEdge(
    HalfEdgeQuadratic* a_opposite) {
  opposite_m = a_opposite;
}
template <class VertexType>
inline HalfEdgeQuadratic<VertexType>*
HalfEdgeQuadratic<VertexType>::getOppositeHalfEdge(void) {
  return opposite_m;
}
template <class VertexType>
inline const HalfEdgeQuadratic<VertexType>*
HalfEdgeQuadratic<VertexType>::getOppositeHalfEdge(void) const {
  return opposite_m;
}

template <class VertexType>
inline void HalfEdgeQuadratic<VertexType>::setFace(
    FaceQuadratic<HalfEdgeQuadratic<VertexType>>* a_face) {
  corresponding_face_m = a_face;
}
template <class VertexType>
inline FaceQuadratic<HalfEdgeQuadratic<VertexType>>*
HalfEdgeQuadratic<VertexType>::getFace(void) {
  return corresponding_face_m;
}
template <class VertexType>
inline const FaceQuadratic<HalfEdgeQuadratic<VertexType>>*
HalfEdgeQuadratic<VertexType>::getFace(void) const {
  return corresponding_face_m;
}

template <class VertexType>
inline void HalfEdgeQuadratic<VertexType>::setVertex(VertexType* a_vertex) {
  end_point_m = a_vertex;
}
template <class VertexType>
inline VertexType* HalfEdgeQuadratic<VertexType>::getVertex(void) {
  return end_point_m;
}
template <class VertexType>
inline const VertexType* HalfEdgeQuadratic<VertexType>::getVertex(void) const {
  return end_point_m;
}

template <class VertexType>
inline VertexType* HalfEdgeQuadratic<VertexType>::getPreviousVertex(void) {
  assert(previous_m != nullptr);
  return previous_m->getVertex();
}
template <class VertexType>
inline const VertexType* HalfEdgeQuadratic<VertexType>::getPreviousVertex(
    void) const {
  assert(previous_m != nullptr);
  return previous_m->getVertex();
}

template <class PtType>
inline VertexQuadratic<PtType>::VertexQuadratic(void)
    : vertex_location_m(PtBase<typename PtType::value_type>(
          static_cast<typename PtType::value_type>(0),
          static_cast<typename PtType::value_type>(0),
          static_cast<typename PtType::value_type>(0))),
      half_edge_m(nullptr),
      distance_m(DBL_MAX),
      is_clipped_m(false),
      needs_to_seek_m(false) {}

template <class PtType>
inline VertexQuadratic<PtType>::VertexQuadratic(const PtType& a_location)
    : vertex_location_m(a_location),
      half_edge_m(nullptr),
      distance_m(DBL_MAX),
      is_clipped_m(false),
      needs_to_seek_m(false) {}

template <class PtType>
inline void VertexQuadratic<PtType>::setHalfEdge(
    HalfEdgeQuadratic<VertexQuadratic>* a_half_edge) {
  half_edge_m = a_half_edge;
}
template <class PtType>
inline HalfEdgeQuadratic<VertexQuadratic<PtType>>*
VertexQuadratic<PtType>::getHalfEdge(void) {
  return half_edge_m;
}
template <class PtType>
inline const HalfEdgeQuadratic<VertexQuadratic<PtType>>*
VertexQuadratic<PtType>::getHalfEdge(void) const {
  return half_edge_m;
}

template <class PtType>
inline PtType& VertexQuadratic<PtType>::getLocation(void) {
  return vertex_location_m;
}

template <class PtType>
inline const PtType& VertexQuadratic<PtType>::getLocation(void) const {
  return vertex_location_m;
}

template <class PtType>
inline void VertexQuadratic<PtType>::setLocation(const PtType& a_location) {
  vertex_location_m = a_location;
}

template <class PtType>
inline void VertexQuadratic<PtType>::calculateDistanceToPlane(
    const PlaneBase<typename PtType::value_type>& a_plane) {
  distance_m = a_plane.signedDistanceToPoint(this->getLocation());
}
template <class PtType>
inline typename PtType::value_type VertexQuadratic<PtType>::getDistance(
    void) const {
  return distance_m;
}

template <class PtType>
inline void VertexQuadratic<PtType>::markToBeClipped(void) {
  is_clipped_m = true;
}
template <class PtType>
inline void VertexQuadratic<PtType>::markToBeNotClipped(void) {
  is_clipped_m = false;
}
template <class PtType>
inline void VertexQuadratic<PtType>::setClip(const bool a_is_clipped) {
  is_clipped_m = a_is_clipped;
}

template <class PtType>
inline bool VertexQuadratic<PtType>::isClipped(void) const {
  return is_clipped_m;
}
template <class PtType>
inline bool VertexQuadratic<PtType>::isNotClipped(void) const {
  return !is_clipped_m;
}
template <class PtType>
inline void VertexQuadratic<PtType>::setAsUnnecessaryToSeek(void) {
  needs_to_seek_m = false;
}
template <class PtType>
inline void VertexQuadratic<PtType>::setToSeek(void) {
  needs_to_seek_m = true;
}
template <class PtType>
inline bool VertexQuadratic<PtType>::needsToSeek(void) const {
  return needs_to_seek_m;
}
template <class PtType>
inline bool VertexQuadratic<PtType>::doesNotNeedToSeek(void) const {
  return !needs_to_seek_m;
}

template <class PtType>
bool VertexQuadratic<PtType>::checkValidHalfEdgeCycle(void) const {
  assert(half_edge_m != nullptr);
  const HalfEdgeQuadratic<VertexQuadratic>* current_half_edge =
      this->getHalfEdge();
  // To avoid hanging if enters another complete circuit
  // and can never return to start.
  constexpr UnsignedIndex_t max_loop_size = 1000;
  UnsignedIndex_t number_of_vertices = 0;
  do {
    if (number_of_vertices > max_loop_size) {
      std::cout
          << "Half edge for vertex did not cycle around vertex back to itself. "
          << std::endl;
      return false;
    }
    ++number_of_vertices;
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != this->getHalfEdge());
  return true;
}

template <class HalfEdgeType>
inline FaceQuadratic<HalfEdgeType>::FaceQuadratic(void)
    : Face<HalfEdgeType>(),
      face_plane_m(),
      intersections_m(0),
      edge_parallel_intersections_m(0),
      is_triangle_m(false) {}

template <class HalfEdgeType>
inline FaceQuadratic<HalfEdgeType>::FaceQuadratic(
    HalfEdgeType* a_starting_half_edge)
    : Face<HalfEdgeType>(a_starting_half_edge),
      face_plane_m(),
      intersections_m(0),
      edge_parallel_intersections_m(0),
      is_triangle_m(false) {}

template <class HalfEdgeType>
inline void FaceQuadratic<HalfEdgeType>::setPlane(
    const PlaneBase<typename HalfEdgeType::value_type>& a_plane) {
  face_plane_m = a_plane;
}
template <class HalfEdgeType>
inline const PlaneBase<typename HalfEdgeType::value_type>&
FaceQuadratic<HalfEdgeType>::getPlane(void) const {
  return face_plane_m;
}

template <class HalfEdgeType>
void FaceQuadratic<HalfEdgeType>::clearIntersections(void) {
  intersections_m = 0;
  edge_parallel_intersections_m = 0;
}

template <class HalfEdgeType>
void FaceQuadratic<HalfEdgeType>::addIntersection(void) {
  ++intersections_m;
}

template <class HalfEdgeType>
void FaceQuadratic<HalfEdgeType>::addDoubleIntersection(void) {
  intersections_m += 2;
}

template <class HalfEdgeType>
UnsignedIndex_t FaceQuadratic<HalfEdgeType>::getNumberOfIntersections(
    void) const {
  return intersections_m;
}

template <class HalfEdgeType>
void FaceQuadratic<HalfEdgeType>::addEdgeParallelIntersection(void) {
  ++edge_parallel_intersections_m;
}

template <class HalfEdgeType>
void FaceQuadratic<HalfEdgeType>::addEdgeParallelIntersections(
    const UnsignedIndex_t a_intersections) {
  edge_parallel_intersections_m += a_intersections;
}

template <class HalfEdgeType>
UnsignedIndex_t
FaceQuadratic<HalfEdgeType>::getNumberOfEdgeParallelIntersections(void) const {
  return edge_parallel_intersections_m;
}

template <class HalfEdgeType>
inline void FaceQuadratic<HalfEdgeType>::setAsTriangle(void) {
  is_triangle_m = true;
}

template <class HalfEdgeType>
inline void FaceQuadratic<HalfEdgeType>::setAsNotTriangle(void) {
  is_triangle_m = false;
}

template <class HalfEdgeType>
inline bool FaceQuadratic<HalfEdgeType>::isTriangle(void) const {
  return is_triangle_m;
}


template <class PtType>
inline std::ostream& operator<<(
    std::ostream& out, const VertexQuadratic<PtType>& a_vertex) {
  const auto& st_point = a_vertex.getLocation();

  out << "position " << st_point << '\n';
  out << "distance : " << a_vertex.getDistance() << '\n';
  out << "is_clipped ? " << a_vertex.isClipped() << '\n';
  out << "needs_to_seek_m ? " << a_vertex.needsToSeek() << '\n';
  return out;
}

template <class HalfEdgeType>
inline std::ostream& operator<<(
    std::ostream& out, const FaceQuadratic<HalfEdgeType>& a_face) {
  const auto& a_plane = a_face.getPlane();
  const auto& a_normal = a_plane.normal();

  out << std::setprecision(15);

  out << "face of equation : " << 
         a_normal[0] << " * x + " <<
         a_normal[1] << " * y + " <<
         a_normal[2] << " * z - " <<
         a_plane.distance() << " = 0\n";
  out << "nb intersection " << a_face.getNumberOfIntersections() << '\n';
  out << "nb parra intersection " << a_face.getNumberOfEdgeParallelIntersections() << '\n';
  out << "triangle ? " << a_face.isTriangle() << '\n';
  return out;
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_QUADRATIC_TPP_
