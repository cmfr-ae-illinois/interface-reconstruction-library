// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// NOTE: THIS SHOULD BE WRITTEN WITH JUST THE OTHER HALF EDGE, WHERE
// A NEW FACE CLASS IS MADE AND EVERYTHING TEMPLATED ON IT. THIS
// HOWEVER LEADS TO A RECURSIVE INSTANTIATION OF TEMPLATES, SINCE
// THE HALF EDGE CLASS AND FACE CLASS RELY ON ONE ANOTHER.
// FOR NOW, JUST HARD CODE DIFFERENT CLASSES.

#ifndef IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_QUADRATIC_H_
#define IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_QUADRATIC_H_

#include <float.h>
#include <utility>

#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/half_edge_structures/half_edge.h"

namespace IRL {

// Forward declarations.
template <class PtType>
class VertexQuadratic;

template <class VertexType>
class HalfEdgeQuadratic;

template <class HalfEdgeType>
class FaceQuadratic;

template <class VertexType>
inline void doubleLinkHalfEdges(
    HalfEdgeQuadratic<VertexType>* a_starting_half_edge,
    HalfEdgeQuadratic<VertexType>* a_ending_half_edge);

template <class VertexType>
inline void setMutualOpposites(
    HalfEdgeQuadratic<VertexType>* a_starting_half_edge,
    HalfEdgeQuadratic<VertexType>* a_ending_half_edge);

template <class VertexType>
class HalfEdgeQuadratic {
 public:
  using vertex_type = VertexType;
  using value_type = typename vertex_type::value_type;

  HalfEdgeQuadratic(void);

  HalfEdgeQuadratic(VertexType* a_vertex, HalfEdgeQuadratic* a_previous,
                     HalfEdgeQuadratic* a_next,
                     FaceQuadratic<HalfEdgeQuadratic>* a_face);

  HalfEdgeQuadratic(const HalfEdgeQuadratic& a_other) = default;
  HalfEdgeQuadratic& operator=(const HalfEdgeQuadratic& a_other) = default;

  void setPreviousHalfEdge(HalfEdgeQuadratic* a_previous);
  HalfEdgeQuadratic* getPreviousHalfEdge(void);
  const HalfEdgeQuadratic* getPreviousHalfEdge(void) const;

  void setNextHalfEdge(HalfEdgeQuadratic* a_next);
  HalfEdgeQuadratic* getNextHalfEdge(void);
  const HalfEdgeQuadratic* getNextHalfEdge(void) const;

  void setOppositeHalfEdge(HalfEdgeQuadratic* a_opposite);
  HalfEdgeQuadratic* getOppositeHalfEdge(void);
  const HalfEdgeQuadratic* getOppositeHalfEdge(void) const;

  void setFace(FaceQuadratic<HalfEdgeQuadratic>* a_face);
  FaceQuadratic<HalfEdgeQuadratic>* getFace(void);
  const FaceQuadratic<HalfEdgeQuadratic>* getFace(void) const;

  void setVertex(VertexType* a_vertex);
  VertexType* getVertex(void);
  const VertexType* getVertex(void) const;

  VertexType* getPreviousVertex(void);
  const VertexType* getPreviousVertex(void) const;

  ~HalfEdgeQuadratic(void) = default;

 private:
  HalfEdgeQuadratic* previous_m;
  HalfEdgeQuadratic* next_m;
  HalfEdgeQuadratic* opposite_m;
  VertexType* end_point_m;
  FaceQuadratic<HalfEdgeQuadratic>* corresponding_face_m;
};

template <class PtType>
class VertexQuadratic {
 public:
  using pt_type = PtType;
  using value_type = typename pt_type::value_type;

  VertexQuadratic(void);

  explicit VertexQuadratic(const PtType& a_location);

  void setHalfEdge(HalfEdgeQuadratic<VertexQuadratic>* a_half_edge);
  HalfEdgeQuadratic<VertexQuadratic>* getHalfEdge(void);
  const HalfEdgeQuadratic<VertexQuadratic>* getHalfEdge(void) const;

  pt_type& getLocation(void);
  const pt_type& getLocation(void) const;
  void setLocation(const PtType& a_location);

  void calculateDistanceToPlane(const PlaneBase<value_type>& a_plane);
  value_type getDistance(void) const;
  value_type setDistance(void) const;

  void markToBeClipped(void);
  void markToBeNotClipped(void);
  void setClip(const bool a_is_clipped);
  bool isClipped(void) const;
  bool isNotClipped(void) const;

  void setAsUnnecessaryToSeek(void);
  void setToSeek(void);
  bool needsToSeek(void) const;
  bool doesNotNeedToSeek(void) const;

  inline bool checkValidHalfEdgeCycle(void) const;

 private:
  PtType vertex_location_m;
  HalfEdgeQuadratic<VertexQuadratic>*
      half_edge_m;  // HalfEdgeQuadratic that ends at this vertex
  value_type distance_m;
  bool is_clipped_m = false;
  bool needs_to_seek_m = false;
};

template <class HalfEdgeType>
class FaceQuadratic : public Face<HalfEdgeType> {
 public:
  using half_edge_type = HalfEdgeType;
  using value_type = typename half_edge_type::value_type;

  FaceQuadratic(void);

  explicit FaceQuadratic(HalfEdgeType* a_starting_half_edge);

  void setPlane(const PlaneBase<value_type>& a_plane);
  const PlaneBase<value_type>& getPlane(void) const;

  void clearIntersections(void);
  void addIntersection(void);
  void addDoubleIntersection(void);
  UnsignedIndex_t getNumberOfIntersections(void) const;

  /// \brief This is a subset of the intersections whose edge is directly
  /// parallel to the tangent of the paraboloid at that point.
  void addEdgeParallelIntersection(void);
  void addEdgeParallelIntersections(const UnsignedIndex_t a_intersections);
  UnsignedIndex_t getNumberOfEdgeParallelIntersections(void) const;
  void setAsTriangle(void);
  void setAsNotTriangle(void);
  bool isTriangle(void) const;

 private:
  PlaneBase<value_type> face_plane_m;
  UnsignedIndex_t intersections_m;
  UnsignedIndex_t edge_parallel_intersections_m;
  bool is_triangle_m = false;
};

template <class PtType>
inline std::ostream& operator<<(std::ostream& out,
                                const VertexQuadratic<PtType>& a_vertex);

template <class HalfEdgeType>
inline std::ostream& operator<<(std::ostream& out,
                                const FaceQuadratic<HalfEdgeType>& a_face);

}  // namespace IRL

#include "irl/geometry/half_edge_structures/half_edge_quadratic.tpp"

#endif  // IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_QUADRATIC_H_
