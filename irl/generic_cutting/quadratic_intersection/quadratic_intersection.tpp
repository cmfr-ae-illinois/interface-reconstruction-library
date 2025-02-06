// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_TPP_
#define IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_TPP_

#include <float.h>
#include <cassert>
#include <cmath>
#include <random>
#include <type_traits>
#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_initializer.tpp"
#include "irl/generic_cutting/quadratic_intersection/moment_contributions.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/geometry/half_edge_structures/brep_to_half_edge.h"
#include "irl/helpers/mymath.h"
#include "irl/moments/general_moments.h"
#include "irl/quadratic_reconstruction/rational_bezier_arc.h"
#include "irl/generic_cutting/quadratic_intersection/quadratic_intersection.h"

#define NUMERICAL_INTEGRATION

namespace IRL {

// ScalarType refers to the type of scalar container used. This can be:
// - double
// - Quad_t (i.e., __float128)
// or a scalar with embedded derivatives:
// - ScalarWithGradient<FloatType,GradientType>
//
// The templeted type float_type<ScalarType> is the primitive floating-point
// type used by the scalar container, i.e., double or __float128

template <class ScalarType>
inline ScalarType distance_epsilon(void) {
  using FloatType = float_type<ScalarType>;
  if constexpr (std::is_same_v<FloatType, double>) {
    return ScalarType(1.0e2 * DBL_EPSILON);
  } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
    return ScalarType(1.0e6q * FLT128_EPSILON);
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain distance epsilon for unknown float type");
  }
}

template <class ScalarType>
inline const ScalarType angle_epsilon(void) {
  using FloatType = float_type<ScalarType>;
  if constexpr (std::is_same_v<FloatType, double>) {
    return ScalarType(1.0e6 * DBL_EPSILON);
  } else if constexpr (std::is_same_v<FloatType, Quad_t>) {
    return ScalarType(1.0e10q * FLT128_EPSILON);
  } else {
    static_assert(sizeof(ScalarType) == 0,
                  "Trying to obtain angle epsilon for unknown float type");
  }
}

/**************** Calculate line contribution ******************/
// Starts from an entry, returns the exit that is reached.
template <class ReturnType, class ScalarType, class HalfEdgeType, class PtType,
          class NormalType>
ReturnType computeUnclippedSegmentType1Contribution(
    const PtType& a_ref_pt, const HalfEdgeType a_entry_half_edge,
    HalfEdgeType& a_exit_half_edge, bool* skip_first,
    const NormalType& a_face_normal, const UnsignedIndex_t a_proj_dir) {
  // Defining constants and types
  using ReturnScalarType = typename ReturnType::value_type;
  const ScalarType ZERO = ScalarType(0);

  // Initialise moments to 0
  ReturnType full_moments = ReturnType::fromScalarConstant(ReturnScalarType(0));

  assert(a_entry_half_edge->getPreviousVertex()->isClipped() ||
         a_entry_half_edge->getPreviousVertex()->doesNotNeedToSeek());
  assert(a_entry_half_edge->getNextHalfEdge()->getVertex()->isNotClipped());

  // Traverse face and add type 1 moment contribution of unclipped edge
  auto current_half_edge = a_entry_half_edge->getNextHalfEdge();
  auto vertex = current_half_edge->getVertex();
  auto prev_pt = current_half_edge->getPreviousVertex()->getLocation();

  while (true) {
    vertex = current_half_edge->getVertex();
    const auto& curr_pt = vertex->getLocation();
    full_moments += computeType1Contribution<ReturnType, ScalarType>(
        a_ref_pt, prev_pt, curr_pt, skip_first, false, a_face_normal,
        a_proj_dir);
    if (vertex->needsToSeek()) {
      prev_pt = curr_pt;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } else {
      assert(!(
          current_half_edge->getPreviousVertex()->isClipped() ||
          (current_half_edge->getPreviousVertex()->doesNotNeedToSeek() &&
           current_half_edge->getNextHalfEdge()->getVertex()->isNotClipped())));
      break;
    }
  }
  // Return the last exit intersection encounrtered
  a_exit_half_edge = current_half_edge;
  return full_moments;
}

/******************* Place one intersection on segment **********************/
template <class PtType, class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType>
void placeSingleIntercept(const PtType& a_intersection_location,
                          HalfEdgeType* a_half_edge_with_intersection,
                          SegmentedHalfEdgePolytopeType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope) {
  auto first_intersection_vertex = a_complete_polytope->getNewVertex(
      typename SegmentedHalfEdgePolytopeType::vertex_type(
          a_intersection_location));
  first_intersection_vertex->setToSeek();
  first_intersection_vertex->markToBeNotClipped();
  a_polytope->addVertex(first_intersection_vertex);
  HalfEdgeType* new_half_edge = separateIntersectedHalfEdge(
      first_intersection_vertex, a_half_edge_with_intersection, a_polytope,
      a_complete_polytope);
  createOppositeHalfEdgeFromIntersection(a_half_edge_with_intersection,
                                         new_half_edge, a_complete_polytope);
}

/******************* Place two intersection on segment **********************/
template <class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType, class VertexType>
void placeDoubleIntercept(
    const StackVector<VertexType, 2>& a_intersection_location,
    HalfEdgeType* a_half_edge_with_intersection,
    SegmentedHalfEdgePolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope) {
  assert(a_intersection_location.size() == 2);

  // Need to place furthest vertex first so that
  // a_half_edge_with_intersection remains attached to same vertex
  // it started on. Also want to keep property that the half edge
  // the vertex stores has the previousVertex also unclipped. This
  // means for the first (furthest) vertex, need to reference the
  // opposite half edge of the current one. This enables more
  // efficient face truncation in the next phase of the algorithm,
  // since new vertices will always have a half edge going from
  // unclipped->new->clipped.

  placeSingleIntercept(a_intersection_location[0],
                       a_half_edge_with_intersection, a_polytope,
                       a_complete_polytope);
  a_half_edge_with_intersection->getPreviousVertex()->setHalfEdge(
      a_half_edge_with_intersection->getOppositeHalfEdge());
  placeSingleIntercept(a_intersection_location[1],
                       a_half_edge_with_intersection, a_polytope,
                       a_complete_polytope);
}

// This algorithm is based on the method shared from the website below, and
// from several stackoverflow posts citing that website.
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdgeXY(
    const std::array<ScalarType, 2>& a_test_pt,
    const PtBase<ScalarType>& a_vertex_0,
    const PtBase<ScalarType>& a_vertex_1) {
  if ((a_test_pt[1] > a_vertex_0[1]) == (a_test_pt[1] > a_vertex_1[1])) {
    return false;  // Projected ray never intersects edge.
  }
  const ScalarType location_of_intersection_along_ray =
      (a_vertex_0[0] - a_vertex_1[0]) * (a_test_pt[1] - a_vertex_1[1]) /
          (a_vertex_0[1] - a_vertex_1[1]) +
      a_vertex_1[0];
  // Intersection was to the right if
  // location_of_intersection_along_ray is greater.
  return a_test_pt[0] < location_of_intersection_along_ray;
}

// For the cylinder case, i wont use x-y coordinate but y-z coordinate
template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdgeYZ(
    const std::array<ScalarType, 2>& a_test_pt,
    const PtBase<ScalarType>& a_vertex_0,
    const PtBase<ScalarType>& a_vertex_1) {
  if ((a_test_pt[2] > a_vertex_0[2]) == (a_test_pt[2] > a_vertex_1[2])) {
    return false;  // Projected ray never intersects edge.
  }
  const ScalarType location_of_intersection_along_ray =
      (a_vertex_0[1] - a_vertex_1[1]) * (a_test_pt[1] - a_vertex_1[2]) /
          (a_vertex_0[2] - a_vertex_1[2]) +
      a_vertex_1[1];
  // Intersection was to the right if
  // location_of_intersection_along_ray is greater.
  return a_test_pt[0] < location_of_intersection_along_ray;
}

// This algorithm is based on the method shared from the website below, and
// from several stackoverflow posts citing that website.
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdgeWithComponent(
    const PtBase<ScalarType>& a_test_pt, const PtBase<ScalarType>& a_vertex_0,
    const PtBase<ScalarType>& a_vertex_1, const UnsignedIndex_t a_index) {
  const UnsignedIndex_t id0 = a_index;
  const UnsignedIndex_t id1 = (a_index + 1) % 3;
  if ((a_test_pt[id1] > a_vertex_0[id1]) ==
      (a_test_pt[id1] > a_vertex_1[id1])) {
    return false;  // Projected ray never intersects edge.
  }
  const ScalarType location_of_intersection_along_ray =
      (a_vertex_0[id0] - a_vertex_1[id0]) * (a_test_pt[id1] - a_vertex_1[id1]) /
          (a_vertex_0[id1] - a_vertex_1[id1]) +
      a_vertex_1[id0];
  // Intersection was to the right if
  // location_of_intersection_along_ray is greater.
  return a_test_pt[id0] < location_of_intersection_along_ray;
}

// Note: This essentially abandons the half edges and vertices
// in a_complete_polytope, doing nothing to reclaim their memory usage.
// This should not be a problem unless resetPolyhedron is called
// many times for the same polyhedron intersection operation
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, void>
resetPolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                HalfEdgePolytopeType* a_complete_polytope) {
  UnsignedIndex_t original_verts = 0;
  for (UnsignedIndex_t v = 0; v < a_polytope->getNumberOfVertices(); ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    if (vertex.needsToSeek()) {
      ++original_verts;
    } else {
      // New vertex from intersection, remove it and patch
      // half-edges
      auto current_edge = vertex.getHalfEdge();
      auto original_edge = current_edge->getNextHalfEdge();
      doubleLinkHalfEdges(current_edge->getPreviousHalfEdge(), original_edge);
      original_edge->getFace()->setStartingHalfEdge(original_edge);
      auto opposite_edge = current_edge->getOppositeHalfEdge();
      doubleLinkHalfEdges(
          original_edge->getOppositeHalfEdge()->getPreviousHalfEdge(),
          opposite_edge);
      opposite_edge->getFace()->setStartingHalfEdge(opposite_edge);
      setMutualOpposites(original_edge, opposite_edge);
    }
  }
  a_polytope->setNumberOfVertices(original_verts);
  for (UnsignedIndex_t f = 0; f < a_polytope->getNumberOfFaces(); ++f) {
    (*a_polytope)[f]->clearIntersections();
  }
}

/********************** Move poly to find stable topology ****************/
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType,
          class SurfaceOutputType>
void nudgePolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                     HalfEdgePolytopeType* a_complete_polytope,
                     const UnsignedIndex_t a_nudge_iter,
                     SurfaceOutputType* a_surface) {
  using PtType = typename SegmentedHalfEdgePolyhedronType::vertex_type::pt_type;
  using ScalarType = typename PtType::value_type;
  using FloatType = float_type<ScalarType>;
  static_assert(std::is_same_v<FloatType, Quad_t>);

  // Create a random number generator and seed it with the number of prior
  // nudge iterations (for reproductibility)
  std::random_device rd;
  std::mt19937 gen(a_nudge_iter);
  std::uniform_real_distribution distr(-1.0, 1.0);

  // This is a-hoc but works well
  const ScalarType nudge_epsilon = ScalarType(
      1.0e10q * powq(static_cast<FloatType>(a_nudge_iter + 1), 1.05q) *
      distance_epsilon<FloatType>());

  // Compute polytope center to rotate about it
  auto center = a_polytope->calculateCentroid();
  auto converted_center = PtBase<ScalarType>(
      ScalarType(center[0]), ScalarType(center[1]), ScalarType(center[2]));
  // Get random translation and rotations
  const auto nudge_translation =
      nudge_epsilon * NormalBase<ScalarType>(ScalarType(distr(gen)),
                                             ScalarType(distr(gen)),
                                             ScalarType(distr(gen)));
  const auto nudge_rotation =
      2.0q * nudge_epsilon *
      NormalBase<ScalarType>(ScalarType(distr(gen)), ScalarType(distr(gen)),
                             ScalarType(distr(gen)));

  // Compute rotated reference frame
  ReferenceFrameBase<ScalarType> frame(NormalBase<ScalarType>(1, 0, 0),
                                       NormalBase<ScalarType>(0, 1, 0),
                                       NormalBase<ScalarType>(0, 0, 1));
  UnitQuaternionBase<ScalarType> x_rotation(nudge_rotation[0], frame[0]);
  UnitQuaternionBase<ScalarType> y_rotation(nudge_rotation[1], frame[1]);
  UnitQuaternionBase<ScalarType> z_rotation(nudge_rotation[2], frame[2]);
  auto total_rotation = x_rotation * y_rotation * z_rotation;
  total_rotation.normalize();
  frame = total_rotation * frame;

  // Shake polytope with random translation and rotation. These are well
  // below DBL_EPSILON so should not impact the accuracy of the moments in
  // DP
  for (UnsignedIndex_t v = 0; v < a_polytope->getNumberOfVertices(); ++v) {
    const PtBase<ScalarType> original_pt =
        a_polytope->getVertex(v)->getLocation().getPt() - converted_center;
    PtBase<ScalarType> projected_location(0, 0, 0);
    auto& pt = projected_location.getPt();
    for (UnsignedIndex_t n = 0; n < 3; ++n) {
      pt[n] = frame[n] * original_pt;
    }
    pt += converted_center + nudge_translation;
    a_polytope->getVertex(v)->setLocation(projected_location);
  }

  // Clear surface from any arcs that might have already been added
  if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
    a_surface->clear();
  }
}

// Other helper functions


template <class MappingContainer>
UnsignedIndex_t positionInMapping(const std::vector<MappingContainer> a_mapping,
                                  const MappingContainer a_element) {
  const auto position =
      std::find(a_mapping.begin(), a_mapping.end(), a_element);
  assert(position != a_mapping.end());
  return static_cast<UnsignedIndex_t>(position - a_mapping.begin());
}

template <class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType>
void fullCopyOfCompletePolytope(
    SegmentedHalfEdgePolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    HalfEdgePolytopeType* a_copy_polytope) {
  // Convert polytope type
  using face_type = typename SegmentedHalfEdgePolytopeType::face_type;
  using half_edge_type = typename face_type::half_edge_type;
  using vertex_type = typename half_edge_type::vertex_type;
  using pt_type = typename vertex_type::pt_type;
  using scalar_type = typename pt_type::value_type;

  assert(a_polytope->checkValidHalfEdgeStructure());

  // Calculate number of half edges, vertices and faces
  const UnsignedIndex_t number_of_vertices = a_polytope->getNumberOfVertices();
  const UnsignedIndex_t number_of_faces = a_polytope->getNumberOfFaces();
  UnsignedIndex_t number_of_half_edges = 0;
  std::vector<face_type*> face_mapping;
  face_mapping.resize(number_of_faces);
  std::vector<half_edge_type*> hald_edge_mapping;
  hald_edge_mapping.resize(0);
  std::vector<vertex_type*> vertex_mapping;
  vertex_mapping.resize(number_of_vertices);

  UnsignedIndex_t f = 0;
  for (const auto& face : (*a_polytope)) {
    face_mapping[f++] = face;
    const auto starting_half_edge = face->getStartingHalfEdge();
    auto current_half_edge = starting_half_edge;
    do {
      hald_edge_mapping.push_back(current_half_edge);
      number_of_half_edges++;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != starting_half_edge);
  }
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    vertex_mapping[v] = a_polytope->getVertex(v);
  }

  a_copy_polytope->resize(number_of_half_edges, number_of_vertices,
                               number_of_faces);
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    const auto old_pt = a_polytope->getVertex(v)->getLocation();
    auto& new_pt = a_copy_polytope->getVertex(v).getLocation();
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      new_pt[d] = old_pt[d];
    }
  }

  UnsignedIndex_t e = 0;
  for (const auto& half_edge_DP : hald_edge_mapping) {
    // Pointer to QP halfedge
    half_edge_type* half_edge_QP =
        &a_copy_polytope->getHalfEdge(e++);

    // Find QP face on which it lies
    const auto face_DP = half_edge_DP->getFace();
    const UnsignedIndex_t index_face = positionInMapping(face_mapping, face_DP);
    face_type* face_QP = &a_copy_polytope->getFace(index_face);

    // Find end-point vertex
    const auto vertex_DP = half_edge_DP->getVertex();
    const UnsignedIndex_t index_vertex =
        positionInMapping(vertex_mapping, vertex_DP);
    vertex_type* vertex_QP =
        &a_copy_polytope->getVertex(index_vertex);

    // Find previous, next and opposite QP halfedges
    const auto previous_DP = half_edge_DP->getPreviousHalfEdge();
    const UnsignedIndex_t index_previous =
        positionInMapping(hald_edge_mapping, previous_DP);
    half_edge_type* previous_QP =
        &a_copy_polytope->getHalfEdge(index_previous);
    const auto next_DP = half_edge_DP->getNextHalfEdge();
    const UnsignedIndex_t index_next =
        positionInMapping(hald_edge_mapping, next_DP);
    half_edge_type* next_QP =
        &a_copy_polytope->getHalfEdge(index_next);
    const auto opposite_DP = half_edge_DP->getOppositeHalfEdge();
    const UnsignedIndex_t index_opposite =
        positionInMapping(hald_edge_mapping, opposite_DP);
    half_edge_type* opposite_QP =
        &a_copy_polytope->getHalfEdge(index_opposite);

    *half_edge_QP =
        half_edge_type(vertex_QP, previous_QP, next_QP, face_QP);
    vertex_QP->setHalfEdge(half_edge_QP);
    half_edge_QP->setOppositeHalfEdge(opposite_QP);
    face_QP->setStartingHalfEdge(half_edge_QP);
  }
}

template <class DoubleSegmentedHalfEdgePolytopeType,
          class DoubleHalfEdgePolytopeType, class QuadHalfEdgePolytopeType>
void convertPolytopeFromDoubleToQuadPrecision(
    DoubleSegmentedHalfEdgePolytopeType* a_polytope,
    DoubleHalfEdgePolytopeType* a_complete_polytope,
    QuadHalfEdgePolytopeType* a_converted_polytope) {
  // Convert polytope type
  using face_type = typename DoubleSegmentedHalfEdgePolytopeType::face_type;
  using half_edge_type = typename face_type::half_edge_type;
  using vertex_type = typename half_edge_type::vertex_type;
  using pt_type = typename vertex_type::pt_type;
  using scalar_type = typename pt_type::value_type;
  using converted_scalar_type = convert_to_quad<scalar_type>;
  using converted_pt_type = PtBase<converted_scalar_type>;
  using converted_vertex_type = VertexQuadratic<converted_pt_type>;
  using converted_halfedge_type = HalfEdgeQuadratic<converted_vertex_type>;
  using converted_face_type = FaceQuadratic<converted_halfedge_type>;
  const UnsignedIndex_t converted_kMaxHalfEdges =
      DoubleHalfEdgePolytopeType::maxHalfEdges;
  const UnsignedIndex_t converted_kMaxVertices =
      DoubleHalfEdgePolytopeType::maxVertices;
  const UnsignedIndex_t converted_kMaxFaces =
      DoubleHalfEdgePolytopeType::maxFaces;

  assert(a_polytope->checkValidHalfEdgeStructure());

  // Calculate number of half edges, vertices and faces
  const UnsignedIndex_t number_of_vertices = a_polytope->getNumberOfVertices();
  const UnsignedIndex_t number_of_faces = a_polytope->getNumberOfFaces();

#if 1
  UnsignedIndex_t number_of_half_edges = 0;
  std::vector<face_type*> face_mapping;
  face_mapping.resize(number_of_faces);
  std::vector<half_edge_type*> hald_edge_mapping;
  hald_edge_mapping.resize(0);
  std::vector<vertex_type*> vertex_mapping;
  vertex_mapping.resize(number_of_vertices);

  UnsignedIndex_t f = 0;
  for (const auto& face : (*a_polytope)) {
    face_mapping[f++] = face;
    const auto starting_half_edge = face->getStartingHalfEdge();
    auto current_half_edge = starting_half_edge;
    do {
      hald_edge_mapping.push_back(current_half_edge);
      number_of_half_edges++;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != starting_half_edge);
  }
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    vertex_mapping[v] = a_polytope->getVertex(v);
  }

  a_converted_polytope->resize(number_of_half_edges, number_of_vertices,
                               number_of_faces);

  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    const auto old_pt = a_polytope->getVertex(v)->getLocation();
    auto& new_pt = a_converted_polytope->getVertex(v).getLocation();
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      new_pt[d] = converted_scalar_type(old_pt[d]);
    }
  }

  UnsignedIndex_t e = 0;
  for (const auto& half_edge_DP : hald_edge_mapping) {
    // Pointer to QP halfedge
    converted_halfedge_type* half_edge_QP =
        &a_converted_polytope->getHalfEdge(e++);

    // Find QP face on which it lies
    const auto face_DP = half_edge_DP->getFace();
    const UnsignedIndex_t index_face = positionInMapping(face_mapping, face_DP);
    converted_face_type* face_QP = &a_converted_polytope->getFace(index_face);

    // Find end-point vertex
    const auto vertex_DP = half_edge_DP->getVertex();
    const UnsignedIndex_t index_vertex =
        positionInMapping(vertex_mapping, vertex_DP);
    converted_vertex_type* vertex_QP =
        &a_converted_polytope->getVertex(index_vertex);

    // Find previous, next and opposite QP halfedges
    const auto previous_DP = half_edge_DP->getPreviousHalfEdge();
    const UnsignedIndex_t index_previous =
        positionInMapping(hald_edge_mapping, previous_DP);
    converted_halfedge_type* previous_QP =
        &a_converted_polytope->getHalfEdge(index_previous);
    const auto next_DP = half_edge_DP->getNextHalfEdge();
    const UnsignedIndex_t index_next =
        positionInMapping(hald_edge_mapping, next_DP);
    converted_halfedge_type* next_QP =
        &a_converted_polytope->getHalfEdge(index_next);
    const auto opposite_DP = half_edge_DP->getOppositeHalfEdge();
    const UnsignedIndex_t index_opposite =
        positionInMapping(hald_edge_mapping, opposite_DP);
    converted_halfedge_type* opposite_QP =
        &a_converted_polytope->getHalfEdge(index_opposite);

    *half_edge_QP =
        converted_halfedge_type(vertex_QP, previous_QP, next_QP, face_QP);
    vertex_QP->setHalfEdge(half_edge_QP);
    half_edge_QP->setOppositeHalfEdge(opposite_QP);
    face_QP->setStartingHalfEdge(half_edge_QP);
  }
#else
  // Old implementation: can lead to invalid structures if faces have less
  // than 3 half-edges

  // Convert vertices to QP
  std::vector<PtBase<Quad_t>> pt_list;
  pt_list.resize(number_of_vertices);
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    const auto old_pt = a_polytope->getVertex(v)->getLocation();
    pt_list[v] = PtBase<Quad_t>(static_cast<Quad_t>(old_pt[0]),
                                static_cast<Quad_t>(old_pt[1]),
                                static_cast<Quad_t>(old_pt[2]));
  }

  // Create face -> halfedge mapping
  std::vector<std::vector<UnsignedIndex_t>> face_mapping;
  face_mapping.resize(number_of_faces);
  UnsignedIndex_t f = 0;
  for (const auto& face : (*a_polytope)) {
    const auto starting_half_edge = face->getStartingHalfEdge();
    auto current_half_edge = starting_half_edge;
    do {
      bool found = false;
      for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
        if (current_half_edge->getVertex() == a_polytope->getVertex(v)) {
          face_mapping[f].push_back(v);
          found = true;
          break;
        }
      }
      assert(found);
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != starting_half_edge);
    assert(face_mapping[f].size() >= 3);
    f++;
  }
  assert(f == number_of_faces);

  // Create converted QP polytope
  BREPToHalfEdge<converted_pt_type, converted_vertex_type,
                 converted_halfedge_type, converted_face_type,
                 converted_kMaxHalfEdges, converted_kMaxVertices,
                 converted_kMaxFaces>::setHalfEdgeVersion(pt_list, face_mapping,
                                                          a_converted_polytope);
#endif
}

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType, class ScalarType>
void triangulatePolytopeAndComputeNormals(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ScalarType a_nudge_epsilon, bool* a_requires_nudge) {
  using VertexType = typename HalfEdgePolytopeType::vertex_type;
  using PtType = typename VertexType::pt_type;
  using HalfEdgeType = typename HalfEdgePolytopeType::half_edge_type;
  using FaceType = typename HalfEdgePolytopeType::face_type;
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;

  // Loop over all faces. Completely face independent procedure
  const ScalarType EPSILON_NORMAL_SQ =
      machine_epsilon<ScalarType>() * machine_epsilon<ScalarType>();
  const ScalarType EPSILON_NORMAL_DIFF_SQ = ScalarType(100) * EPSILON_NORMAL_SQ;
  const auto nfaces = a_polytope->getNumberOfFaces();
  UnsignedIndex_t new_faces = 0;
  for (UnsignedIndex_t f = 0; f < nfaces; ++f) {
    auto face = (*a_polytope)[f];
    auto starting_half_edge = face->getStartingHalfEdge();
    const auto& start_location =
        starting_half_edge->getVertex()->getLocation().getPt();
    auto half_edge = starting_half_edge->getNextHalfEdge();
    auto next = half_edge->getNextHalfEdge();

    // Compute normal from starting half edge
    Normal normal = crossProduct(
        half_edge->getVertex()->getLocation().getPt() - start_location,
        next->getVertex()->getLocation().getPt() - start_location);
    ScalarType squaredMag = squaredMagnitude(normal);
    if (squaredMag < EPSILON_NORMAL_SQ) {
      normal = Normal(0, 0, 0);
    } else {
      normal /= sqrt(squaredMag);
    }

    // If face is already a triangle, move to next face
    if (half_edge == starting_half_edge || next == starting_half_edge ||
        next->getNextHalfEdge() == starting_half_edge) {
      face->setPlane(Plane(normal, normal * start_location));
      face->setAsTriangle();
      continue;
    }

    // Check for non-planarity of the face (only in DP case)
    if constexpr (std::is_same_v<ScalarType, double>) {
      bool need_triangulation = false;
      half_edge = next;
      next = next->getNextHalfEdge();
      do {
        Normal new_normal = crossProduct(
            half_edge->getVertex()->getLocation().getPt() - start_location,
            next->getVertex()->getLocation().getPt() - start_location);
        new_normal.normalize();
        ScalarType normal_diff_sq = squaredMagnitude(normal - new_normal);
        if (normal_diff_sq > EPSILON_NORMAL_DIFF_SQ) {
          need_triangulation = true;
          break;
        }
        half_edge = next;
        next = next->getNextHalfEdge();
      } while (next != starting_half_edge);
      if (!need_triangulation) {
        face->setPlane(Plane(normal, normal * start_location));
        face->setAsNotTriangle();
        continue;
      }
    }

// Triangulate face
#if 0
    // Old implementation that does not introduce a Steiner vertex
    // In some peculiar cases (associated with localizer cuts) this
    // can lead to invalid half-edge structures and therefore wrong
    // moments and/or segfaults
    face->setAsTriangle();
    half_edge = starting_half_edge->getNextHalfEdge();
    next = half_edge->getNextHalfEdge();
    auto next_next = next->getNextHalfEdge();
    while (next_next != starting_half_edge) {
      // Create new triangular face
      a_polytope->addFace(a_complete_polytope->getNewFace(FaceType(half_edge)));
      auto new_face = (*a_polytope)[nfaces + new_faces++];
      new_face->setPlane(Plane(normal, normal * start_location));
      new_face->setAsTriangle();

      // Creating new half-edge to close new triangular face
      auto new_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
          starting_half_edge->getVertex(), next, half_edge, new_face));
      half_edge->setPreviousHalfEdge(new_half_edge);
      half_edge->setFace(new_face);
      next->setNextHalfEdge(new_half_edge);
      next->setFace(new_face);

      // Creating new opposite half-edge
      auto new_opposite_half_edge = a_complete_polytope->getNewHalfEdge(
          HalfEdgeType(next->getVertex(), starting_half_edge, next_next, face));
      setMutualOpposites(new_half_edge, new_opposite_half_edge);
      starting_half_edge->setNextHalfEdge(new_opposite_half_edge);
      next_next->setPreviousHalfEdge(new_opposite_half_edge);

      // Closing old face
      half_edge = new_opposite_half_edge;
      next = next_next;
      next_next = next->getNextHalfEdge();

      // Compute normal of next triangle
      normal = crossProduct(
          half_edge->getVertex()->getLocation().getPt() - start_location,
          next->getVertex()->getLocation().getPt() - start_location);
      squaredMag = squaredMagnitude(normal);
      if (squaredMag < EPSILON_NORMAL_SQ) {
        normal = Normal(0, 0, 0);
      } else {
        normal /= sqrt(squaredMag);
      }

      assert(starting_half_edge->getVertex()->checkValidHalfEdgeCycle());
      assert(new_half_edge->getVertex()->checkValidHalfEdgeCycle());
    }

    face->setPlane(Plane(normal, normal * start_location));
#else
    // New implementation introducing a Steiner point (the approximate face
    // centroid) This is more expensive but prevent creating half-edge
    // duplicates
    Pt barycenter(0, 0, 0);
    half_edge = starting_half_edge;
    std::size_t npts = 0;
    do {
      barycenter += half_edge->getVertex()->getLocation().getPt();
      ++npts;
      half_edge = half_edge->getNextHalfEdge();
    } while (half_edge != starting_half_edge);
    assert(npts >= 3);
    barycenter *= (ScalarType(1) / ScalarType(npts));
    PtType face_center(barycenter);

    // Add new vertex to the polytope
    auto center_vert =
        a_complete_polytope->getNewVertex(VertexType(face_center));
    center_vert->setToSeek();
    center_vert->markToBeNotClipped();
    a_polytope->addVertex(center_vert);

    // Loop over face and triangulate now
    std::array<HalfEdgeType*, 2> new_he{{nullptr, nullptr}};
    // First triangle
    const auto previous_half_edge = half_edge->getPreviousHalfEdge();
    new_he[0] = a_complete_polytope->getNewHalfEdge(
        HalfEdgeType(center_vert, half_edge, nullptr, face));
    center_vert->setHalfEdge(new_he[0]);
    new_he[1] = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
        half_edge->getPreviousVertex(), new_he[0], half_edge, face));
    new_he[0]->setNextHalfEdge(new_he[1]);
    half_edge->setNextHalfEdge(new_he[0]);
    half_edge->setPreviousHalfEdge(new_he[1]);

    // Compute normal of first triangle
    normal = crossProduct(new_he[0]->getVertex()->getLocation().getPt() -
                              half_edge->getVertex()->getLocation().getPt(),
                          new_he[1]->getVertex()->getLocation().getPt() -
                              new_he[0]->getVertex()->getLocation().getPt());
    squaredMag = squaredMagnitude(normal);
    if (squaredMag < EPSILON_NORMAL_SQ) {
      normal = Normal(0, 0, 0);
    } else {
      normal /= sqrt(squaredMag);
    }
    face->setPlane(Plane(normal, normal * barycenter));
    face->setAsTriangle();

    // Move to previous half-edge and loop in reversed order
    half_edge = previous_half_edge;
    while (half_edge != starting_half_edge) {
      const auto following_half_edge = half_edge->getPreviousHalfEdge();

      // Middle face
      auto new_face = a_complete_polytope->getNewFace(FaceType(half_edge));
      a_polytope->addFace(new_face);
      half_edge->setFace(new_face);
      new_he[0] = a_complete_polytope->getNewHalfEdge(
          HalfEdgeType(center_vert, half_edge, nullptr, new_face));
      setMutualOpposites(new_he[0], new_he[1]);
      new_he[1] = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
          half_edge->getPreviousVertex(), new_he[0], half_edge, new_face));
      new_he[0]->setNextHalfEdge(new_he[1]);
      half_edge->setNextHalfEdge(new_he[0]);
      half_edge->setPreviousHalfEdge(new_he[1]);

      // Compute normal of triangle
      normal = crossProduct(new_he[0]->getVertex()->getLocation().getPt() -
                                half_edge->getVertex()->getLocation().getPt(),
                            new_he[1]->getVertex()->getLocation().getPt() -
                                new_he[0]->getVertex()->getLocation().getPt());
      squaredMag = squaredMagnitude(normal);
      if (squaredMag < EPSILON_NORMAL_SQ) {
        normal = Normal(0, 0, 0);
      } else {
        normal /= sqrt(squaredMag);
      }
      new_face->setPlane(Plane(normal, normal * barycenter));
      new_face->setAsTriangle();

      // Traverse in reversed order
      half_edge = following_half_edge;
    }
    // Final link
    auto first_new_edge = center_vert->getHalfEdge();
    setMutualOpposites(new_he[1], first_new_edge);
#endif
  }

  assert(a_polytope->checkValidHalfEdgeStructure());
  assert(checkTriangleFlags(a_polytope));
}

template <class HalfEdgePolytopeType>
inline bool checkTriangleFlags(HalfEdgePolytopeType* a_polytope) {
  for (UnsignedIndex_t f = 0; f < a_polytope->getNumberOfFaces(); ++f) {
    auto face = (*a_polytope)[f];
    const auto starting_half_edge = face->getStartingHalfEdge();
    auto three_after = starting_half_edge->getNextHalfEdge();
    if (face->isTriangle() && three_after == starting_half_edge) {
      continue;
    }
    if (!face->isTriangle() && three_after == starting_half_edge) {
      return false;
    }
    three_after = three_after->getNextHalfEdge();
    if (face->isTriangle() && three_after == starting_half_edge) {
      continue;
    }
    if (!face->isTriangle() && three_after == starting_half_edge) {
      return false;
    }
    three_after = three_after->getNextHalfEdge();
    if (face->isTriangle() && three_after != starting_half_edge) {
      return false;
    }
    if (!face->isTriangle() && three_after == starting_half_edge) {
      return false;
    }
  }
  return true;
}

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_QUADRATIC_INTERSECTION_QUADRATIC_INTERSECTION_TPP_