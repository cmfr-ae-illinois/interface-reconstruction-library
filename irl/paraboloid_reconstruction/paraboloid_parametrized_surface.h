// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_PARAMETRIZED_SURFACE_H_
#define IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_PARAMETRIZED_SURFACE_H_

#include <vector>
#include "irl/quadratic_reconstruction/parametrized_surface.h"

// #define IRL_USE_EARCUT
// #define IRL_USE_TRIANGLE
#define IRL_USE_CGAL
// #define IRL_USE_GEOGRAM

#ifdef IRL_USE_EARCUT
#include "external/earcut.hpp/include/mapbox/earcut.hpp"
#elif defined IRL_USE_TRIANGLE
#include "external/triangle/triangle.h"
#elif defined IRL_USE_CGAL
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#elif defined IRL_USE_GEOGRAM
#include "external/geogram.psm.Delaunay/Delaunay_psm.h"
#endif

#include "irl/geometry/general/normal.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/quadratic_reconstruction/ellipse.h"
#include "irl/quadratic_reconstruction/rational_bezier_arc.h"
#include "irl/surface_mesher/triangulated_surface.h"

namespace IRL {

template <class C>
struct has_paraboloid_surface : std::false_type {};

template <class C>
struct has_paraboloid_surface<const C> : has_paraboloid_surface<C> {};

template <class MomentType>
struct has_paraboloid_surface<AddSurfaceOutput<MomentType, ParaboloidParametrizedSurfaceOutput>>
    : std::true_type {};

/// \brief Parametrized surface defined by coeffs A,B of paraboloid + list of
/// rational Bézier arcs
class ParaboloidParametrizedSurfaceOutput : public ParametrizedSurfaceOutput {
 public:
  /// \brief Default constructor.
  ParaboloidParametrizedSurfaceOutput(void);
  ParaboloidParametrizedSurfaceOutput(const Paraboloid& a_paraboloid);

  ParaboloidParametrizedSurfaceOutput(
      const ParaboloidParametrizedSurfaceOutput& a_rhs);
  ParaboloidParametrizedSurfaceOutput(
      ParaboloidParametrizedSurfaceOutput&& a_rhs);

  ParaboloidParametrizedSurfaceOutput& operator=(
      const ParaboloidParametrizedSurfaceOutput& a_rhs);
  ParaboloidParametrizedSurfaceOutput& operator=(
      ParaboloidParametrizedSurfaceOutput&& a_rhs);

  void setParaboloid(const Paraboloid& a_paraboloid);

  const Paraboloid& getParaboloid(void) const;
  inline double getSurfaceArea(void);
  inline double getMeanCurvatureIntegral(void);
  inline double getGaussianCurvatureIntegral(void);
  inline Normal getAverageNormal(void);
  inline Normal getAverageNormalNonAligned(void);
  inline Normal getNormalAligned(const Pt a_pt);
  inline Normal getNormalNonAligned(const Pt a_pt);
  inline double getMeanCurvatureAligned(const Pt a_pt);
  inline double getMeanCurvatureNonAligned(const Pt a_pt);
  inline double getGaussianCurvatureAligned(const Pt a_pt);
  inline double getGaussianCurvatureNonAligned(const Pt a_pt);

  void triangulate_fromPtr(
      const double a_length_scale = -1.0, const UnsignedIndex_t a_nsplit = 5,
      TriangulatedSurfaceOutput* a_surface = nullptr) const;

  TriangulatedSurfaceOutput triangulate(
      const double a_length_scale = -1.0,
      const UnsignedIndex_t a_nsplit = 5) const;

  ~ParaboloidParametrizedSurfaceOutput(void);

 private:
  Paraboloid paraboloid_m;
};

inline std::ostream& operator<<(std::ostream& out,
                                const ParaboloidParametrizedSurfaceOutput&
                                    a_paraboloid_parametrized_surface);

}  // namespace IRL

#include "irl/paraboloid_reconstruction/paraboloid_parametrized_surface.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_PARAMETRIZED_SURFACE_H_
