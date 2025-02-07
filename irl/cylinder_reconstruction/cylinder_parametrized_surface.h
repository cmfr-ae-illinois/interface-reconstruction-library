// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_CYLINDER_RECONSTRUCTION_CYLINDER_PARAMETRIZED_SURFACE_H_
#define IRL_CYLINDER_RECONSTRUCTION_CYLINDER_PARAMETRIZED_SURFACE_H_

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

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/geometry/general/normal.h"
#include "irl/quadratic_reconstruction/ellipse.h"
#include "irl/quadratic_reconstruction/rational_bezier_arc.h"
#include "irl/surface_mesher/triangulated_surface.h"

namespace IRL {

template <class C>
struct has_cylinder_surface : std::false_type {};

template <class C>
struct has_cylinder_surface<const C> : has_cylinder_surface<C> {};

template <class MomentType>
struct has_cylinder_surface<AddSurfaceOutput<MomentType, CylinderParametrizedSurfaceOutput>>
    : std::true_type {};

/// \brief Parametrized surface defined by coeffs B and R of Cylinder + list of
/// rational Bézier arcs
class CylinderParametrizedSurfaceOutput : public ParametrizedSurfaceOutput {
 public:
  /// \brief Default constructor.
  CylinderParametrizedSurfaceOutput(void);
  CylinderParametrizedSurfaceOutput(const Cylinder& a_cylinder);

  CylinderParametrizedSurfaceOutput(
      const CylinderParametrizedSurfaceOutput& a_rhs);
  CylinderParametrizedSurfaceOutput(CylinderParametrizedSurfaceOutput&& a_rhs);

  CylinderParametrizedSurfaceOutput& operator=(
      const CylinderParametrizedSurfaceOutput& a_rhs);
  CylinderParametrizedSurfaceOutput& operator=(
      CylinderParametrizedSurfaceOutput&& a_rhs);

  void setCylinder(const Cylinder& a_cylinder);
  void resetCylinder(void);
  void setScale(double a_scale);

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

  // ~CylinderParametrizedSurfaceOutput(void);

 private:
  std::vector<int> indexes_of_flip;
  std::vector<Cylinder> cylinder_m;
  double scale_m;
};

inline std::ostream& operator<<(
    std::ostream& out,
    const CylinderParametrizedSurfaceOutput& a_cylinder_parametrized_surface);

}  // namespace IRL

#include "irl/cylinder_reconstruction/cylinder_parametrized_surface.tpp"

#endif  // IRL_CYLINDER_RECONSTRUCTION_CYLINDER_PARAMETRIZED_SURFACE_H_
