// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "gtest/gtest.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace {

using namespace IRL;

TEST(ParaboloidPolygonIntersection, AreaAndArcLengthPolygonParabola) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;

  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, 0.5, 0.5);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(0, 0, -1);
  polygon[1] = Pt(1, 0, -1);
  polygon[2] = Pt(1, 0, 1);
  polygon[3] = Pt(0, 0, 1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = 5.0 / 6.0;
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that routine returns a single arc for the intersection
    EXPECT_EQ(arcs.size(), 1);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    const auto& arc = arcs[0];
    double arc_length = 0.0;
    for (int k = 0; k < QuadRuleOrder; k++) {
      const double t = 0.5 * (1.0 + abscissea[k]);
      const double w = 0.5 * weights[k];
      const Pt Cp = arc.derivative(t);
      const double Jac =
          std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
      arc_length += w * Jac;
    }

    // Exact expected value
    const double expected_arc_length =
        std::sqrt(2.0) / 2.0 + std::log(1.0 + std::sqrt(2.0)) / 2.0;

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }

  // Now reverse order of polygon and do the same
  {
    polygon[0] = Pt(0, 0, -1);
    polygon[3] = Pt(1, 0, -1);
    polygon[2] = Pt(1, 0, 1);
    polygon[1] = Pt(0, 0, 1);
    polygon.calculateAndSetPlaneOfExistence();
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = 5.0 / 6.0;
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that routine returns a single arc for the intersection
    EXPECT_EQ(arcs.size(), 1);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    const auto& arc = arcs[0];
    double arc_length = 0.0;
    for (int k = 0; k < QuadRuleOrder; k++) {
      const double t = 0.5 * (1.0 + abscissea[k]);
      const double w = 0.5 * weights[k];
      const Pt Cp = arc.derivative(t);
      const double Jac =
          std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
      arc_length += w * Jac;
    }

    // Exact expected value
    const double expected_arc_length =
        std::sqrt(2.0) / 2.0 + std::log(1.0 + std::sqrt(2.0)) / 2.0;

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }
}

TEST(ParaboloidPolygonIntersection, AreaAndArcLengthTriangleParabola) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;

  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));

  const auto paraboloid = Paraboloid(datum, frame, 0.5, 0.5);

  Tri triangle({Pt(0, 0, -1), Pt(2, 0, 0), Pt(0, 0, 1)});
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(triangle, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    EXPECT_EQ(arcs.size(), 1);

    // Verify that the area of the intersection is as expected
    const double expected_area = 7.0 / 12.0;
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    const auto& arc = arcs[0];
    double arc_length = 0.0;
    for (int k = 0; k < QuadRuleOrder; k++) {
      const double t = 0.5 * (1.0 + abscissea[k]);
      const double w = 0.5 * weights[k];
      const Pt Cp = arc.derivative(t);
      const double Jac =
          std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
      arc_length += w * Jac;
    }

    // Exact expected value
    const double expected_arc_length =
        std::sqrt(2.0) / 2.0 + std::log(1.0 + std::sqrt(2.0)) / 2.0;

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }

  // Now reverse order of polygon and do the same
  {
    triangle = Tri({Pt(0, 0, -1), Pt(0, 0, 1), Pt(2, 0, 0)});
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(triangle, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    EXPECT_EQ(arcs.size(), 1);

    // Verify that the area of the intersection is as expected
    const double expected_area = 7.0 / 12.0;
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    const auto& arc = arcs[0];
    double arc_length = 0.0;
    for (int k = 0; k < QuadRuleOrder; k++) {
      const double t = 0.5 * (1.0 + abscissea[k]);
      const double w = 0.5 * weights[k];
      const Pt Cp = arc.derivative(t);
      const double Jac =
          std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
      arc_length += w * Jac;
    }

    // Exact expected value
    const double expected_arc_length =
        std::sqrt(2.0) / 2.0 + std::log(1.0 + std::sqrt(2.0)) / 2.0;

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }
}

TEST(ParaboloidPolygonIntersection, AreaAndArcLengthWithNudgePolygonParabola) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;

  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, 1.0, 1.0);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(-1, 0, -1);
  polygon[1] = Pt(1, 0, -1);
  polygon[2] = Pt(1, 0, 1);
  polygon[3] = Pt(-1, 0, 1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = 4.0 / 3.0;
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length =
        std::sqrt(5.0) + std::log(2.0 + std::sqrt(5.0)) / 2.0;

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }

  // Now reverse order of polygon and do the same
  {
    polygon[0] = Pt(-1, 0, -1);
    polygon[3] = Pt(1, 0, -1);
    polygon[2] = Pt(1, 0, 1);
    polygon[1] = Pt(-1, 0, 1);
    polygon.calculateAndSetPlaneOfExistence();
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = 4.0 / 3.0;
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length =
        std::sqrt(5.0) + std::log(2.0 + std::sqrt(5.0)) / 2.0;

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }
}

TEST(ParaboloidPolygonIntersection, AreaAndArcLengthPolygonCircleArc) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;

  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, 1.0, 1.0);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(0, 0, -1);
  polygon[1] = Pt(2, 0, -1);
  polygon[2] = Pt(2, 2, -1);
  polygon[3] = Pt(0, 2, -1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = machine_pi<double>() / 4.0;
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length = machine_pi<double>() / 2.0;

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }

  // Now reverse order of polygon and do the same
  polygon[0] = Pt(0, 0, -1);
  polygon[3] = Pt(2, 0, -1);
  polygon[2] = Pt(2, 2, -1);
  polygon[1] = Pt(0, 2, -1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = machine_pi<double>() / 4.0;
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length = machine_pi<double>() / 2.0;

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }
}

TEST(ParaboloidPolygonIntersection, AreaPolygonEllipseArc) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;

  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, 1.0, 2.0);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(0, 0, -1);
  polygon[1] = Pt(2, 0, -1);
  polygon[2] = Pt(2, 2, -1);
  polygon[3] = Pt(0, 2, -1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = machine_pi<double>() / 4.0 / std::sqrt(2.0);
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Note: there is no closed form for the arc length of an ellipse.
  }

  // Now reverse order of polygon and do the same
  polygon[0] = Pt(0, 0, -1);
  polygon[3] = Pt(2, 0, -1);
  polygon[2] = Pt(2, 2, -1);
  polygon[1] = Pt(0, 2, -1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = machine_pi<double>() / 4.0 / std::sqrt(2.0);
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);
  }
}

TEST(ParaboloidPolygonIntersection, AreaPolygonEllipseArcReverse) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;

  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, -1.0, -2.0);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(0, 0, 1);
  polygon[1] = Pt(2, 0, 1);
  polygon[2] = Pt(2, 2, 1);
  polygon[3] = Pt(0, 2, 1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area =
        4.0 - machine_pi<double>() / 4.0 / std::sqrt(2.0);
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Note: there is no closed form for the arc length of an ellipse.
  }

  // Now reverse order of polygon and do the same
  polygon[0] = Pt(0, 0, 1);
  polygon[3] = Pt(2, 0, 1);
  polygon[2] = Pt(2, 2, 1);
  polygon[1] = Pt(0, 2, 1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area =
        4.0 - machine_pi<double>() / 4.0 / std::sqrt(2.0);
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);
  }
}

TEST(ParaboloidPolygonIntersection, AreaAndArcLengthCircleInscribedInSquare) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;

  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, 1.0, 1.0);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(-1, -1, -1);
  polygon[1] = Pt(1, -1, -1);
  polygon[2] = Pt(1, 1, -1);
  polygon[3] = Pt(-1, 1, -1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = machine_pi<double>();
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length = 2.0 * machine_pi<double>();

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-7);
  }

  // Now reverse order of polygon and do the same
  polygon[0] = Pt(-1, -1, -1);
  polygon[3] = Pt(1, -1, -1);
  polygon[2] = Pt(1, 1, -1);
  polygon[1] = Pt(-1, 1, -1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = machine_pi<double>();
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length = 2.0 * machine_pi<double>();

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-7);
  }
}

TEST(ParaboloidPolygonIntersection, AreaAndArcLengthCircleFaceOnly) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;

  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, 1.0, 1.0);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(-2, -2, -1);
  polygon[1] = Pt(2, -2, -1);
  polygon[2] = Pt(2, 2, -1);
  polygon[3] = Pt(-2, 2, -1);
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = machine_pi<double>();
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length = 2.0 * machine_pi<double>();

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }

  // Now reverse order of polygon and do the same
  polygon[0] = Pt(-2, -2, -1);
  polygon[3] = Pt(2, -2, -1);
  polygon[2] = Pt(2, 2, -1);
  polygon[1] = Pt(-2, 2, -1);
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = machine_pi<double>();
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length = 2.0 * machine_pi<double>();

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }
}

TEST(ParaboloidPolygonIntersection, AreaAndArcLengthCircleFaceOnlyReverse) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;

  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, -1.0, -1.0);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(-2, -2, 1);
  polygon[1] = Pt(2, -2, 1);
  polygon[2] = Pt(2, 2, 1);
  polygon[3] = Pt(-2, 2, 1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = 16.0 - machine_pi<double>();
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length = 2.0 * machine_pi<double>();

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }

  // Now reverse order of polygon and do the same
  polygon[0] = Pt(-2, -2, 1);
  polygon[3] = Pt(2, -2, 1);
  polygon[2] = Pt(2, 2, 1);
  polygon[1] = Pt(-2, 2, 1);
  polygon.calculateAndSetPlaneOfExistence();
  {
    auto area_and_surf =
        intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid);
    const auto arcs = area_and_surf.getSurface().getArcs();

    // Verify that the area of the intersection is as expected
    const double expected_area = 16.0 - machine_pi<double>();
    EXPECT_NEAR(area_and_surf.getMoments().area(), expected_area, 1.0e-14);

    // Verify that the arc length of the intersection is as expected
    constexpr int QuadRuleOrder = 100;
    const auto& abscissea = AbscissaeGauss<double, QuadRuleOrder>();
    const auto& weights = WeightsGauss<double, QuadRuleOrder>();
    double arc_length = 0.0;
    for (const auto& arc : arcs) {
      for (int k = 0; k < QuadRuleOrder; k++) {
        const double t = 0.5 * (1.0 + abscissea[k]);
        const double w = 0.5 * weights[k];
        const Pt Cp = arc.derivative(t);
        const double Jac =
            std::sqrt(Cp[0] * Cp[0] + Cp[1] * Cp[1] + Cp[2] * Cp[2]);
        arc_length += w * Jac;
      }
    }

    // Exact expected value
    const double expected_arc_length = 2.0 * machine_pi<double>();

    EXPECT_NEAR(arc_length, expected_arc_length, 1.0e-14);
  }
}

}  // namespace
