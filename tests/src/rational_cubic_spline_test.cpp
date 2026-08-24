#include "gtest/gtest.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/geometry/spline/rational_cubic_bezier_arc.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/quadratic_reconstruction/rational_bezier_arc.h"

namespace {
using namespace IRL;

TEST(RationalCubicSpline, PointEvaluation) {
  Pt start = Pt(0.0, 0.0, 0.0);
  Pt control1 = Pt(0.0, 1.0, 1.0);
  Pt control2 = Pt(1.0, 1.0, 0.0);
  Pt end = Pt(1.0, 0.0, 1.0);
  double weight1 = 1.0;
  double weight2 = 1.0;
  // Evaluation Points
  double t1 = 0.3;
  double t2 = 0.5;
  double t3 = 0.0;
  double t4 = 1.0;
  // Expected Results
  Pt expected1 = Pt(0.216, 0.63, 0.468);
  Pt expected2 = Pt(0.5, 0.75, 0.5);
  Pt expected3 = Pt(0.0, 0.0, 0.0);
  Pt expected4 = Pt(1.0, 0.0, 1.0);
  // Create Nonrational Cubic Bezier Arc
  RationalCubicBezierArc arc(start, control1, control2, end, weight1, weight2);
  // evaluate at each point and compare to expected results
  Pt result1 = arc.point(t1);
  Pt result2 = arc.point(t2);
  Pt result3 = arc.point(t3);
  Pt result4 = arc.point(t4);
  // Compare results
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result1[i], expected1[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Point Evaluation Fail for Point 1 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result2[i], expected2[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Point Evaluation Fail for Point 2 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result3[i], expected3[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Point Evaluation Fail for Point 3 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result4[i], expected4[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Point Evaluation Fail for Point 4 index" << i;
  }
  // Rational Case
  weight1 = 2.0;
  weight2 = 0.5;
  RationalCubicBezierArc rational_arc(start, control1, control2, end, weight1,
                                      weight2);
  // Expected Results
  expected1 =
      Pt(0.090233939844040098, 0.725213516524322331, 0.675083549944300088);
  expected2 =
      Pt(0.263157894736842078, 0.789473684210526291, 0.736842105263157864);
  expected3 = Pt(0.0, 0.0, 0.0);
  expected4 = Pt(1.0, 0.0, 1.0);
  // Evaluate
  result1 = rational_arc.point(t1);
  result2 = rational_arc.point(t2);
  result3 = rational_arc.point(t3);
  result4 = rational_arc.point(t4);
  // Compare results
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result1[i], expected1[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Weighted Point Evaluation Fail for Point 1 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result2[i], expected2[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Weighted Point Evaluation Fail for Point 2 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result3[i], expected3[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Weighted Point Evaluation Fail for Point 3 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result4[i], expected4[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Weighted Point Evaluation Fail for Point 4 index" << i;
  }
  SUCCEED();
}

TEST(RationalCubicSpline, DerivativeEvaluation) {
  Pt start = Pt(0.0, 0.0, 0.0);
  Pt control1 = Pt(0.0, 1.0, 1.0);
  Pt control2 = Pt(1.0, 1.0, 0.0);
  Pt end = Pt(1.0, 0.0, 1.0);
  double weight1 = 1.0;
  double weight2 = 1.0;
  // Evaluation Points
  double t1 = 0.3;
  double t2 = 0.5;
  double t3 = 0.0;
  double t4 = 1.0;
  // Expected Results
  Pt expected1 = Pt(1.26, 1.2, 0.48);
  Pt expected2 = Pt(1.5, 0.0, 0.0);
  Pt expected3 = Pt(0.0, 3.0, 3.0);
  Pt expected4 = Pt(0.0, -3.0, 3.0);
  // Create Nonrational Cubic Bezier Arc
  RationalCubicBezierArc arc(start, control1, control2, end, weight1, weight2);
  // evaluate at each point and compare to expected results
  Pt result1 = arc.derivative(t1);
  Pt result2 = arc.derivative(t2);
  Pt result3 = arc.derivative(t3);
  Pt result4 = arc.derivative(t4);
  // Compare results
  double Scale = 1000.0;
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result1[i], expected1[i],
                Scale * std::numeric_limits<double>::epsilon())
        << "Derivative Evaluation Fail for Point 1 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result2[i], expected2[i],
                Scale * std::numeric_limits<double>::epsilon())
        << "Derivative Evaluation Fail for Point 2 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result3[i], expected3[i],
                Scale * std::numeric_limits<double>::epsilon())
        << "Derivative Evaluation Fail for Point 3 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result4[i], expected4[i],
                Scale * std::numeric_limits<double>::epsilon())
        << "Derivative Evaluation Fail for Point 4 index" << i;
  }
  // Rational Case
  weight1 = 2.0;
  weight2 = 0.5;
  RationalCubicBezierArc rational_arc(start, control1, control2, end, weight1,
                                      weight2);
  // Expected Results
  Pt expectedWeighted1 =
      Pt(0.58723852421504003, 0.83303813754877972, 0.655327747296045651);
  Pt expectedWeighted2 =
      Pt(1.196675900277008303, -0.199445983379501267, 0.066481994459833833);
  Pt expectedWeighted3 = Pt(0.0, 6.0, 6.0);
  Pt expectedWeighted4 = Pt(0.0, -1.5, 1.5);
  // Evaluate
  result1 = rational_arc.derivative(t1);
  result2 = rational_arc.derivative(t2);
  result3 = rational_arc.derivative(t3);
  result4 = rational_arc.derivative(t4);
  // Compare results

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result1[i], expectedWeighted1[i],
                Scale * std::numeric_limits<double>::epsilon())
        << "Weighted Derivative Evaluation Fail for Point 1 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result2[i], expectedWeighted2[i],
                Scale * std::numeric_limits<double>::epsilon())
        << "Weighted Derivative Evaluation Fail for Point 2 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result3[i], expectedWeighted3[i],
                Scale * std::numeric_limits<double>::epsilon())
        << "Weighted Derivative Evaluation Fail for Point 3 index" << i;
  }
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(result4[i], expectedWeighted4[i],
                Scale * std::numeric_limits<double>::epsilon())
        << "Weighted Derivative Evaluation Fail for Point 4 index" << i;
  }
  SUCCEED();
}

TEST(RationalCubicSpline, GettersTest) {
  Pt start = Pt(0.0, 0.0, 0.0);
  Pt control1 = Pt(0.0, 1.0, 1.0);
  Pt control2 = Pt(1.0, 1.0, 0.0);
  Pt end = Pt(1.0, 0.0, 1.0);
  double weight1 = 2.0;
  double weight2 = 0.5;
  RationalCubicBezierArc arc(start, control1, control2, end, weight1, weight2);
  // Test Getters
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(arc.start_point()[i], start[i],
                std::numeric_limits<double>::epsilon())
        << "Start Point Getter Fail for index " << i;
    EXPECT_NEAR(arc.control_point_1()[i], control1[i],
                std::numeric_limits<double>::epsilon())
        << "Control Point 1 Getter Fail for index " << i;
    EXPECT_NEAR(arc.control_point_2()[i], control2[i],
                std::numeric_limits<double>::epsilon())
        << "Control Point 2 Getter Fail for index " << i;
    EXPECT_NEAR(arc.end_point()[i], end[i],
                std::numeric_limits<double>::epsilon())
        << "End Point Getter Fail for index " << i;
  }
  EXPECT_NEAR(arc.weight(1), weight1, std::numeric_limits<double>::epsilon())
      << "Weight 1 Getter Fail";
  EXPECT_NEAR(arc.weight(2), weight2, std::numeric_limits<double>::epsilon())
      << "Weight 2 Getter Fail";
  SUCCEED();
}

TEST(RationalCubicSpline, EnergyMinimizationTest) {
  Pt start = Pt(0.0, 0.0, 0.0);
  Normal tangent1 = Normal(1.0, 0.0, 0.0);
  Normal tangent2 = Normal(0.6, 0.8, 0.0);
  Pt end = Pt(1.0, 0.0, 0.0);

  double expectedAlpha1 = 0.450549450549451;  // Double Check these values
  double expectedAlpha2 = 0.164835164835165;
  Pt expectedControl1 = start + expectedAlpha1 * tangent1;
  Pt expectedControl2 = end - expectedAlpha2 * tangent2;

  RationalCubicBezierArc arc(start, tangent1, end, tangent2);
  // check location of control points
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(arc.control_point_1()[i], expectedControl1[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Energy Minimization Fail for Control Point 1 index " << i;
    EXPECT_NEAR(arc.control_point_2()[i], expectedControl2[i],
                100 * std::numeric_limits<double>::epsilon())
        << "Energy Minimization Fail for Control Point 2 index " << i;
  }
  SUCCEED();
}

TEST(RationalCubicSpline, OutputToVTKTest) {
  Pt start = Pt(0.0, 0.0, 0.0);
  Pt control1 = Pt(0.0, 1.0, 1.0);
  Pt control2 = Pt(1.0, 1.0, 0.0);
  Pt end = Pt(1.0, 0.0, 1.0);

  RationalCubicBezierArc arc(start, control1, control2, end);
  std::string filename = "rational_cubic_bezier_arc";
  arc.saveToVTK(filename);
  // Check that file was created manually inspect file to ensure it is correct.
  SUCCEED();
}

TEST(RationalCubicSpline, BreakageTest) {
  Normal tangent1 = Normal(0.0, 1.0, 0.0);
  Normal tangent2 = Normal(-1.0, 0.0, 0.0);
  // Evaluate a point with alpha_0=-1 and alpha_1=-1 which should cause breakage
  Pt start = Pt(0.0, 0.0, 0.0);
  Pt end = Pt(2.0, -2.0, 0.0);
  RationalCubicBezierArc arc(start, tangent1, end, tangent2);
  std::string filename = "BreakingSpline";
  arc.saveToVTK(filename);
  SUCCEED();
}
// Levelset Function Tests =======================================

// =====  Circle ====
void circleLevelset(const double radius, const Pt& pt, double* levelset) {
  *levelset = pt[0] * pt[0] + pt[1] * pt[1] - radius * radius;
}
void circleLevelset(const double radius, const Pt& pt,
                    std::pair<double, Eigen::Vector3d>* retVal) {
  double F;
  Eigen::Vector3d gradF;
  // Implementation for circle levelset (Cylinder)
  F = pt[0] * pt[0] + pt[1] * pt[1] - radius * radius;
  gradF = Eigen::Vector3d(2.0 * pt[0], 2.0 * pt[1], 0.0);
  // Return
  *retVal = std::make_pair(F, gradF);
}
void circleLevelset(const double radius, const Pt& pt,
                    std::pair<double, Normal>* retVal) {
  double F;
  Normal tangent, gradF;
  // Implementation for circle levelset (Cylinder)
  F = pt[0] * pt[0] + pt[1] * pt[1] - radius * radius;
  gradF = Normal(2.0 * pt[0], 2.0 * pt[1], 0.0);
  tangent = Normal(-gradF[1], gradF[0], 0.0);
  // Return
  *retVal = std::make_pair(F, tangent);
}

TEST(RationalCubicSpline, CircleLevelsetTest) {
  double radius = 1.0;
  // Pick two points on circle - No need to solve bisection method for this
  // test.
  Pt pt1 = Pt(1.0, 0.0, 0.0);
  Pt pt2 = Pt(0.0, 1.0, 0.0);
  // Sample Value and Gradient
  std::pair<double, Normal> retVal1;
  std::pair<double, Normal> retVal2;
  circleLevelset(radius, pt1, &retVal1);
  circleLevelset(radius, pt2, &retVal2);
  // Check that value is near zero and gradient is correct
  double levelset1 = retVal1.first;
  double levelset2 = retVal2.first;
  Normal tangent1 = retVal1.second;
  Normal tangent2 = retVal2.second;
  tangent1.normalize();
  tangent2.normalize();
  std::cout << " Tangent 1: " << tangent1 << std::endl;
  std::cout << " Tangent 2: " << tangent2 << std::endl;
  EXPECT_NEAR(levelset1, 0.0, std::numeric_limits<double>::epsilon())
      << "Levelset Value Fail for Point 1";
  EXPECT_NEAR(levelset2, 0.0, std::numeric_limits<double>::epsilon())
      << "Levelset Value Fail for Point 2";
  tangent1.normalize();
  tangent2.normalize();

  RationalCubicBezierArc arc(pt1, tangent1, pt2, tangent2);
  std::cout << "THE OUTPUT ARC IS =====================" << std::endl;
  std::cout << arc;
  std::string filename = "Energy_Minimizing_RCBA_Circle_Levelset";
  arc.saveToVTK(filename);
  // Check that file was created manually inspect file to ensure it is correct.
  SUCCEED();
}

TEST(RationalCubicSpline, BasicallyEllipseTest) {
  // Pick two points on circle - No need to solve bisection method for this
  // test.
  Pt pt1 = Pt(1.0, 0.0, 0.0);
  Pt pt2 = Pt(0.0, 0.5, 0.0);
  Normal tangent1 = Normal(0.0, 1.0, 0.0);
  Normal tangent2 = Normal(-1.0, 0.0, 0.0);

  tangent1.normalize();
  tangent2.normalize();

  RationalCubicBezierArc arc(pt1, tangent1, pt2, tangent2);
  std::cout << "THE OUTPUT ARC IS =====================" << std::endl;
  std::cout << arc;
  std::string filename = "Energy_Minimizing_RCBA_Basically_Ellipse";
  arc.saveToVTK(filename);
  // Check that file was created manually inspect file to ensure it is correct.
  SUCCEED();
}

}  // namespace