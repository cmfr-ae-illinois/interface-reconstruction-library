#ifndef IRL_PARTITION_OF_UNITY_SURFACE_TENSION_TPP_
#define IRL_PARTITION_OF_UNITY_SURFACE_TENSION_TPP_

#include <limits>
#include <tuple>
#include <vector>

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/quadratic_intersection/moment_contributions.h"
#include "irl/geometry/spline/rational_cubic_bezier_arc.h"
#include "irl/helpers/wendland.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/moments/general_moments.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

namespace IRL {

// Constructors
template <class CellType>
PUST<CellType>::PUST(void) : PU<CellType>() {}

template <class CellType>
PUST<CellType>::PUST(const PUNeighborhood<CellType>& stencil_)
    : PU<CellType>(stencil_, 2.0) {}

template <class CellType>
PUST<CellType>::PUST(const PUNeighborhood<CellType>& stencil_,
                     const double kernel_size)
    : PU<CellType>(stencil_, kernel_size) {}

// Solve Edge
template <class CellType>
Normal PUST<CellType>::solveEdge(const double STin, const Pt& P0, const Pt& P1,
                                 const double delta, const double Pressure,
                                 const Normal& Marangoni) {
  // The Marangoni normal object holds the Xgradient, then the Y gradient, then
  // the temperature gradient of ST Marangoni = [Gx,Gy,sigma_T] Make Implicit
  // Surface std::cout << "In Solve Edge\n";
  double intersection_threshold_m = 1e-6;  // This is the value of the wendland
                                           // function at 0.5 (if radius is 1).
  // Something is up with paraboloids because they take like 8x the time to run.
  double STCoeff = STin;
  // The Pressure Option tells us if we should include the pressure terms or
  // not.
  // Find Intersection with edge
  bool blocked = false;
  std::vector<Pt> intersections =
      this->intersectEdge(P0, P1, 10, intersection_threshold_m, blocked);

  // Calculate some edge properties
  Pt dP = P1 - P0;
  Normal OutwardsNormal = {dP[1], -dP[0], 0};
  OutwardsNormal.normalize();
  double D = std::sqrt(dP[0] * dP[0] + dP[1] * dP[1]);
  double denom = 1.0 / (safelyEpsilon(D));
  // Give space for working variables
  Normal tangent;
  Normal gradient;
  Normal total = {0.0, 0.0, 0.0};  // Total force

  // Loop over intersections and add tangents to force
  if (intersections.size() > 0) {
    for (int j = 0; j < intersections.size(); j++) {
      Normal N = this->getNormal(intersections[j]);
      Normal Z = {0.0, 0.0, 1.0};
      tangent = IRL::crossProduct(Z, N);
      tangent.normalize();
      // Here we apply the Marangoni surface tension. To do this, we assume
      // STCoeff = STCoeff + gamma_T(T-T_0). Letting gamma_T=-0.002 (Ratio for
      // water). We also pick T-T0=Gx to be the form we use. This will give us
      // th surface tension at this point We pick G = 10 for the stake of
      // simplicity
      STCoeff = STin + Marangoni[2] * (Marangoni[0] * intersections[j][0] +
                                       Marangoni[1] * intersections[j][1]);
      if (Marangoni[2] == 0.0 &&
          (Marangoni[0] != 0.0 ||
           Marangoni[1] != 0.0)) {  // Force droplet breakup
        double gamma0 = 1.0;        // Base surface tension coefficient
        double R = 1.0;             // Characteristic length scale

        // Linear decrease in surface tension with x-coordinate - Al-Saud
        // Style
        STCoeff =
            gamma0 * std::max(1 - 1.25 * std::abs(intersections[j][0]), 0.1);
      }
      // If facing inside, multiply by negative
      // If facing outside, leave the same
      double Scale = OutwardsNormal * tangent;
      Scale = Scale / safelyEpsilon(std::abs(Scale));
      // Add
      total = total + Scale * STCoeff * denom * tangent;
      // Next, we need to calculate pressure contribution, if wanted
      // To determine if an intersection is going into or out of the fluid, we
      // can use the dot product of dP with the gradient (which is outwards
      // pointing). This product will contain orientation information also.
      if (Pressure >= 0.5) {
        gradient =
            Normal(tangent[1], -tangent[0],
                   0.0);  // Use this to get gradient so we don't have torecalc
        Normal dPN = {dP[0], dP[1], 0};  // Make dP into a vector
        dPN.normalize();
        Scale = gradient * dPN;                          // Dot Product
        Scale = Scale / safelyEpsilon(std::abs(Scale));  // Sign Only
        // This scale tells me if, going along the edge from P0 to P1, we are
        // going into the fluid (Scale >0) or out of the fluid (Scale <0).

        // Now we  have the sign. The domain that each intersection influences
        // is from the centerpoint out. This means that if the distance from P0
        // to the intersection if < 0.5D, then we look at the distance
        // between the intersection and P0. If it is >0.5D, then look look
        // at the distance between it and P1, which will be the same as
        // D-[the initial distance].

        Pt dPI = intersections[j] - P0;
        double L = std::sqrt(dPI[0] * dPI[0] +
                             dPI[1] * dPI[1]);  // Calculated Distance

        double Sx = L * denom;  // Also definitely positive
        if (L > 0.5 * D) {      // If More than 0.5D
          L = L - D;            // Switch to other direction (Center Out)
        }

        // Now that we have the length of section, we just need to calculat
        // the curvature, then calculat the pressure jump, then multiply it
        // all together and add

        double curv = this->getMeanCurvature(intersections[j]);
        // This also contains some direction information, but I am unsure how
        // this needs to be positive or not. I am going to leave it for now,
        // and we will see if that is a problem. If there seems to be issues,
        // try to make this stritly positive.

        // Now we just add this all together. Note, the pressure force
        // alwaysacts inwards, so we multiply by the negative of the outwards
        // pointingnormal
        if (std::abs(Scale) >=
            1e-10) {  // Ensure that in the undefined areas, we
          // // do not add pressure terms.
          total = total + Scale * STCoeff * curv * L * denom * OutwardsNormal;
          // std::cout << "Surface Tension Coeff: " << STCoeff << "\n"
          //           << "Curvature: " << curv << "\n"
          //           << "Length Segment: " << L << "\n"
          //           << "Pressure Contribution Scale: " << Scale << "\n"
          //           << "Intersection Point: " << intersections[j] << "\n"
          //           << "Point 1: " << P0 << "\n"
          //           << "Point 2: " << P1 << "\n"
          //           << "--------------------------------------\n";
        }
      }
    }
  }
  // if (std::abs(total[0]) >= 1e-6 || std::abs(total[1]) >= 1e-6) {
  //   std::cout << "total: " << total << "\n";
  // }
  return total;
}

template <class CellType>
Normal PUST<CellType>::solveFace(const double STin, const Pt& P0, const Pt& P1,
                                 const Pt& P2, const Pt& P3, const double delta,
                                 const double Pressure,
                                 const Normal& Marangoni) {
  // std::cout << "=========================== Solve Face \n" << std::endl;
  const double EPSILON = machine_epsilon<double>();
  // The Marangoni normal object holds the Xgradient, then the Y gradient, then
  // the temperature gradient of ST Marangoni = [Gx,Gy,sigma_T] Make Implicit
  // Surface std::cout << "In Solve Edge\n";

  // Something is up with paraboloids because they take like 8x the time to run.
  double STCoeff = STin;
  double intersection_threshold_m = 1e-6;  // This is the value of the wendland
                                           // function at 0.5 (if radius is 1).
  // The Pressure Option tells us if we should include the pressure terms or
  // not.
  double val0, val1, val2, val3, val4;
  // Copy Points for now to get around constant constraints
  Pt P0_copy = P0;
  Pt P1_copy = P1;
  Pt P2_copy = P2;
  Pt P3_copy = P3;
  // Calculate Corner Values
  val0 = this->getPU(P0);
  val1 = this->getPU(P1);
  val2 = this->getPU(P2);
  val3 = this->getPU(P3);
  val4 = this->getPU(0.25 * (P0 + P1 + P2 + P3));
  // Calculate Corner Weights
  double weight0, weight1, weight2, weight3, weight4;
  weight0 = this->getTotalWeight(P0);
  weight1 = this->getTotalWeight(P1);
  weight2 = this->getTotalWeight(P2);
  weight3 = this->getTotalWeight(P3);
  weight4 = this->getTotalWeight(0.25 * (P0 + P1 + P2 + P3));
  // If any of the weights are below A threshold, return 0 force since that
  // means we are not near the interface
  double factor = 100.0;
  if (weight0 < factor * EPSILON || weight1 < factor * EPSILON ||
      weight2 < factor * EPSILON || weight3 < factor * EPSILON) {
    // std::cout << "Low Weight Detected: " << weight0 << "," << weight1 << ","
    //           << weight2 << "," << weight3 << "\n";
    return Normal(0.0, 0.0, 0.0);
  }
  // std::cout << "Corner Values: " << val0 << "," << val1 << "," << val2 << ","
  // << val3 << "\n";
  // std::cout << "Corner Weights: " << weight0 << "," << weight1 << "," <<
  // weight2
  //           << "," << weight3 << "\n";
  Pt dP = P1 - P0;
  double D = std::sqrt(dP[0] * dP[0] + dP[1] * dP[1] + dP[2] * dP[2]);
  double denom = 1.0 / (safelyEpsilon(D));
  // Calculate case:
  int caseValue =
      1 * (-val0 > 0) + 2 * (-val1 > 0) + 4 * (-val2 > 0) + 8 * (-val3 > 0);
  // If all corners are outside or inside, return 0 force. Otehrwise, pair
  // intersections
  // std::cout << "Case Value: " << caseValue << "\n";
  std::vector<std::vector<int>> pairs = {{-1, -1}, {-1, -1}};
  // std::cout << "================== Case Value: " << caseValue << "\n";
  switch (caseValue) {
    case 0:
      return Normal(0.0, 0.0, 0.0);
      break;
    case 1:
      pairs[0] = {0, 3};
      break;
    case 2:
      pairs[0] = {1, 0};
      break;
    case 3:
      pairs[0] = {1, 3};
      break;
    case 4:
      pairs[0] = {2, 1};
      break;
    case 5:
      if (val4 > 0) {
        pairs[0] = {0, 1};
        pairs[1] = {2, 3};
      } else {
        pairs[0] = {0, 3};
        pairs[1] = {2, 1};
      }
      break;
    case 6:
      pairs[0] = {2, 0};
      break;
    case 7:
      pairs[0] = {2, 3};
      break;
    case 8:
      pairs[0] = {3, 2};
      break;
    case 9:
      pairs[0] = {0, 2};
      break;
    case 10:
      if (val4 > 0) {
        pairs[0] = {1, 2};
        pairs[1] = {3, 0};
      } else {
        pairs[0] = {1, 0};
        pairs[1] = {3, 2};
      }
      break;
    case 11:
      pairs[0] = {1, 2};
      break;
    case 12:
      pairs[0] = {3, 1};
      break;
    case 13:
      pairs[0] = {0, 1};
      break;
    case 14:
      pairs[0] = {3, 0};
      break;
    case 15:
      return Normal(0.0, 0.0, 0.0);
      break;
    default:
      break;
  }
  std::vector<std::vector<Pt>> intersectionsSet =
      std::vector<std::vector<Pt>>(4, std::vector<Pt>{});
  // Find all Intersections with face edges
  bool blocked01, blocked12, blocked23, blocked30;
  // std::cout << "Intersections01" << std::endl;
  std::vector<Pt> intersections01 =
      this->intersectEdge(P0, P1, 1, intersection_threshold_m, blocked01);
  // std::cout << "Intersections12" << std::endl;
  std::vector<Pt> intersections12 =
      this->intersectEdge(P1, P2, 1, intersection_threshold_m, blocked12);
  // std::cout << "Intersections23" << std::endl;
  std::vector<Pt> intersections23 =
      this->intersectEdge(P2, P3, 1, intersection_threshold_m, blocked23);
  // std::cout << "Intersections30" << std::endl;
  std::vector<Pt> intersections30 =
      this->intersectEdge(P3, P0, 1, intersection_threshold_m, blocked30);
  if (blocked01 || blocked12 || blocked23 || blocked30) {
    // std::cout << "Blocked Intersection Detected" << std::endl;
    return Normal(0.0, 0.0, 0.0);
  }
  intersectionsSet[0] = intersections01;
  intersectionsSet[1] = intersections12;
  intersectionsSet[2] = intersections23;
  intersectionsSet[3] = intersections30;

  // Calculate some edge and plane properties
  Pt dP1 = P1 - P0;
  Pt dP2 = P2 - P1;
  Normal edge1 = {dP1[0], dP1[1], dP1[2]};
  Normal edge2 = {dP2[0], dP2[1], dP2[2]};
  Normal faceNormal = IRL::crossProduct(edge1, edge2);
  faceNormal.normalize();

  // Give space for working variables
  Normal normal1, normal2, tangentStart, tangentEnd;
  Pt startPoint, endPoint, controlPoint;
  double weight = 1;  // Make Nonrational for now, but can add later if we want
  Normal total = {0.0, 0.0, 0.0};  // Total force

  // Loop over pairs, make splines, integrate, add forces to total.
  for (int j = 0; j < 2; j++) {
    if (pairs[j][0] != -1) {  // If there is a pair
      startPoint = intersectionsSet[pairs[j][0]][0];
      endPoint = intersectionsSet[pairs[j][1]][0];
      normal1 = this->getNormal(startPoint);
      normal1.normalize();
      normal2 = this->getNormal(endPoint);
      normal2.normalize();
      // Make tangent by crossing normal with face  normal
      tangentStart = IRL::crossProduct(faceNormal, normal1);
      tangentStart.normalize();
      tangentEnd = IRL::crossProduct(faceNormal, normal2);
      tangentEnd.normalize();

      RationalCubicBezierArc arc = RationalCubicBezierArc(
          startPoint, tangentStart, endPoint, tangentEnd);
      // std::cout
      //     << "===================== In Face Solver =====================\n";
      // std::cout << "Start Point: " << startPoint << "\n";
      // std::cout << "End Point: " << endPoint << "\n";
      // std::cout << "Tangent Start: " << tangentStart << "\n";
      // std::cout << "Tangent End: " << tangentEnd << "\n";
      // std::cout << "Control Point 1: " << arc.control_point_1() << "\n";
      // std::cout << "Control Point 2: " << arc.control_point_2() << "\n";
      // std::cout << "weight 1: " << arc.weight(1) << "\n";
      // std::cout << "weight 1: " << arc.weight(2) << "\n";
      int QuadRuleOrder = 50;
      const auto& abscissea = AbscissaeGauss<double, 50>();
      const auto& weights = WeightsGauss<double, 50>();
      for (int j = 0; j < QuadRuleOrder; ++j) {
        const double t = 0.5 * (1.0 + abscissea[j]);
        const double w = weights[j];

        // Evaluation Points
        Pt pt = arc.point(t);
        // get normal,tangent, and surface tension at each point.
        Pt tangent = arc.derivative(t);
        Normal tan = Normal(tangent[0], tangent[1], tangent[2]);
        double speed =
            std::sqrt(tan[0] * tan[0] + tan[1] * tan[1] + tan[2] * tan[2]);
        tan.normalize();

        Normal normal = this->getNormal(pt);
        normal.normalize();

        // here is where we get the surface tension coefficient
        // double ST = stencil_m.getScalar(pt);
        double ST = STCoeff + Marangoni[2] *
                                  (Marangoni[0] * pt[0] + Marangoni[1] * pt[1]);
        if (Marangoni[2] == -1.0 &&
            (Marangoni[0] != 0.0 ||
             Marangoni[1] != 0.0)) {  // Force droplet breakup
          double gamma0 = 1.0;        // Base surface tension coefficient
          double R = 0.5;             // Characteristic length scale

          // Linear decrease in surface tension with x-coordinate - Al-Saud
          // Style
          ST = gamma0 * std::max(1 - 1.25 * std::abs(pt[0]) / R, 0.1);
        }
        Normal f = IRL::crossProduct(tan, normal);
        f.normalize();
        f = f * ST;
        total = total + 0.5 * w * f * denom * speed;
      }

      // Debug Prints
      if (false) {  // False added as a toggle for the debug. Make true to
                    // turn on debug.
        // Arc Setup Info
        std::cout << "Start Point: " << startPoint << "\n";
        std::cout << "End Point: " << endPoint << "\n";
        std::cout << "Control Point 1: " << arc.control_point_1() << "\n";
        std::cout << "Control Point 2: " << arc.control_point_2() << "\n";
        std::cout << "pairs[j][0]: " << pairs[j][0] << "\n";
        std::cout << "pairs[j][1]: " << pairs[j][1] << "\n";
        for (int i = 0; i < intersectionsSet.size(); i++) {
          std::cout << "intersectionsSet[" << i << "]: ";
          if (!intersectionsSet[i].empty()) {
            std::cout << "(" << intersectionsSet[i][0] << ")";
          } else {
            std::cout << "(empty)";
          }
          std::cout << "\n";
        }
        std::cout << "P0: " << P0 << "\n";
        std::cout << "P1: " << P1 << "\n";
        std::cout << "P2: " << P2 << "\n";
        std::cout << "P3: " << P3 << "\n\n";
        // Output Info
        std::cout << "Total: " << total << "\n";
        std::cout << "===================================================== \n";
      }
    }
  }
  return -total;
}
// Print Solver
template <class CellType>
void PUST<CellType>::printSolver(void) {
  this->printSurface();
}

}  // namespace IRL

#endif  // IRL_PARTITION_OF_UNITY_SURFACE_TENSION_TPP_