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

// Solve Edge
template <class CellType>
Normal PUST<CellType>::solveEdge(const double STCoeff, const Pt& P0,
                                 const Pt& P1, const double delta,
                                 const double Pressure,
                                 const Normal& Marangoni) {
  // The Marangoni normal object holds the Xgradient, then the Y gradient, then
  // the temperature gradient of ST Marangoni = [Gx,Gy,sigma_T] Make Implicit
  // Surface std::cout << "In Solve Edge\n";

  // Something is up with paraboloids because they take like 8x the time to run.
  double STCoeff = STin;
  // The Pressure Option tells us if we should include the pressure terms or
  // not.
  PUImplicitSurface s = this->neighborhoodToImplicitSurface(delta);
  // Find Intersection with edge
  bool blocked = false;
  std::vector<Pt> intersections =
      s.intersectEdge(P0, P1, 10, intersection_threshold_m, blocked);

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
      tangent = s.getTangent(intersections[j]);
      tangent.normalize();
      // Here we apply the Marangoni surface tension. To do this, we assume
      // STCoeff = STCoeff + gamma_T(T-T_0). Letting gamma_T=-0.002 (Ratio for
      // water). We also pick T-T0=Gx to be the form we use. This will give us
      // th surface tension at this point We pick G = 10 for the stake of
      // simplicity
      STCoeff = STCoeff + Marangoni[2] * (Marangoni[0] * intersections[j][0] +
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

        double curv = s.getCurvature(intersections[j]);
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

}  // namespace IRL