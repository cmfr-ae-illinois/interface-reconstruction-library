#include <array>
#include <iostream>
#include "irl/splines/Spline.h"
#include <math.h>
#include <nlopt.hpp>
#include <random>

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty()) {
        for(int i=0; i < x.size(); i++) {
          grad[i] = 2*x[i];
        }
    }

    double obj = 0;
    for(int i=0; i<x.size(); i++) {
      obj += x[i]*x[i];
    }
    std::cout << obj << "\n";
    return obj;
}



int main() {
  
  // Define Rectangle for Clipping
  std::vector<std::vector<double>> square = {
      {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.0, 0.0}};

  // Pick Interpolating Points and Tangents (Vibes)
  std::vector<std::vector<double>> InterpolatingPoints = {
    {-1.0, 0.0}, {0.5, -1.0}, {1.5, 0.5},
    {0.75, 0.6}, {0.25, 0.5}, {-1.0, 0.0}};
  std::vector<std::vector<double>> Tangents = {
    {-0.555, -0.832}, {1.0, 0},        {-0.196, 0.981},
    {-0.707, -0.707}, {-0.707, 0.707}, {-0.555, -0.832}};
    
  // Circle
  double R = 0.999999;
  std::vector<std::vector<double>> Circle =
  {{R,0.0},{0.0,R},{-R,0.0},{0.0,-R},{R,0.0}};
  std::vector<std::vector<double>> CircleT =
  {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};


  // InterpolatingPoints = Circle;
  // Tangents = CircleT;
  // Interpolate
  Spline s = Spline<double>::LocalRQuadInterp(InterpolatingPoints, Tangents);
  s.saveToVTK("FullSection");
  // Calculate Parameter Loops
  std::vector<std::vector<double>> ret = s.getParameterLoop(square);
  std::vector<double> parameter =
      ret[0];  // Parameter value, dependent on indictator
  std::vector<double> indicator =
      ret[1];  // 1 = corner (outside), 2 = corner (inside), 3 = Curve hit
               // (entry), 4 = curve hit (Exit)

  // We want to follow Bezier curve when we go from exit to entry (Indicator
  // goes from 4 to 3). We will exact this part of the curve
  std::vector<std::vector<double>> Points;
  std::vector<std::vector<double>> tans;
  auto breakpoints = s.getBreakpoints();
  auto spans = s.getSpans();
  double u1,u2;
  for (int i = 0; i < parameter.size() - 1; i++) {
    double ind1 = indicator[i];      // Indicator at start
    double ind2 = indicator[i + 1];  // Indicator at end

    if (fabs(ind1 - 4) < 1e-8 &&
        fabs(ind2 - 3) < 1e-8) {  // 4 -> 3 case we want to look at

      // Inside this loop is taken directly from the integrateSplineSquare
      // function, but just changed a bit such that it makes vectors containing
      // the points and weights instead of the area contribution of this segment

      // Get Parameter Values
      u1 = parameter[i];
      u2 = parameter[i + 1];
      Points = {s.makeRationalQuadCurve(
          {u1})[0]};              // Add Starting Point to Points (Exit Point)
      tans = {s.getTangent(u1)};  // Add starting point Tangent (Exit Point)

      // Find Spans, in order
      int spanIndex1 = s.findSpan(u1);
      int spanIndex2 = s.findSpan(u2);
      // std::cout << "spanIndex1 = " << spanIndex1 << "\n";
      // std::cout << "spanIndex2 = " << spanIndex2 << "\n";
      // if u1,u2 in the same span, we can go direct.
      // If u1,u2 are in different spans, we have to segment through each span
      if (spanIndex1 == spanIndex2) {  // Same Span
        // I know this won't happen so I didn't code it for now
        continue;
      } else {  // Different Spans

        int numSpans = spans.size();

        int tempSpanIndex2 = spanIndex2;
        if (spanIndex2 < spanIndex1) {
          tempSpanIndex2 = spanIndex2 + numSpans;
        }

        // Get Set of Breakpoint Indices
        std::vector<int> breakIndexSet = {spanIndex1 + 1};
        for (int j = spanIndex1 + 2; j <= tempSpanIndex2; j++) {
          breakIndexSet.insert(breakIndexSet.end(), j);
        }

        // Mod them into range
        for (int j = 0; j < breakIndexSet.size(); j++) {
          breakIndexSet[i] = breakIndexSet[i] % numSpans;
        }
        // Calculate Break Values
        std::vector<double> breaks = {u1};
        for (int j = 0; j < breakIndexSet.size(); j++) {
          breaks.insert(breaks.end(), breakpoints[breakIndexSet[j]]);
        }

        breaks.insert(breaks.end(),
                      u2);  // This contains the parameter value of the start,
                            // all the breakpoints between
        // segments, then the end.

        for (int j = 0; j < breaks.size() - 1; j++) {
          // For each of these points, add the intersection and end point, along
          // with the previous tangent
          std::vector<double> Pend =
              s.makeRationalQuadCurve({breaks[j + 1]})[0];
          std::vector<double> Tend = s.getTangent(breaks[j + 1]);

          std::vector<double> Pstart = s.makeRationalQuadCurve({breaks[j]})[0];
          std::vector<double> Tstart = s.getTangent(breaks[j]);
          // std::cout << "j = " << j << "==================\n";
          // std::cout << "Pstart = " << Pstart[0] << "," << Pstart[1] << "\n";
          // std::cout << "Pend = " << Pend[0] << "," << Pend[1] << "\n";
          // std::cout << "Tstart = " << Tstart[0] << "," << Tstart[1] << "\n";
          // std::cout << "Tend = " << Tend[0] << "," << Tend[1] << "\n";

          // Calculate Intersection
          std::vector<std::vector<double>> solution =
              Spline<double>::solvePointTangentIntersection(Pstart, Pend, Tstart, Tend);
          std::vector<double> inter = solution[0];

          // Now, add intersection and end point to array
          Points.insert(Points.end(), inter);
          Points.insert(Points.end(), Pend);

          tans.insert(tans.end(),Tend);  // Put end tangent to back of this to move forward.
        }
      }
    }
  }
  // Print Points - These are the control points for the bezier splines. They go
  // start point, intersection point, end point.
  // std::cout << "\n======== Internal Points Result ======== \n";
  std::vector<std::vector<double>> P2 = s.clippedBezier(u1,u2);
  std::cout << "\n======== Points Result ======== \n";
  for (int i = 0; i < Points.size(); i++) {
    std::cout << Points[i][0] << "," << Points[i][1] << "\n";
  }
  
  // From here, we will use this information to calculate the  weights.
  std::vector<double> weights = {
      1};  // First weight is always one (every other one will be too);
  for (int i = 0; i < Points.size() - 1;
       i += 2) {  // We go every 2 since the points go in triples (0,1,2 then
                  // 2,3,4 then 4,5,6 etc.)

    // Here is an example of how to use the Points array, and what the ordering
    // of it is
    std::vector<double> P0 = Points[i];      // Start
    std::vector<double> P1 = Points[i + 1];  // Intersection
    std::vector<double> P2 = Points[i + 2];  // End

    // Since I don't fully understand the curvature method in paper yet, I am
    // just quickly doing this to get weights Since I know they are well
    // behaved. I will understand the method in the paper soon and then that can
    // be used In general case.

    // First, calculate midpoint of P0,P2
    std::vector<double> mid = {(P0[0] + P2[0]) / 2, (P0[1] + P2[1]) / 2};
    // Next, calculate intersection with spline (Shoulder Point)
    double uHit = s.lineCurveIntersection(mid, P1)[0];
    std::vector<double> S = s.makeRationalQuadCurve({uHit})[0];
    // Calculate Midpoint to Shoulder Point Distance
    double MS = sqrt(pow(S[0] - mid[0], 2) + pow(S[1] - mid[1], 2));
    // Calculate Shoulder Point to intersection distance
    double SP1 = sqrt(pow(S[0] - P1[0], 2) + pow(S[1] - P1[1], 2));
    // Calculate Weight (eq 7.32)
    double w = MS / SP1;
    // Add Weight
    weights.insert(weights.end(), w);
    weights.insert(weights.end(), 1);  // Always 1 at the start and end
    // std::cout << "weight = " << w << "\n";
    // std::cout << "weight should = " << s.makeWeight(P0,P1,P2) << "\n";
  }
  
  // Printing Results
  std::cout << "\n======= Weight Result ========== \n";
  // Notice that, starting at the first, every other weight is 1. This is just a
  // choice. The other weights (index 1,3,5,...) are the ones which should be
  // used when converting to the rational Bezier curves of the form in the
  // paper.
  for (int i = 0; i < weights.size(); i++) {
    std::cout << weights[i] << ",";
  }
  std::cout << "\n";

  // std::cout << "\n======== Points 2 Result ======== \n";
  // for (int i = 0; i < P2.size(); i++) {
  //   std::cout << P2[i][0] << "," << P2[i][1] << "," << P2[i][2]<< "\n";
  // }

  double ClipArea = s.integrateSplineSquare(square);
  // Finally, we will print out some properties of the curve
  std::cout << "\n========== Curve Properties ========== \n";
  std::cout << "Total Curve Area: " << s.getArea() << "\n";
  std::cout << "Clipped Area: " << ClipArea << "\n";
  std::cout << "Arc Length: " << s.getArcLength() << "\n";
  std::cout << "Surface Energy: " << s.getSurfaceEnergy() << "\n";

  // Here we will make the spline for a visualization Test
  // We just need to make a knot vector
  // Make Knot Vector
  
  /////////////////////////////////////////// TEST BY FABIEN
  // First let's close boundary of area to compute.
  // We'll define straight lines are quadratic Bezier curves with 0 weight.
  // That's an overkill but quicker to handle for now.
  const auto Ps = std::vector<double>({Points[0][0], Points[0][1]}); // Start Point
  const auto Pe = std::vector<double>( // End Point
      {Points[Points.size() - 1][0], Points[Points.size() - 1][1]});
  const auto P00 = std::vector<double>({0.0, 0.0}); // Bottom Left
  const auto P10 = std::vector<double>({1.0, 0.0}); // Bottom Right

  weights.push_back(0.0);
  Points.push_back(
      std::vector<double>({0.5 * (Pe[0] + P00[0]), 0.5 * (Pe[1] + P00[1])}));

  weights.push_back(1.0);
  Points.push_back(P00);

  weights.push_back(0.0);
  Points.push_back(
      std::vector<double>({0.5 * (P00[0] + P10[0]), 0.5 * (P00[1] + P10[1])}));

  weights.push_back(1.0);
  Points.push_back(P10);

  weights.push_back(0.0);
  Points.push_back(
      std::vector<double>({0.5 * (Ps[0] + P10[0]), 0.5 * (Ps[1] + P10[1])}));


  std::cout << "\n======== Closed contour result ======== \n";
  for (int i = 0; i < Points.size(); i += 2) {
    std::cout << std::scientific << std::setprecision(2) << "P0 = ("
              << Points[i][0] << "," << Points[i][1] << "); ";
    std::cout << "P1 = (" << Points[i + 1][0] << "," << Points[i+1][1] << "); ";
    std::cout << "P2 = (" << Points[(i + 2) % Points.size()][0] << ","
              << Points[(i + 2) % Points.size()][1] << "); ";
    std::cout << " w = " << weights[i + 1] << ")\n";
  }
  std::vector<double> U = {0,0,0};
  std::cout << "Point Size = " << Points.size() << "\n";
  int ubarSize = (Points.size())/2;
  for(int i = 1;i < ubarSize;i++) {
      U.insert(U.end(),{double(i+1)/ubarSize,double(i+1)/ubarSize});
  }
  U.insert(U.end(),1);
  std::vector<std::vector<double>> PControl = Points;
  PControl.insert(PControl.end(),Points[0]);
  Spline clipped = Spline<double>(PControl,U,weights);
  std::cout<< "Clipped Area 2 = " << clipped.getArea() << "\n"; 
  clipped.saveToVTK("ClippedSection");
  

  std::cout << "\n======== Computing area with exact formula ======== \n";
  double A1 = 0.0;
  for (int i = 0; i < Points.size(); i += 2) {
    const auto x0 = Points[i][0], y0 = Points[i][1];
    const auto x1 = Points[i + 1][0], y1 = Points[i+1][1];
    const auto x2 = Points[(i + 2) % Points.size()][0],
               y2 = Points[(i + 2) % Points.size()][1];
    const auto w = weights[i + 1];

    auto K = Spline<double>::coeffsAreaExact(w);
    A1 -= (x0 * y2 - x2 * y0) * K[0] +
          (x1 * y0 + x2 * y1 - x0 * y1 - x1 * y2) * K[1];
  }
  std::cout << "A = " << A1 << std::endl;

  std::cout << "\n======== Computing area with exact + Taylor series formula "
               "======== \n";
  double A2 = 0.0;
  for (int i = 0; i < Points.size(); i += 2) {
    const auto x0 = Points[i][0], y0 = Points[i][1];
    const auto x1 = Points[i + 1][0], y1 = Points[i+1][1];
    const auto x2 = Points[(i + 2) % Points.size()][0],
               y2 = Points[(i + 2) % Points.size()][1];
    const auto w = weights[i + 1];
    double dA;
    
    // std::cout << "P0 = " << x0 << "," << y0 << ", w = " << weights[i] << "\n";
    //             std::cout << "P1 = " << x1 << "," << y1 << ", w = " << weights[i+1]<< "\n";
    //             std::cout << "P2 = " << x2 << "," << y2 << ", w = " << weights[i+2]<< "\n";
    if (w < 0.35 || w > 1.7) {
      auto K = Spline<double>::coeffsAreaExact(w);
      dA = (x0 * y2 - x2 * y0) * K[0] +
            (x1 * y0 + x2 * y1 - x0 * y1 - x1 * y2) * K[1];
    } else {
      auto K = Spline<double>::coeffsAreaSeries(w);
      dA = (x0 * y2 - x2 * y0) * K[0] +
            (x1 * y0 + x2 * y1 - x0 * y1 - x1 * y2) * K[1];
    }
    dA = -dA;
    A2 += dA;
    // std::cout << "dA = " << dA <<"\n";
    // std::cout << "weight = " << w << "\n";
  }
  std::cout << "A = " << A2 << std::endl;

  double Exact = M_PI*R*R/4.0;
  std::cout << "\n======== Error ======== \n";
  std::cout << "Diff = " << (A1 - A2) / A2 << std::endl;
  std::cout << "Original - A2= " << (ClipArea-A2) / A2 << std::endl;
  std::cout << "Original - A1= " << (ClipArea-A1) / A1 << std::endl;
  std::cout << "Original - Exact= " << (ClipArea-Exact) / Exact << std::endl;
  std::cout << "A1 - Exact= " << (A1-Exact) / Exact << std::endl;
  std::cout << "A2 - Exact= " << (A2-Exact) / Exact << std::endl;

  /////////////////////////////////////////// TEST BY FABIEN



  // NLOPT TESTING
  int N = 500000000;
  nlopt::opt opt(nlopt::LD_MMA, N);
  opt.set_min_objective(myvfunc, NULL);
  opt.set_xtol_rel(1e-4);
  opt.set_ftol_abs(1e-12);
  std::vector<double> x(N);
  for(int i = 0; i < N; i++) {
    x[i] = 0.5 * pow(-1,i);
  }
double minf;
std::cout <<"\n";
try{
    nlopt::result result = opt.optimize(x, minf);
    std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
        << std::setprecision(10) << minf << std::endl;
    std::cout << "X Result: \n" << x.size();
}
catch(std::exception &e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
}
  return 0;
}