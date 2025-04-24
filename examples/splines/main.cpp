#include <iostream>
#include "Spline.h"

int main() {
    // std::vector<std::vector<double>> V1 = {{1,2},{0.7239,3.3806},{2.1078,2},{3.6852,0.4262},{4,2}};
    // std::vector<double> KV = {0,0,0,0.6667,0.6667,1,1,1};
    // std::vector<double> W = {1,0.3361,1,0.5025,1};

    // std::vector<std::vector<double>> Circle = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    // std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    // std::vector<std::vector<double>> InterpTester = {{1.0,0.0},{0.0,1.0},{-1.0,0.0},{-3.0,0.0},{1.0,0.0}};
    // std::vector<std::vector<double>> InterpTesterT = {{0.0,1.0},{-1.0,0.0},{-1.0,0.0},{-1.0,0.0},{0.0,1.0}};
    
    // Spline spline = Spline(V1,KV,W);
    // Spline circle = Spline::LocalRQuadInterp(Circle,CircleT);

    
    // std::cout << "================== Being Circle test ===========================\n";
    // std::cout << "Total Area = " <<circle.getArea() << "\n";
    // std::cout << "Total Arc Length = " << circle.getArcLength() << "\n";
    // std::cout << "Surface Energy = " << circle.getSurfaceEnergy() << "\n";
    // std::vector<std::vector<double>> square = {{{0,0},{1.5,0},{1.5,1.5},{0,1.5},{0,0}}};
    // double A = circle.integrateSplineSquare(square);
    // std::cout << "Clipped Area = " << A << "\n";

    // circle.saveToVTK("circle_test");
    // // Circle Tests Working and Matching


    // std::cout << "\n================== Being Weird shape test ===========================\n";
    // Spline InterpTest = Spline::LocalRQuadInterp(InterpTester,InterpTesterT);

    // InterpTest.saveToVTK("InterpTester_test");
    // std::cout << "Total Area = " <<InterpTest.getArea() << "\n";
    // std::cout << "Total Arc Length = " << InterpTest.getArcLength() << "\n";
    // std::cout << "Surface Energy = " << InterpTest.getSurfaceEnergy() << "\n";
    // square = {{-4,-2},{-2,-2},{-2,0.5},{-4,0.5},{-4,-2}};
    // A = InterpTest.integrateSplineSquare(square);
    // std::cout << "Clipped Area = " << A << "\n";

    // Define Rectange for Clipping
    std::vector<std::vector<double>> square = {{0.0,0.0},{1.0,0.0},{1.0,1.0},{0.0,1.0},{0.0,0.0}}; 

    // Pick Interpolating Points and Tangents (Vibes)
    std::vector<std::vector<double>> InterpolatingPoints = {{-1.0,0.0},{0.5,-1.0},{1.5,0.5},{0.75,0.6},{0.25,0.5},{-1.0,0.0}};
    std::vector<std::vector<double>> Tangents = {{-0.555,-0.832},{1.0,0},{-0.196,0.981},{-0.707,-0.707},{-0.707,0.707},{-0.555,-0.832}};

    // Interpolate
    Spline s = Spline::LocalRQuadInterp(InterpolatingPoints,Tangents);

    // Calculate Parameter Loops
    std::vector<std::vector<double>> ret = s.getParameterLoop(square);
    std::vector<double> parameter = ret[0]; // Parameter value, dependent on indictator
    std::vector<double> indicator = ret[1]; // 1 = corner (outside), 2 = corner (inside), 3 = Curve hit (entry), 4 = curve hit (Exit)

    // We want to follow Bezier curve when we go from exit to entry (Indicator goes from 4 to 3). We will exact this part of the curve
    std::vector<std::vector<double>> Points;
    std::vector<std::vector<double>> tans;
    auto breakpoints = s.getBreakpoints();
    auto spans = s.getSpans();
    for(int i = 0; i < parameter.size()-1;i++) {
        double ind1 = indicator[i]; // Indicator at start
        double ind2 = indicator[i+1]; // Indicator at end

        if(fabs(ind1-4) < 1e-8 && fabs(ind2 - 3) < 1e-8) { // 4 -> 3 case we want to look at

            // Inside this loop is taken directly from the integrateSplineSquare function, but just changed a bit 
            // such that it makes vectors containing the points and weights instead of the area contribution of this segment

            // Get Parameter Values
            double u1 = parameter[i];
            double u2 = parameter[i+1];
            Points = {s.makeRationalQuadCurve({u1})[0]}; // Add Starting Point to Points (Exit Point)
            tans = {s.getTangent(u1)}; // Add starting point Tangent (Exit Point)

            // Find Spans, in order
            int spanIndex1 = s.findSpan(u1);
            int spanIndex2 = s.findSpan(u2);
            // std::cout << "spanIndex1 = " << spanIndex1 << "\n";
            // std::cout << "spanIndex2 = " << spanIndex2 << "\n";
            // if u1,u2 in the same span, we can go direct.
            // If u1,u2 are in different spans, we have to segment through each span
            if(spanIndex1 == spanIndex2) { // Same Span
                // I know this won't happen so I didn't code it for now
                continue;
            } else { // Different Spans
                
                int numSpans = spans.size();

                int tempSpanIndex2 = spanIndex2;
                if(spanIndex2 < spanIndex1) {
                    tempSpanIndex2 = spanIndex2 + numSpans;
                }

                // Get Set of Breakpoint Indices
                std::vector<int> breakIndexSet = {spanIndex1+1};
                for(int j = spanIndex1+2;j<=tempSpanIndex2; j++) {
                    breakIndexSet.insert(breakIndexSet.end(),j);
                }
                
                // Mod them into range
                for(int j = 0; j < breakIndexSet.size();j++ ) {
                    breakIndexSet[i] = breakIndexSet[i]%numSpans;
                }
                // Calculate Break Values
                std::vector<double> breaks = {u1};
                for(int j = 0; j < breakIndexSet.size();j++) {
                    breaks.insert(breaks.end(),breakpoints[breakIndexSet[j]]);
                }
                
                breaks.insert(breaks.end(),u2); // This contains the parameter value of the start, all the breakpoints between
                // segments, then the end.
                
                for(int j = 0; j < breaks.size()-1; j++) {
                    // For each of these points, add the intersection and end point, along with the previous tangent
                    std::vector<double> Pend = s.makeRationalQuadCurve({breaks[j+1]})[0];
                    std::vector<double> Tend = s.getTangent(breaks[j+1]);
                    
                    std::vector<double> Pstart = s.makeRationalQuadCurve({breaks[j]})[0];
                    std::vector<double> Tstart = s.getTangent(breaks[j]);
                    // std::cout << "j = " << j << "==================\n";
                    // std::cout << "Pstart = " << Pstart[0] << "," << Pstart[1] << "\n";
                    // std::cout << "Pend = " << Pend[0] << "," << Pend[1] << "\n";
                    // std::cout << "Tstart = " << Tstart[0] << "," << Tstart[1] << "\n";
                    // std::cout << "Tend = " << Tend[0] << "," << Tend[1] << "\n";

                    // Calculate Intersection
                    std::vector<std::vector<double>> solution = Spline::solvePointTangentIntersection(Pstart,Pend,Tstart,Tend);
                    std::vector<double> inter = solution[0];

                    // Now, add intersection and end point to array
                    Points.insert(Points.end(),inter);
                    Points.insert(Points.end(),Pend);

                    tans.insert(tans.end(),Tend); // Put end tangent to back of this to move forward.
                }
            }
        }
    }
    // Print Points - These are the control points for the bezier splines. They go start point, intersection point, end point.
    std::cout << "\n======== Points Result ======== \n";
    for(int i = 0; i < Points.size(); i++) {
        std::cout << Points[i][0] << "," << Points[i][1] << "\n";
    }

    // From here, we will use this information to calculate the  weights.
    std::vector<double> weights = {1}; // First weight is always one (every other one will be too);
    for(int i = 0; i < Points.size()-1; i+= 2) { // We go every 2 since the points go in triples (0,1,2 then 2,3,4 then 4,5,6 etc.)
        
        // Here is an example of how to use the Points array, and what the ordering of it is
        std::vector<double> P0 = Points[i]; // Start
        std::vector<double> P1 = Points[i+1]; // Intersection
        std::vector<double> P2 = Points[i+2]; // End

        // Since I don't fully understand the curvature method in paper yet, I am just quickly doing this to get weights 
        // Since I know they are well behaved. I will understand the method in the paper soon and then that can be used
        // In general case.

        // First, calculate midpoint of P0,P2
        std::vector<double> mid = {(P0[0]+P2[0])/2,(P0[1]+P2[1])/2};
        // Next, calculate intersection with spline (Shoulder Point)
        double uHit= s.lineCurveIntersection(mid,P1)[0]; 
        std::vector<double> S = s.makeRationalQuadCurve({uHit})[0];
        // Calculate Midpoint to Shoulder Point Distance
        double MS = sqrt(pow(S[0]-mid[0],2)+pow(S[1]-mid[1],2));
        // Calculate Shoulder Point to intersection distance
        double SP1 = sqrt(pow(S[0]-P1[0],2)+pow(S[1]-P1[1],2));
        // Calculate Weight (eq 7.32)
        double w = MS/SP1;
        // Add Weight
        weights.insert(weights.end(),w);
        weights.insert(weights.end(),1); // Always 1 at the start and end
    }

    // Printing Results
    std::cout << "\n======= Weight Result ========== \n";
    for(int i = 0; i < weights.size(); i++) {
        std::cout << weights[i] << ",";
    }
    std::cout << "\n";

    // Finally, we will print out some properties of the curve
    std::cout << "\n========== Curve Properties ========== \n";
    std::cout << "Total Curve Area: " << s.getArea() << "\n";
    std::cout << "Clipped Area: " << s.integrateSplineSquare(square) << "\n";
    std::cout << "Arc Length: " << s.getArcLength() << "\n";
    std::cout << "Surface Energy: " << s.getSurfaceEnergy() << "\n";

    return 0;
}
