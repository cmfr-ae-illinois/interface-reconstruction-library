#include <iomanip>
#include <iostream>
#include <vector>
#include <math.h>
#include <nlopt.hpp>
// #include "examples/spline_optimization/Spline.h"
#include "irl/splines/Spline.h"
#include "examples/spline_optimization/basic_mesh.h"
#include "examples/spline_optimization/irl2d.h"
#include "examples/spline_optimization/data.h"
#include "irl/splines/MarchingSquares.h"


// This is data for our constraints
// typedef struct {
//     double A;
//     int N;
// } my_data;

typedef struct {
    double A;
} my_constraint_data;

typedef struct {
    std::vector<std::vector<double>> square;
    double VOF;
} my_square_data;

typedef struct {
    std::vector<std::vector<std::vector<double>>> squares;
    double tolerance;
} my_empty_data;

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
    std::cout << "Call my Func\n";
    my_constraint_data *d = (my_constraint_data *) my_func_data;
    double A = d->A;
    std::cout << "A = " << A << "\n";
    // std::cout << "Unpacking\n";
    std::vector<std::vector<std::vector<double>>> ret = Spline<double>::unpack(x);
    // std::cout << "Interpolating \n";
    Spline s = Spline<double>::LocalRQuadInterp(ret[0],ret[1]);
    // std::cout << "Interpolated \n";
    // std::cout << "Getting Surface Energy \n";
    double E = s.getSurfaceEnergy();
    // std::cout << "Got Surface Energy \n";
    // std::cout << "Getting Arc Length \n";
    double Al = s.getArcLength();
    // std::cout << "Got Arc Length \n";
    double AlFactor = 2*std::sqrt(M_PI*A);
    double Test = 1.0/AlFactor;
    std::cout << "ALFactor = " << AlFactor << "\n";
    std::cout << "E ===" << E/(2*M_PI) << "\n";
    std::cout << "Al ===" << Al << "\n";
    std::cout << "Test ===" << Test << "\n";
    std::cout << "AlContrib ===" << Al*Test << "\n";

    s.saveToVTK("LastStepResult");
    std::cout << "======" << E/(2*M_PI) + Al/AlFactor << "\n";
    double objective = E/(2*M_PI) + Al/AlFactor;
    if(objective != objective) {
        s.saveToVTK("nanObjective");
        std::cout << "nanObjective\n";
        s.printControlPoints();
        objective = 1000;
    }
    // Calculate Gradient with Finite Differences

    // std::cout <<"Grad Size" << grad.size() << "\n";
    std::vector<double> xCopy(x.size());
    double nudge = 1e-6;
    if(!grad.empty()) {
        std::cout <<"In Gradient\n"; 
        // Loop Over Gradient Components
        for(int j = 0; j < grad.size(); j++) {
            // Copy X
            for(int i = 0; i < xCopy.size(); i++) {
                if(i == j) {
                    xCopy[i] = x[i] + nudge;
                } else {
                    xCopy[i] = x[i];
                }
            }

            // Unpack, Interpolate
            ret = Spline<double>::unpack(xCopy);
            Spline s2 = Spline<double>::LocalRQuadInterp(ret[0],ret[1]);
            // Get Properties
            double E2 = s2.getSurfaceEnergy();
            double Al2 = s2.getArcLength();

            // Objective
            double objective2 = E2/(2*M_PI) + Al2/AlFactor;
            if(objective2 != objective2) {
                objective2 = 1000;
            }
            // Calculate Derivative
            // std::cout << "Grad " << j << " = " << (objective2-objective)/nudge << "\n";
            grad[j] = (objective2-objective)/nudge;
        }
    }
    // Show x
    std::cout << "x = [";
    for(int i = 0; i < x.size(); i++) {
        std::cout << x[i] << ",";
    }
    std::cout << "]\n";
    std::cout << "End my Func = " << objective << "\n";
    
    return objective;
}

// This is the constraint function
double myVOFconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data) { // VOF Constraint
    // std::cout << "Call VOF Constraint \n";
    my_square_data *d = (my_square_data *) data;
    std::vector<std::vector<double>> square = d->square; // 5 Point Rectangle
    double VOF = d->VOF;
    // std::cout << "dx, dy get\n";
    double dx = fabs(square[0][0]-square[2][0]);
    // std::cout << "dx done\n";
    double dy = fabs(square[0][1]-square[2][1]);
    // std::cout << "dy done\n";
    
    std::vector<std::vector<std::vector<double>>> ret = Spline<double>::unpack(x);
    std::vector<std::vector<double>> Q = ret[0];
    std::vector<std::vector<double>> T = ret[1];
    
    Spline s = Spline<double>::LocalRQuadInterp(Q,T);
    // std::cout << "Finding Intersection\n";
    // s.saveToVTK("OptimizationSplineTest");
    double Acurr = s.integrateSplineSquare(square);
    // std::cout << "Found\n";
    // std::cout << "Done\n";
    // std::cout <<"===========================" << (Acurr/(dx*dy) - VOF) << "\n";
    double obj = (Acurr/(dx*dy) - VOF);

    double nudge = 1e-6;
    std::vector<double> xCopy(x.size());
    if(!grad.empty()) {
        // Loop Over Gradient Components
        for(int j = 0; j < grad.size(); j++) {
            // Copy X
            for(int i = 0; i < xCopy.size(); i++) {
                if(i == j) {
                    xCopy[i] = x[i] + nudge;
                } else {
                    xCopy[i] = x[i];
                }
            }

            // Unpack, Interpolate
            ret = Spline<double>::unpack(xCopy);
            Spline s2 = Spline<double>::LocalRQuadInterp(ret[0],ret[1]);
            // Get Properties
            double A2 = s2.integrateSplineSquare(square);

            // Objective
            double obj2 = (A2/(dx*dy) - VOF);

            // Calculate Derivative
            grad[j] = (obj2-obj)/nudge;
        }
    }

    // std::cout << "End VOF Constraint \n";

    return std::abs(obj); // Calculated VOF - Exact VOF
}

double myEmptyconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data) { // VOF Constraint
    std::cout << "Call Empty Constraint \n";
    my_empty_data *d = (my_empty_data *) data;
    std::vector<std::vector<std::vector<double>>> squares = d->squares; // A vector of 5 point rectangles
    double tol = d->tolerance;
    
    std::vector<std::vector<std::vector<double>>> ret = Spline<double>::unpack(x);
    std::vector<std::vector<double>> Q = ret[0];
    std::vector<std::vector<double>> T = ret[1];
    
    Spline s = Spline<double>::LocalRQuadInterp(Q,T);
    double result = 0;
    for(int i = 0; i < squares.size(); i++) {
        std::vector<std::vector<double>> square = squares[i];
        // std::cout << "dx, dy get\n";
        double dx = fabs(square[0][0]-square[2][0]);
        // std::cout << "dx done\n";
        double dy = fabs(square[0][1]-square[2][1]);

        double A = s.integrateSplineSquare(square);

        result += A;
    }

    double nudge = 1e-6;
    std::vector<double> xCopy(x.size());
    if(!grad.empty()) {
        // Loop Over Gradient Components
        for(int j = 0; j < grad.size(); j++) {
            // Copy X
            for(int i = 0; i < xCopy.size(); i++) {
                if(i == j) {
                    xCopy[i] = x[i] + nudge;
                } else {
                    xCopy[i] = x[i];
                }
            }

            // Unpack, Interpolate
            ret = Spline<double>::unpack(xCopy);
            Spline s2 = Spline<double>::LocalRQuadInterp(ret[0],ret[1]);
            // Get Properties
            // double A2 = s2.integrateSplineSquare(square);

            // Objective
            double obj2 = 0;

            // Calculate Derivative
            grad[j] = (obj2-0)/nudge;
        }
    }

    std::cout << "End VOF Constraint \n";

    return result; // Calculated VOF - Exact VOF
}

double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data) { // Total Area
    std::cout << "Call Constraint \n";
    my_constraint_data *d = (my_constraint_data *) data;
    double A = d->A;
    
    std::vector<std::vector<std::vector<double>>> ret = Spline<double>::unpack(x);
    std::vector<std::vector<double>> Q = ret[0];
    std::vector<std::vector<double>> T = ret[1];
    Spline s = Spline<double>::LocalRQuadInterp(Q,T);
    double Acurr = s.getArea();
    double obj = (Acurr-A);

    double nudge = 1e-6;
    std::vector<double> xCopy(x.size());
    if(!grad.empty()) {
        // Loop Over Gradient Components
        for(int j = 0; j < grad.size(); j++) {
            // Copy X
            for(int i = 0; i < xCopy.size(); i++) {
                if(i == j) {
                    xCopy[i] = x[i] + nudge;
                } else {
                    xCopy[i] = x[i];
                }
            }

            // Unpack, Interpolate
            ret = Spline<double>::unpack(xCopy);
            Spline s2 = Spline<double>::LocalRQuadInterp(ret[0],ret[1]);
            // Get Properties
            double A2 = s2.getArea();

            // Objective
            double obj2 = (A2 - A);

            // Calculate Derivative
            grad[j] = (obj2-obj)/nudge;
        }
    }

    std::cout << "End Constraint \n";
    return -fabs(Acurr-A);
}


// ================================== BEGIN MAIN ===========================================================
int main() {

// Ignore this, WIP
    // Field Initialization
    const int a_nx = 25;
    std::vector<std::vector<double>> VOFarray(a_nx, std::vector<double>(a_nx,-1));
    const int GC = 3;
    // Set up Mesh
    double xmin = -0.5;
    double xmax = 0.5;
    double ymin = -0.5;
    double ymax = 0.5;
    BasicMesh mesh(a_nx, a_nx, GC);
    IRL2D::Vec my_lower_domain(xmin, ymin);
    IRL2D::Vec my_upper_domain(xmax, ymax);
    mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
    // auto dx = mesh.dx();
    std::cout << "========== Mesh Properties =============\n";
    std::cout << "dx = " << mesh.dx() << "\n";
    std::cout << "dy = " << mesh.dy() << "\n";
    for(int i =0; i <= a_nx; i++) {
        std::cout << mesh.x(i) << ",";
    }
    std::cout << "\n";
    // Initialize
    const auto circle_center = IRL2D::Vec(0,0);
    const double circle_radius = 0.15;

    // Data<IRL2D::Parabola> a_interface;
    // a_interface->updateBorder();
    // correctInterfaceBorders(a_interface);
  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
        const IRL2D::Vec lower_cell_pt(mesh.x(i), mesh.y(j));
        const IRL2D::Vec upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1));
        const IRL2D::Vec mid_pt = 0.5 * (lower_cell_pt + upper_cell_pt);
        IRL2D::Vec disp = mid_pt - circle_center;
        const auto mag = disp.magnitude();
        if (mag < circle_radius - 2.0 * mesh.dx()) {
            // (a_interface)(i, j).markAsAlwaysAbove();
            VOFarray[i][j] = 1;
        } else if (mag > circle_radius + 2.0 * mesh.dx()) {
            // (a_interface)(i, j).markAsAlwaysBelow();
            VOFarray[i][j] = 0;
        } else {
            auto circle_normal = IRL2D::Vec(disp);
            circle_normal.normalize();
            const IRL2D::Vec datum = circle_center + circle_radius * circle_normal;
            const IRL2D::ReferenceFrame frame(
                IRL2D::Vec(circle_normal[1], -circle_normal[0]), circle_normal);
            const double coeff = 1.0 / (2.0 * circle_radius);
            auto temp = IRL2D::Parabola(datum, frame, coeff);
            for (int ii = 0; ii < 10; ii++) {
                auto cell = IRL2D::RectangleFromBounds(lower_cell_pt, upper_cell_pt);
                auto moments = IRL2D::ComputeMoments(cell, temp);
                const double liquid_volume_fraction =
                    moments.m0() / IRL2D::ComputeArea(cell);
                VOFarray[i][j] = liquid_volume_fraction;
                if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
                    liquid_volume_fraction <= IRL::global_constants::VF_HIGH) { // Whats happening here???
                    const auto new_disp = moments.m1() / moments.m0() - circle_center;
                    auto new_circle_normal = IRL2D::Vec(new_disp);
                    new_circle_normal.normalize();
                    const IRL2D::Vec new_datum =
                        circle_center + circle_radius * new_circle_normal;
                    const IRL2D::ReferenceFrame new_frame(
                        IRL2D::Vec(new_circle_normal[1], -new_circle_normal[0]),
                        new_circle_normal);
                    // (a_interface)(i, j) = IRL2D::Parabola(new_datum, new_frame, coeff);
                } else {
                    break;
                }
            }
        }
      }
    }

    // Print Result
    // std::cout << "========== VOF INITIAL =========\n";
    // for(int i = 0; i < a_nx; i++) { // Note that i corresponds to x, j corresponds to y
    //     for(int j = 0; j < a_nx; j++) {
    //         std::cout << VOFarray[i][j] << ",           ";
    //     }
    //     std::cout << "\n";
    // }
    // At this point, we have the volume fractions. Now we need to set up the array of square
    // While Doing this, we will also seed our initial conditions.
    // The initial conditions on this first pass will be the midpoints of the cells.
    // The tangents will be from one point to the next.

    std::vector<my_square_data> VOFconstraintSet = {};
    std::vector<std::vector<std::vector<double>>> emptyConstraintSet = {};
    std::vector<std::vector<bool>> addArray(a_nx, std::vector<bool>(a_nx,false));
    std::vector<std::vector<double>> interpPoints = {};

    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
            // Make variable to see if we will add point
            bool add = false;
            // If any of the adjacent points, or itself, are nonzero, add the current square
            add = (fabs(VOFarray[i][j]) > IRL::global_constants::VF_LOW);
            // Here we take a moment to add the interpolation Points 
            // This process will need to be refined as we get more complex shapes.
            if(add) {
                double xmid = mesh.xm(i);
                double ymid = mesh.ym(j);

                interpPoints.push_back({xmid,ymid});
            }
            // add VOF constraints
            if(add) {
                // Get Square
                double xmin = mesh.x(i);
                double xmax = mesh.x(i+1);
                double ymin = mesh.y(i);
                double ymax = mesh.y(i+1);
                
                // Make Square
                std::vector<std::vector<double>> squ = {{xmin,ymin},{xmax,ymin}
                                                        ,{xmax,ymax},{xmin,ymax},{xmin,ymin}};
                // Get VOF
                double VOF = VOFarray[i][j];

                // Make Struct
                my_square_data temp = {squ,VOF};
                VOFconstraintSet.push_back(temp);
            }

            // Look for Boundary Cells that are empty
            add = false;
            if(i != 0) { // Check Left 
                add = ((fabs(VOFarray[i-1][j]) > IRL::global_constants::VF_LOW) && (fabs(VOFarray[i][j]) <= 1e-6)) || add;
            }
            if(i != mesh.imax()) { // Check Right
                add = ((fabs(VOFarray[i+1][j]) > IRL::global_constants::VF_LOW) && (fabs(VOFarray[i][j]) <= 1e-6)) || add;
            }
            if(j != 0) { // Check Up
                add = ((fabs(VOFarray[i][j-1]) > IRL::global_constants::VF_LOW) && (fabs(VOFarray[i][j]) <= 1e-6)) || add;
            }
            if(j != mesh.jmax()) { // Check Down
                add = ((fabs(VOFarray[i][j+1]) > IRL::global_constants::VF_LOW) && (fabs(VOFarray[i][j]) <= 1e-6)) || add;
            }

            if(add) {
                // Get Square
                double xmin = mesh.x(i);
                double xmax = mesh.x(i+1);
                double ymin = mesh.y(i);
                double ymax = mesh.y(i+1);
                
                // Make Square
                std::vector<std::vector<double>> squ = {{xmin,ymin},{xmax,ymin}
                                                        ,{xmax,ymax},{xmin,ymax},{xmin,ymin}};

                // Make Struct
                emptyConstraintSet.push_back(squ);
            }
            // addArray[i][j] = add; // For Debugging
            
            // if add, then add to contraints
        }
    }
    my_empty_data emp = {emptyConstraintSet,1e-6};
    // Before Making Periodic, Change Order Such that we order them by the closest Points.
    std::vector<std::vector<double>> copyInterp = interpPoints;
    std::vector<int> traverseOrder = {0}; // Makes sure that multiple points cannot go back to the same point
    for(int i = 0; i < interpPoints.size()-1; i++) {
        // std::cout << " i = " << i << "\n";

        std::vector<double> Qcurr = copyInterp[traverseOrder[i]];
        double dmin = 10.0*mesh.dx();
        int minIndex = 0;
        for(int j = 0; j < interpPoints.size(); j++) {
            //  std::cout << " j = " << j << "\n";
            if(std::find(traverseOrder.begin(), traverseOrder.end(), j) != traverseOrder.end()) { // If in traverseOrder, do nothing
                continue;
            } else { // Have not used this point yet, check order
                std::vector<double> Qcheck = copyInterp[j];
                double dist = std::sqrt((Qcheck[0]-Qcurr[0])*(Qcheck[0]-Qcurr[0]) 
                                     + (Qcheck[1]-Qcurr[1])*(Qcheck[1]-Qcurr[1]));
                if(dist < dmin) { // If new minimum, set it
                    dmin = dist;
                    minIndex = j;
                }
            }
        }
        traverseOrder.push_back(minIndex);
    }
    std::cout << "======= Traverse Result ======\n";
    for(int i = 0; i < traverseOrder.size(); i++) { // Note that i corresponds to x, j corresponds to y
        std::cout << traverseOrder[i] << ",";
    }
    // Make Periodic
    traverseOrder.push_back(traverseOrder[0]);
    // std::cout << "Num Constraints = " << VOFconstraintSet.size() << "\n";
    // Print Result
    // std::cout << "========== Interp INITIAL =========\n";
    // for(int i = 0; i < interpPoints.size(); i++) { // Note that i corresponds to x, j corresponds to y
    //     std::cout << interpPoints[i][0] << "," << interpPoints[i][1] << "\n";
    // }
    // Interp Post Traverse
    copyInterp = interpPoints;
    copyInterp.push_back({0,0});
    for(int i = 0; i < traverseOrder.size(); i++) {
        copyInterp[i] = {interpPoints[traverseOrder[i]]};
    }
    interpPoints = copyInterp;
    // std::cout << "========== Interp Final =========\n";
    // for(int i = 0; i < interpPoints.size(); i++) { // Note that i corresponds to x, j corresponds to y
    //     std::cout << interpPoints[i][0] << "," << interpPoints[i][1] << "\n";
    // }
    // Calculate Tangents

    double tresh = 0.5;
    
    MarchingSquares<double> testMS(&mesh,tresh,2);
    std::vector<std::vector<double>> ret = testMS.vertexPoints(VOFarray); // Note this does not have the repeating start point
    interpPoints = ret;

    std::cout << "GOT MARCHING CUBES RESULT\n";
    std::vector<std::vector<double>> tangents = {};
    for(int i = 0; i < interpPoints.size(); i++) {
        int iNext = (i+1);
        int iPrev = (i-1);
        while(iNext >= interpPoints.size()) {
            iNext -= interpPoints.size();
        }
        while(iPrev < 0) {
            iPrev += interpPoints.size();
        }
        // Tangent points from previous point to next point 
        std::vector<double> Pprev = interpPoints[iPrev];
        std::vector<double> Pcurr = interpPoints[i];
        std::vector<double> Pnext = interpPoints[iNext];

        std::vector<double> tCurr = {Pnext[0]-Pprev[0],Pnext[1]-Pprev[1]};
        //Normalize
        double tCurrMag = std::sqrt(tCurr[0]*tCurr[0]+tCurr[1]*tCurr[1]);
        tCurr[0] = tCurr[0]/tCurrMag;
        tCurr[1] = tCurr[1]/tCurrMag;
        tangents.push_back(tCurr);
    }

    std::vector<std::vector<double>> ret2 = {ret[0]};
    std::vector<std::vector<double>> tangents2 = {tangents[0]};
    interpPoints = ret;
    Spline first = Spline<double>::LocalRQuadInterp(interpPoints,tangents);
    first.saveToVTK("MarchingCubesInitialGuess");

    // std::vector<std::vector<double>> Blob = {{1.5,-0.1},{1.0,1.2},{-1.5,-0.1},{1.0,-0.5},{1.5,-0.1}};
    // std::vector<std::vector<double>> BlobT = {{0.0,1.0},{-1.0,0.0},{0.0,-1},{1.0,0},{0.0,1.0}};
    // interpPoints = Blob;

    // tangs = BlobT;
    // Now we have an initial guess, we can go onto the optimization
    std::vector<double> xInitial = Spline<double>::pack(interpPoints,tangents);
    int variables = xInitial.size();
    std::cout << "\nVariables = " << variables <<"\n";
    std::cout << "Constraints = " << VOFconstraintSet.size() << "\n";
    // variables = 3*(interpPoints.size()-1);
    std::cout << "\nMake Optimizer\n";
    nlopt::opt splineOptim(nlopt::LD_MMA,variables);
    std::cout << "Made Optimizer\n";
    double initialArea = first.getArea();
    my_constraint_data InitialData = {initialArea};
    
    splineOptim.set_min_objective(myfunc,&InitialData);
    std::vector<double> stepSizes(variables);
    const std::vector<double> xTemp = xInitial;
    splineOptim.get_initial_step(xInitial,stepSizes);

    // std::cout << "====== Default Step Sizes =======\n Step Sizes = [";
    // for(int i = 0; i < stepSizes.size(); i++) {
    //     std::cout << stepSizes[i] << ",";
    // }
    // std::cout << "]\n";

    // for(int i = 0; i < stepSizes.size(); i+=3) {
    //     stepSizes[i] = mesh.dx()/4-(1e-6);
    //     stepSizes[i+1] = mesh.dy()/4-(1e-6);
    //     stepSizes[i+2] = M_PI/4;
    // }

    // for(int i = 0; i < stepSizes.size(); i+=3) {
    //     stepSizes[i] = 0.5;
    //     stepSizes[i+1] = 0.5;
    //     stepSizes[i+2] = M_PI/4;
    // }
    // std::cout << "stepSizes Size = " << stepSizes.size() << "\n";
    // splineOptim.add_equality_constraint(myconstraint,&InitialData,1e-12);
    // splineOptim.set_initial_step(stepSizes);
    // Set xtol
    // opt.set_ftol_rel(1e-12);
    splineOptim.set_ftol_abs(1e-4);
    // opt.set_xtol_abs(1e-16);
    splineOptim.set_xtol_rel(1e-4);
    // opt.set_maxeval(100000);
    splineOptim.set_maxtime(100); 
    // Constraint Data
    for(int i = 0; i < VOFconstraintSet.size();i++) {
        my_square_data tem = VOFconstraintSet[i];
        splineOptim.add_inequality_constraint(myVOFconstraint,&(VOFconstraintSet[i]),1e-3);
        // std::cout << "i = " << i+1 << "\n";
    }    
    // std::cout << "num Variables = " << variables << "\n";
    double minres;
    std::cout << "Start Optimizer\n";
    try {
        std::cout << "Start Optimizer2\n";
        nlopt::result res = splineOptim.optimize(xInitial,minres);
        std::cout << "Result = " << res << "\n";
        std::cout << "Found Minimum at f(" << xInitial[0] << "," << xInitial[1] << ") = " 
            << std::setprecision(10) << minres <<std::endl;
        
        // Output Spline
        std::vector<std::vector<std::vector<double>>> ret = Spline<double>::unpack(xInitial);
        std::vector<std::vector<double>> Q = ret[0];
        std::vector<std::vector<double>> T = ret[1];
        for(int i = 0; i < Q.size(); i++) {
            std::cout << "Q = (" << Q[i][0] << "," <<Q[i][1] << ") , ";
            std::cout << "T = (" << T[i][0] << "," <<T[i][1] << ")\n";
        }
        Spline out = Spline<double>::LocalRQuadInterp(Q,T);
        
        out.saveToVTK("AutomatedOptimizationInitial");

    } catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    std::cout << "End Optimizer\n";
// End WIP.


    // // Begin Optimization (Not Updated)
    // int Nd = 12;
    // std::cout << "Start\n";
    // nlopt::opt opt(nlopt::LN_COBYLA,Nd); // Initialize Optimizer, with algorithm and Dimensionality
    
    // // Initial Guess
    // // std::vector<std::vector<double>> Circle = {{1.000,0.0},{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0}};
    // // std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    // // std::vector<std::vector<double>> InterpTester = {{1.000,0.0},{0.0,1.000},{-1.000,0.0},{-3.000,0.0},{1.000,0.0}};
    // // std::vector<std::vector<double>> InterpTesterT = {{0.0,1.000},{-1.000,0.0},{-1.000,0.0},{-1.000,0.0},{0.0,1.000}};

    // std::vector<std::vector<double>> Blob = {{1.5,-0.1},{1.0,1.2},{-1.5,-0.1},{1.0,-0.5},{1.5,-0.1}};
    // std::vector<std::vector<double>> BlobT = {{0.0,1.0},{-1.0,0.0},{0.0,-1},{1.0,0},{0.0,1.0}};

    // // std::cout << "Interpolation\n";
    // Spline initial = Spline<double>::LocalRQuadInterp(Blob,BlobT);
    // // std::cout << "Interpolated\n";
    // std::vector<double> x = Spline<double>::pack(Blob,BlobT);
    // // std::cout << "Size x: " << x.size() << std::endl;
    
    // // Total Area Constraint Data
    // double Ainitial = initial.getArea();
    // // std::cout << "Aread\n";
    // my_constraint_data data = {Ainitial};

    // // VOF Constraint Data
    // std::vector<std::vector<double>> square = {{0.0,0.0},{1.1,0.0},{1.1,1.1},{0.0,1.1},{0.0,0.0}};
    // std::vector<std::vector<double>> square2 = {{0.0,0.0},{0.0,1.1},{-1.1,1.1},{-1.1,0.0},{0.0,0.0}};
    // std::vector<std::vector<double>> square3 = {{0.0,0.0},{-1.1,0.0},{-1.1,-1.1},{0.0,-1.1},{0.0,0.0}};
    // std::vector<std::vector<double>> square4 = {{0.0,0.0},{0.0,-1.1},{1.1,-1.1},{1.1,0.0},{0.0,0.0}};

    // my_square_data dat = {square,0.5};
    // my_square_data dat2 = {square2,0.2};
    // my_square_data dat3 = {square3,0.2};
    // my_square_data dat4 = {square4,0.8};

    // // Bound on objective function
    // opt.set_min_objective(myfunc,&data);
    // std::vector<double> steps(Nd);
    // for(int i = 0; i < steps.size(); i+=3) {
    //     steps[i] = 0.5;
    //     steps[i+1] = 0.5;
    //     steps[i+2] = M_PI/4;
    // }
    // opt.set_initial_step(steps);
    // // Set xtol
    // // opt.set_ftol_rel(1e-12);
    // opt.set_ftol_abs(1e-16);
    // // opt.set_xtol_abs(1e-16);
    // opt.set_xtol_rel(1e-16);
    // // opt.set_maxeval(100000);
    // opt.set_maxtime(10);
    // // Constraint Data
    // opt.add_equality_constraint(myconstraint,&data,1e-12);
    // // opt.add_inequality_constraint(myconstraint,&data,1e-12);
    // // opt.add_equality_constraint(myVOFconstraint,&dat,1e-6);
    // opt.add_equality_constraint(myVOFconstraint,&dat2,1e-6);
    // // opt.add_equality_constraint(myVOFconstraint,&dat3,1e-6);
    // // opt.add_equality_constraint(myVOFconstraint,&dat4,1e-6);

    // initial.saveToVTK("StartingPoint");
    // double minf;
    // std::cout << "x = [";
    // for(int i = 0; i < x.size(); i++) {
    //     std::cout << x[i] << ",";
    // }
    // std::cout <<"]\n";
    // try {
    //     nlopt::result result = opt.optimize(x,minf);
    //     std::cout << "Result = " << result << "\n";
    //     std::cout << "Found Minimum at f(" << x[0] << "," << x[1] << ") = " 
    //         << std::setprecision(10) << minf <<std::endl;
        
    //     std::cout << "x = [";
    //     for(int i = 0; i < x.size(); i++) {
    //         std::cout << x[i] << ",";
    //     }
    //     std::cout <<"]\n";
    //     // Output Spline
    //     std::vector<std::vector<std::vector<double>>> ret = Spline<double>::unpack(x);
    //     std::vector<std::vector<double>> Q = ret[0];
    //     std::vector<std::vector<double>> T = ret[1];
    //     // for(int i = 0; i < Q.size(); i++) {
    //     //     std::cout << "Q = (" << Q[i][0] << "," <<Q[i][1] << ") , ";
    //     //     std::cout << "T = (" << T[i][0] << "," <<T[i][1] << ")\n";
    //     // }
    //     Spline out = Spline<double>::LocalRQuadInterp(Q,T);


    //     std::vector<double> xFinal = 
    //         {0.678274,-0.31701,2.11799,
    //             -0.0897338,0.383153,2.84261,
    //             -1.47507,0.607397,3.41242,
    //             -2.11324,0.327595,3.72911};
    //     // Other Spline
    //     ret = Spline<double>::unpack(xFinal);
    //     Q = ret[0];
    //     T = ret[1];
    //     Spline out2 = Spline<double>::LocalRQuadInterp(Q,T);
    //     out2.saveToVTK("PrintedFinalX");    
    //     initial.printControlPoints();
    //     out.printControlPoints();
    //     out.saveToVTK("DebuggerCircle");
    //     std::cout << "\nInitial Area " << initial.getArea();
    //     std::cout << "\nFinal Area " << out.getArea();
        
    // } catch(std::exception &e) {
    //     std::cout << "nlopt failed: " << e.what() << std::endl;
    // }

    return 1;
    }   