#include <iomanip>
#include <iostream>
#include <vector>
#include <math.h>
#include <nlopt.hpp>
#include "examples/spline_optimization/Spline.h"
// #include "irl/splines/Spline.h"
#include "examples/spline_optimization/basic_mesh.h"
#include "examples/spline_optimization/irl2d.h"
#include "examples/spline_optimization/data.h"



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

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
    std::cout << "Call my Func\n";
    my_constraint_data *d = (my_constraint_data *) my_func_data;
    double A = d->A;
    // std::cout <<"Grad Size" << grad.size() << "\n";
    if(!grad.empty()) {
        grad[0] = 1.0;
        grad[1] = 1.0;
    }
    // std::cout << "Unpacking\n";
    std::vector<std::vector<std::vector<double>>> ret = Spline::unpack(x);
    // std::cout << "Interpolating \n";
    Spline s = Spline::LocalRQuadInterp(ret[0],ret[1]);
    // std::cout << "Interpolated \n";
    // std::cout << "Getting Surface Energy \n";
    double E = s.getSurfaceEnergy();
    // std::cout << "Got Surface Energy \n";
    // std::cout << "Getting Arc Length \n";
    double Al = s.getArcLength();
    // std::cout << "Got Arc Length \n";
    double AlFactor = 2*sqrt(M_PI*A);
    std::cout << "E ===" << E/(2*M_PI) << "\n";
    std::cout << "Al ===" << Al/AlFactor << "\n";

    std::cout << "======" << E/(2*M_PI) + Al/AlFactor << "\n";
    std::cout << "End my Func\n";
    return E/(2*M_PI) + Al/AlFactor;
}

// This is the constraint function
double myVOFconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data) { // VOF Constraint
    std::cout << "Call VOF Constraint \n";
    my_square_data *d = (my_square_data *) data;
    std::vector<std::vector<double>> square = d->square; // 5 Point Rectangle
    double VOF = d->VOF;
    // std::cout << "dx, dy get\n";
    double dx = fabs(square[0][0]-square[2][0]);
    // std::cout << "dx done\n";
    double dy = fabs(square[0][1]-square[2][1]);
    // std::cout << "dy done\n";
    if(!grad.empty()) {
        grad[0] = -2*x[0];
        grad[1] = -2*x[1];
    }
    std::vector<std::vector<std::vector<double>>> ret = Spline::unpack(x);
    std::vector<std::vector<double>> Q = ret[0];
    std::vector<std::vector<double>> T = ret[1];
    
    Spline s = Spline::LocalRQuadInterp(Q,T);
    std::cout << "Finding Intersection\n";
    // s.saveToVTK("OptimizationSplineTest");
    double Acurr = s.integrateSplineSquare(square);
    std::cout << "Found\n";
    std::cout << "Done\n";
    std::cout <<"===========================" <<fabs(Acurr/(dx*dy) - VOF) << "\n";
    std::cout << "End VOF Constraint \n";

    return fabs(Acurr/(dx*dy) - VOF); // Calculated VOF - Exact VOF
}

double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data) { // Total Area
    std::cout << "Call Constraint \n";
    my_constraint_data *d = (my_constraint_data *) data;
    double A = d->A;
    

    std::vector<std::vector<std::vector<double>>> ret = Spline::unpack(x);
    std::vector<std::vector<double>> Q = ret[0];
    std::vector<std::vector<double>> T = ret[1];
    Spline s = Spline::LocalRQuadInterp(Q,T);
    double Acurr = s.getArea();
    if(!grad.empty()) {
        Q = ret[0];
        T = ret[1];
        for(int i = 0; i < x.size(); i++) {
                
        }

        grad[0] = -2*x[0];
        grad[1] = -2*x[1];
    }
    
    // std::cout << " ============= \n";
    // for(int i = 0; i < Q.size(); i++) {
    //     std::cout << "( " <<Q[i][0] << "," << Q[i][1] << ") at ";
    //     std::cout << "[ " <<T[i][0] << "," << T[i][1] << "]\n";
    // }
    // std::cout << "\nx size = " << x.size();
    // std::cout << "\nQ size = " << Q.size();
    // std::cout << "\nT size = " << T.size();
    // std::cout << "Interpolating \n";
    
    // std::cout << "Interpolated \n";
    // std::cout << "Getting Area\n";
    
    // std::cout << "Area Got\n";
    // std::cout << fabs(Acurr-A) << "\n";
    std::cout << "End Constraint \n";
    return (Acurr-A);
}


// ================================== BEGIN MAIN ===========================================================
int main() {

// Ignore this, WIP
    // Field Initialization
//     const int a_nx = 10;
//     const int GC = 3;
//     // Set up Mesh
//     BasicMesh mesh(a_nx, a_nx, GC);
//     IRL2D::Vec my_lower_domain(-0.5, -0.5);
//     IRL2D::Vec my_upper_domain(0.5, 0.5);
//     mesh.setCellBoundaries(my_lower_domain, my_upper_domain);

//     // Initialize
//     const auto circle_center = IRL2D::Vec(0,0);
//     const double circle_radius = 0.15;

//     Data<IRL2D::Parabola> a_interface;
//     // a_interface->updateBorder();
//     // correctInterfaceBorders(a_interface);
//   // Loop over cells in domain. Skip if cell is not mixed phase.
//   for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
//     for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
//       const IRL2D::Vec lower_cell_pt(mesh.x(i), mesh.y(j));
//       const IRL2D::Vec upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1));
//       const IRL2D::Vec mid_pt = 0.5 * (lower_cell_pt + upper_cell_pt);
//       IRL2D::Vec disp = mid_pt - circle_center;
//       const auto mag = disp.magnitude();
//       if (mag < circle_radius - 2.0 * mesh.dx()) {
//         (a_interface)(i, j).markAsAlwaysAbove();
//       } else if (mag > circle_radius + 2.0 * mesh.dx()) {
//         (a_interface)(i, j).markAsAlwaysBelow();
//       } else {
//         auto circle_normal = IRL2D::Vec(disp);
//         circle_normal.normalize();
//         const IRL2D::Vec datum = circle_center + circle_radius * circle_normal;
//         const IRL2D::ReferenceFrame frame(
//             IRL2D::Vec(circle_normal[1], -circle_normal[0]), circle_normal);
//         const double coeff = 1.0 / (2.0 * circle_radius);
//         (a_interface)(i, j) = IRL2D::Parabola(datum, frame, coeff);
//         for (int ii = 0; ii < 10; ii++) {
//           auto cell = IRL2D::RectangleFromBounds(lower_cell_pt, upper_cell_pt);
//           auto moments = IRL2D::ComputeMoments(cell, (a_interface)(i, j));
//           const double liquid_volume_fraction =
//               moments.m0() / IRL2D::ComputeArea(cell);
//           if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
//               liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
//             const auto new_disp = moments.m1() / moments.m0() - circle_center;
//             auto new_circle_normal = IRL2D::Vec(new_disp);
//             new_circle_normal.normalize();
//             const IRL2D::Vec new_datum =
//                 circle_center + circle_radius * new_circle_normal;
//             const IRL2D::ReferenceFrame new_frame(
//                 IRL2D::Vec(new_circle_normal[1], -new_circle_normal[0]),
//                 new_circle_normal);
//             (a_interface)(i, j) = IRL2D::Parabola(new_datum, new_frame, coeff);
//           } else {
//             break;
//           }
//         }
//         }
//       }
//     }
// End WIP.


    


















    // Begin Optimization (Not Updated)
    int Nd = 12;
    std::cout << "Start\n";
    nlopt::opt opt(nlopt::LN_COBYLA,Nd); // Initialize Optimizer, with algorithm and Dimensionality
    
    // Initial Guess
    std::vector<std::vector<double>> Circle = {{1.0001,0.0},{0.0,1.01},{-0.9999,0.0},{0.0,-1.0},{1.0001,0.0}};
    std::vector<std::vector<double>> CircleT = {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};

    std::vector<std::vector<double>> InterpTester = {{1.000,0.0},{0.0,1.000},{-1.000,0.0},{-3.000,0.0},{1.000,0.0}};
    std::vector<std::vector<double>> InterpTesterT = {{0.0,1.000},{-1.000,0.0},{-1.000,0.0},{-1.000,0.0},{0.0,1.000}};

    std::vector<std::vector<double>> Blob = {{1.5,-0.1},{1.0,1.2},{-1.5,-0.1},{1.0,-0.5},{1.5,-0.1}};
    std::vector<std::vector<double>> BlobT = {{0.0,1.0},{-1.0,0.0},{0.0,-1},{1.0,0},{0.0,1.0}};

    std::cout << "Interpolation\n";
    Spline initial = Spline::LocalRQuadInterp(InterpTester,InterpTesterT);
    std::cout << "Interpolated\n";
    std::vector<double> x = Spline::pack(InterpTester,InterpTesterT);
    std::cout << "Size x: " << x.size() << std::endl;
    
    // Total Area Constraint Data
    double Ainitial = initial.getArea();
    std::cout << "Aread\n";
    my_constraint_data data = {Ainitial};

    // VOF Constraint Data
    std::vector<std::vector<double>> square = {{0.0,0.0},{1.1,0.0},{1.1,1.1},{0.0,1.1},{0.0,0.0}};
    std::vector<std::vector<double>> square2 = {{0.0,0.0},{0.0,1.1},{-1.1,1.1},{-1.1,0.0},{0.0,0.0}};
    std::vector<std::vector<double>> square3 = {{0.0,0.0},{-1.1,0.0},{-1.1,-1.1},{0.0,-1.1},{0.0,0.0}};
    std::vector<std::vector<double>> square4 = {{0.0,0.0},{0.0,-1.1},{1.1,-1.1},{1.1,0.0},{0.0,0.0}};

    my_square_data dat = {square,0.5};
    my_square_data dat2 = {square2,0.2};
    my_square_data dat3 = {square3,0.2};
    my_square_data dat4 = {square4,0.8};

    // Bound on objective function
    opt.set_min_objective(myfunc,&data);
    
    // Set xtol
    opt.set_ftol_rel(1e-6);
    opt.set_maxtime(10);
    
    // Constraint Data
    opt.add_equality_constraint(myconstraint,&data,1e-6);
    // opt.add_equality_constraint(myVOFconstraint,&dat,1e-6);
    // opt.add_equality_constraint(myVOFconstraint,&dat2,1e-6);
    // opt.add_equality_constraint(myVOFconstraint,&dat3,1e-6);
    // opt.add_equality_constraint(myVOFconstraint,&dat4,1e-6);

    initial.saveToVTK("StartingPoint");
    double minf;
    std::cout << "x = [";
    for(int i = 0; i < x.size(); i++) {
        std::cout << x[i] << ",";
    }
    std::cout <<"]\n";
    try {
        nlopt::result result = opt.optimize(x,minf);
        std::cout << "Found Minimum at f(" << x[0] << "," << x[1] << ") = " 
            << std::setprecision(10) << minf <<std::endl;
        for(int i = 0; i < x.size(); i+=3) {
            std::cout << x[i] << "," << x[i+1] << "," << x[i+2] << "\n";
        }
        std::cout << "x = [";
        for(int i = 0; i < x.size(); i++) {
            std::cout << x[i] << ",";
        }
        std::cout <<"]\n";
        // Output Spline
        std::vector<std::vector<std::vector<double>>> ret = Spline::unpack(x);
        std::vector<std::vector<double>> Q = ret[0];
        std::vector<std::vector<double>> T = ret[1];
        for(int i = 0; i < Q.size(); i++) {
            std::cout << "Q = (" << Q[i][0] << "," <<Q[i][1] << ") , ";
            std::cout << "T = (" << T[i][0] << "," <<T[i][1] << ")\n";
        }
        Spline out = Spline::LocalRQuadInterp(Q,T);
        
        initial.printControlPoints();
        out.printControlPoints();
        out.saveToVTK("DebuggerCircle");
        std::cout << "\nInitial Area " << initial.getArea();
        std::cout << "\nFinal Area " << out.getArea();
        
    } catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }

    return 1;
    }   