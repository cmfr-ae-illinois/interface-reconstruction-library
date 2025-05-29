#include <array>
#include <iostream>
#include <math.h>
#include <fstream>

#include "examples/spline_optimization/basic_mesh.h"
#include "examples/spline_optimization/irl2d.h"
#include "examples/spline_optimization/data.h"

int main() {
    // Field Initialization
    const int a_nx = 6;
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
    std::cout << "========== VOF INITIAL =========\n";
    for(int i = 0; i < a_nx; i++) { // Note that i corresponds to x, j corresponds to y
        for(int j = 0; j < a_nx; j++) {
            std::cout.width(15); std::cout << std::left << std::setprecision(3) << VOFarray[i][j];
        }
        std::cout << "\n";   
    }
    std::vector<std::vector<double>> xValues(a_nx, std::vector<double>(a_nx,false));
    std::vector<std::vector<double>> yValues(a_nx, std::vector<double>(a_nx,false));
    std::cout << "========== x INITIAL =========\n";
    for(int i = 0; i < a_nx; i++) { // Note that i corresponds to x, j corresponds to y
        for(int j = 0; j < a_nx; j++) {
            xValues[i][j] = mesh.xm(i);
            yValues[i][j] = mesh.ym(j);
            std::cout.width(10); std::cout << std::left << std::setprecision(3) << xValues[i][j];
        }
        std::cout << "\n";
    }
    std::cout << "========== y INITIAL =========\n";
    for(int i = 0; i < a_nx; i++) { // Note that i corresponds to x, j corresponds to y
        for(int j = 0; j < a_nx; j++) {
            std::cout.width(10); std::cout << std::left << std::setprecision(3) << yValues[i][j];
        }
        std::cout << "\n";
    }

    // Treshold
    double tresh = 0.2;
    // First, translate VOF to be in corners. This is done by taking an average of the 4 cells touching a corner.
    // If on a boundary, use an average of the two adjacent cells intead. 
    // If a corner, use just the cell value;
    std::cout << "========== Vertex INITIAL =========\n";
    std::vector<std::vector<double>> vertexVOF(a_nx+1,std::vector<double>(a_nx+1,0));
    for(int i=0;i<a_nx+1;i++){
        for(int j=0;j<a_nx+1;j++) {
            double sum = 0;
            double count = 0;
            if(i != 0) { 
                if(j != a_nx) {
                    sum += VOFarray[i-1][j];
                    count++;
                }
                if(j != 0) {
                    sum += VOFarray[i-1][j-1];
                    count++;
                }
            }
            if(i != a_nx) {
                if(j != a_nx) {
                    sum += VOFarray[i][j];
                    count++;
                }
                if(j != 0) {
                    sum += VOFarray[i][j-1];
                    count++;
                }
            }
            // Average
            vertexVOF[i][j] = sum/count;
            std::cout.width(15); std::cout << std::left << std::setprecision(3) << vertexVOF[i][j];
        }
        std::cout << "\n";
    }
    // Apply Tresholding to get Binary Value
    std::vector<std::vector<bool>> meetsTreshold(a_nx+1, std::vector<bool>(a_nx+1,false));
    for(int i = 0; i < a_nx+1; i++) { // Note that i corresponds to x, j corresponds to y
        for(int j = 0; j < a_nx+1; j++) {
            meetsTreshold[i][j] = vertexVOF[i][j] >= tresh;
        }
    }
    std::cout << "========== Tresh INITIAL =========\n";
    for(int i = 0; i < meetsTreshold.size(); i++) { // Note that i corresponds to x, j corresponds to y
        for(int j = 0; j < meetsTreshold.size(); j++) {
            std::cout.width(5); std::cout << std::left << std::setprecision(1) << meetsTreshold[i][j];
        }
        std::cout << "\n";
    }

    std::cout << "========== Case Values =========\n";
    std::ofstream myfile;
    myfile.open("marching_cubes_output/caseOutput.txt");
    std::vector<std::vector<int>> cases(a_nx,std::vector<int>(a_nx,-1));
    for(int i = 0; i < cases.size(); i++){
        for(int j = 0; j < cases.size();j++) {
            cases[i][j] = meetsTreshold[i][j] + 2*meetsTreshold[i+1][j]+4*meetsTreshold[i+1][j+1]+8*meetsTreshold[i][j+1];
            std::cout.width(5); std::cout << std::left << std::setprecision(1) << cases[i][j];
            myfile << cases[i][j];
            if(j != cases.size()-1) {
                myfile <<",";
            } 
        }
        myfile << "\n";
        std::cout << "\n";
    }
    myfile.close();

    std::cout << "========== Edge Values =========\n";
    // Here we make a cvs that contains the edges in the form (x1,y1,x2,y2). This is to make plotting easier.
    std::vector<std::vector<double>> edges = {}; // This contains sets of points for edge vectors.
    // That is to say it is [[P1x,P1y,P2x,P2y],...]
    double xin,yin,xout,yout,xmid,ymid;
    for(int i = 0; i < a_nx-1; i++) {
        for(int j =0; j < a_nx-1; j++) {
            switch(cases[i][j]) {
                case 0:
                    break; // No Lines
                case 1: // One in corner
                    // // In Corners
                    xin = mesh.x(i);
                    yin = mesh.y(j);

                    // Out Corners
                    xout = mesh.x(i);
                    yout = mesh.y(j+1);

                    // Midpoints
                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    // // Connect to make a line
                    edges.push_back({xmid,yin,xin,ymid});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xmid << ",";
                    std::cout << yin << ",";
                    std::cout << xin << ",";
                    std::cout << ymid << "\n";
                    break;
                case 2:
                    // // In Corners
                    xin = mesh.x(i+1);
                    yin = mesh.y(j);
                    
                    // Out Corners
                    xout = mesh.x(i);
                    yout = mesh.y(j+1);

                    // Midpoints
                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    // // Connect to make a line
                    edges.push_back({xin,ymid,xmid,yin});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xin << ",";
                    std::cout << ymid << ",";
                    std::cout << xmid << ",";
                    std::cout << yin << "\n";
                    break;
                case 3:
                    yin = mesh.y(j);
                    yout = mesh.y(j+1);

                    xin = mesh.x(i); // Both x values are in, but I don't want to make new variables
                    xout = mesh.x(i+1);

                    ymid = (yin+yout)/2;
                    edges.push_back({xout,ymid,xin,ymid});
                    std::cout << xout << ",";
                    std::cout << ymid << ",";
                    std::cout << xin << ",";
                    std::cout << ymid << "\n";
                    break;
                case 4:
                    // // In Corners
                    xin = mesh.x(i+1);
                    yin = mesh.y(j+1);
                    
                    // Out Corners
                    xout = mesh.x(i);
                    yout = mesh.y(j);

                    // Midpoints
                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    // // Connect to make a line
                    edges.push_back({xmid,yin,xin,ymid});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xmid << ",";
                    std::cout << yin << ",";
                    std::cout << xin << ",";
                    std::cout << ymid << "\n";
                    break;
                case 5:
                    xin = mesh.x(i); 
                    yin = mesh.y(j);

                    xout = mesh.x(i+1);
                    yout = mesh.y(j+1); 
                    // Both Corners are in, but once again I did not want to make more variables
                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    // Make First Edge
                    edges.push_back({xmid,yin,xout,ymid});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xmid << ",";
                    std::cout << yin << ",";
                    std::cout << xout << ",";
                    std::cout << ymid << "\n";

                    // Make Second Edge
                    edges.push_back({xmid,yout,xin,ymid});
                    std::cout << xmid << ",";
                    std::cout << yout << ",";
                    std::cout << xin << ",";
                    std::cout << ymid << "\n";
                    break;
                case 6:
                    xin = mesh.x(i+1);
                    xout = mesh.x(i);
                    // In this case, both y values are in
                    yin = mesh.y(j);
                    yout = mesh.y(j+1);

                    xmid = (xin+xout)/2;
                    edges.push_back({xmid,yout,xmid,yin});
                    std::cout << xmid << ",";
                    std::cout << yout << ",";
                    std::cout << xmid << ",";
                    std::cout << yin << "\n";
                    break;
                case 7: // This is case 8, but in and out are flipped
                    // // Out Corners
                    xout = mesh.x(i);
                    yout = mesh.y(j+1);
                    
                    // In Corners
                    xin = mesh.x(i+1);
                    yin = mesh.y(j);

                    // Midpoints
                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    // // Connect to make a line
                    edges.push_back({xmid,yout,xout,ymid});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xmid << ",";
                    std::cout << yout << ",";
                    std::cout << xout << ",";
                    std::cout << ymid << "\n";
                    break;
                case 8:
                    // // In Corners
                    xin = mesh.x(i);
                    yin = mesh.y(j+1);
                    
                    // Out Corners
                    xout = mesh.x(i+1);
                    yout = mesh.y(j);

                    // Midpoints
                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    // // Connect to make a line
                    edges.push_back({xin,ymid,xmid,yin});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xin << ",";
                    std::cout << ymid << ",";
                    std::cout << xmid << ",";
                    std::cout << yin << "\n";
                    break;
                case 9:
                    xin = mesh.x(i);
                    xout = mesh.x(i+1);
                    
                    yin = mesh.y(j);
                    yout = mesh.y(j+1);

                    xmid = (xin+xout)/2;

                    edges.push_back({xmid,yin,xmid,yout});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xmid << ",";
                    std::cout << yin << ",";
                    std::cout << xmid << ",";
                    std::cout << yout << "\n";
                    break;
                case 10:
                    xin = mesh.x(i+1);
                    xout = mesh.x(i);

                    yin = mesh.y(j+1);
                    yout = mesh.y(j);

                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    // First Edge
                    edges.push_back({xin,ymid,xmid,yin});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xin << ",";
                    std::cout << ymid << ",";
                    std::cout << xmid << ",";
                    std::cout << yin << "\n";

                    // Second Edge
                    edges.push_back({xout,ymid,xmid,yout});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xout << ",";
                    std::cout << ymid << ",";
                    std::cout << xmid << ",";
                    std::cout << yout << "\n";
                    break;
                case 11:
                    xout = mesh.x(i+1);
                    yout = mesh.y(j+1);

                    xin = mesh.x(i);
                    yin = mesh.y(j);

                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    edges.push_back({xout,ymid,xmid,yout});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xout << ",";
                    std::cout << ymid << ",";
                    std::cout << xmid << ",";
                    std::cout << yout << "\n";
                    break;
                case 12:
                    xin = mesh.x(i);
                    xout = mesh.x(i+1);

                    yin = mesh.y(j);
                    yout = mesh.y(j+1);

                    ymid = (yin+yout)/2;

                    edges.push_back({xin,ymid,xout,ymid});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xin << ",";
                    std::cout << ymid << ",";
                    std::cout << xout << ",";
                    std::cout << ymid << "\n";
                    break;
                case 13:
                    xout = mesh.x(i+1);
                    yout = mesh.y(j);

                    xin = mesh.x(i);
                    yin = mesh.y(j+1);

                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    edges.push_back({xmid,yout,xout,ymid});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xmid << ",";
                    std::cout << yout << ",";
                    std::cout << xout << ",";
                    std::cout << ymid << "\n";
                    break;
                case 14:
                    xout = mesh.x(i);
                    yout = mesh.y(j);

                    xin = mesh.x(i+1);
                    yin = mesh.y(j+1);

                    xmid = (xin+xout)/2;
                    ymid = (yin+yout)/2;

                    edges.push_back({xout,ymid,xmid,yout});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << xout << ",";
                    std::cout << ymid << ",";
                    std::cout << xmid << ",";
                    std::cout << yout << "\n";
                    break;
                case 15:
                    break;
                default:
                    break;
            }
        }
    }    

    // Loop over all quads
    return 0;
}