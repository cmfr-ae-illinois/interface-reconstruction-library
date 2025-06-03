#include <array>
#include <iostream>
#include <math.h>
#include <fstream>

#include "examples/spline_optimization/basic_mesh.h"
#include "examples/spline_optimization/irl2d.h"
#include "examples/spline_optimization/data.h"
#include "irl/splines/Spline.h"

int main() {
    // Field Initialization
    const int a_nx = 10;
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
    double tresh = 0.5;
    // First, translate VOF to be in corners. This is done by taking an average of the 4 cells touching a corner.
    // If on a boundary, use an average of the two adjacent cells intead. 
    // If a corner, use just the cell value;
    // std::cout << "========== Vertex INITIAL =========\n";
    // std::vector<std::vector<double>> vertexVOF(a_nx+1,std::vector<double>(a_nx+1,0));
    // for(int i=0;i<a_nx+1;i++){
    //     for(int j=0;j<a_nx+1;j++) {
    //         double sum = 0;
    //         double count = 0;
    //         if(i != 0) { 
    //             if(j != a_nx) {
    //                 sum += VOFarray[i-1][j];
    //                 count++;
    //             }
    //             if(j != 0) {
    //                 sum += VOFarray[i-1][j-1];
    //                 count++;
    //             }
    //         }
    //         if(i != a_nx) {
    //             if(j != a_nx) {
    //                 sum += VOFarray[i][j];
    //                 count++;
    //             }
    //             if(j != 0) {
    //                 sum += VOFarray[i][j-1];
    //                 count++;
    //             }
    //         }
    //         // Average
    //         vertexVOF[i][j] = sum/count;
    //         std::cout.width(15); std::cout << std::left << std::setprecision(3) << vertexVOF[i][j];
    //     }
    //     std::cout << "\n";
    // }
    // Apply Tresholding to get Binary Value
    std::vector<std::vector<bool>> meetsTreshold(a_nx, std::vector<bool>(a_nx,false));
    for(int i = 0; i < meetsTreshold.size(); i++) { // Note that i corresponds to x, j corresponds to y
        for(int j = 0; j < meetsTreshold.size(); j++) {
            meetsTreshold[i][j] = VOFarray[i][j] >= tresh;
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
    std::vector<std::vector<int>> cases(a_nx-1,std::vector<int>(a_nx-1,-1));
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
    for(int i = 0; i < cases.size(); i++) {
        for(int j =0; j < cases.size(); j++) {
            // Endpoints
            double xLeft = mesh.x(i);
            double xRight = mesh.x(i+1);
            double yBot = mesh.y(j);
            double yTop = mesh.y(j+1);
            // Calculate Edge Values
            std::vector<std::vector<double>> Points = {{0,0},{0,0},{0,0},{0,0}};
            std::vector<std::vector<double>> Shift = {{0,0},{1,0},{1,1},{0,1},{0,0}};
            for(int k = 0; k < 4; k++) {
                // Extra Gamma Values and Coordinates of endpoints
                double gamma0 = VOFarray[i+Shift[k][0]][j+Shift[k][1]];
                double gamma1 = VOFarray[i+Shift[k+1][0]][j+Shift[k+1][1]];
                double x0 = mesh.x(i+Shift[k][0]);
                double x1 = mesh.x(i+Shift[k+1][0]);
                double y0 = mesh.y(j+Shift[k][1]);
                double y1 = mesh.y(j+Shift[k+1][1]);

                // Make sure they are ordered properly (ie, gamma0 < gamma1)
                double a1 = 0;
                double a2 = 0;
                std::vector<double> P1 = {0,0};
                std::vector<double> P2 = {0,0};
                if(gamma0 < gamma1) {
                    a1 = gamma0;
                    a2 = gamma1;
                    P1 = {x0,y0};
                    P2 = {x1,y1};
                } else {
                    a1 = gamma1;
                    a2 = gamma0;
                    P1 = {x1,y1};
                    P2 = {x0,y0};
                }
                // Now that we know a1 < a2, we can check that a1 < 0.5, a2 > 0.5
                // If that is not true, then we do not have an edge that interesects that face, 
                // so we will just default to the midpoint
                double mu = 0.5;
                if(a1 > 0.5 || a2 < 0.5) {
                    Points[k][0] = (P1[0]+P2[0])/2;
                    Points[k][1] = (P1[1]+P2[1])/2;
                    mu = 0.5;
                } else { // Here we are within the 4 cases to check
                    if(a2 <= 3*a1 || 3*a2 <= a1+2) { // Either or both True
                        if(a2 <= 3*a1 && 3*a2 <= a1+2) { // Both True, Case 1
                            mu = (a1-0.5)/(a1-a2);
                        } else { // Only 1 true, check to see which case. Do case 2 by default
                            // Case 2 Default
                            mu = 1.5-a1-a2; 
                            // Case 3 Check
                            if((2*a1+2*a2-1)*(2*a1+2*a2-1) < 4*a1*(a1+a2)) {
                                mu = 1 - (2*a2-1)/(8*a1+4*a2-8*std::sqrt(a1*(a1+a2)));
                            }
                            // Case 4 Check
                            double a1b = 1-a2;
                            double a2b = 1-a1; 
                            double a1f = 1-a1;
                            double a2f = 1-a2;
                            if((2*a1b+2*a2b-1)*(2*a1b+2*a2b-1) < 4*a1b*(a1b+a2b)) {
                                // Manson Method
                                mu = (2-a2b-1)/(8*a1b+4*a2b-8*std::sqrt(a1b*(a1b+a2b)));
                                // Fabien Method
                                std::cout << " Fabien Method\n"; 
                                mu = 1 - (2*a2f-1)/(8*a1f+4*a2f-8*std::sqrt(a1f*(a1f+a2f)));
                            }
                        }
                    } else { // Neither True, Case 2
                        mu = 1.5-a1-a2; 
                    }
                }
                std::cout << std::setprecision(5);
                // Now that we have the interpolation Parameter mu, we can calculate the point on the edge
                Points[k][0] = P1[0]*(1-mu)+P2[0]*mu;
                Points[k][1] = P1[1]*(1-mu)+P2[1]*mu;
            }
            switch(cases[i][j]) {
                case 0:
                    break; // No Lines
                case 1: // One in corner
                    edges.push_back({Points[0][0],Points[0][1],Points[3][0],Points[3][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 2:
                    // edges.push_back(Points[1]);
                    // edges.push_back(Points[0]);
                    edges.push_back({Points[1][0],Points[1][1],Points[0][0],Points[0][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 3:
                    // edges.push_back(Points[1]);
                    // edges.push_back(Points[3]);
                    edges.push_back({Points[1][0],Points[1][1],Points[3][0],Points[3][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 4:
                    // edges.push_back(Points[2]);
                    // edges.push_back(Points[1]);
                    edges.push_back({Points[2][0],Points[2][1],Points[1][0],Points[1][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 5:
                    // edges.push_back(Points[0]);
                    // edges.push_back(Points[1]);
                    edges.push_back({Points[0][0],Points[0][1],Points[1][0],Points[1][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";

                    // edges.push_back(Points[2]);
                    // edges.push_back(Points[3]);
                    edges.push_back({Points[2][0],Points[2][1],Points[3][0],Points[3][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 6:
                    // edges.push_back(Points[2]);
                    // edges.push_back(Points[0]);
                    edges.push_back({Points[2][0],Points[2][1],Points[0][0],Points[0][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 7: 
                    // edges.push_back(Points[2]);
                    // edges.push_back(Points[3]);
                    edges.push_back({Points[2][0],Points[2][1],Points[3][0],Points[3][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 8:
                    // edges.push_back(Points[3]);
                    // edges.push_back(Points[2]);
                    edges.push_back({Points[3][0],Points[3][1],Points[2][0],Points[2][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 9:
                    // edges.push_back(Points[0]);
                    // edges.push_back(Points[2]);
                    edges.push_back({Points[0][0],Points[0][1],Points[2][0],Points[2][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 10:
                    // First Edge
                    // edges.push_back(Points[1]);
                    // edges.push_back(Points[2]);
                    edges.push_back({Points[1][0],Points[1][1],Points[2][0],Points[2][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    // Second Edge
                    // edges.push_back(Points[3]);
                    // edges.push_back(Points[0]);
                    edges.push_back({Points[3][0],Points[3][1],Points[0][0],Points[0][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 11:
                    // edges.push_back(Points[1]);
                    // edges.push_back(Points[2]);
                    edges.push_back({Points[1][0],Points[1][1],Points[2][0],Points[2][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 12:
                    // edges.push_back(Points[3]);
                    // edges.push_back(Points[1]);
                    edges.push_back({Points[3][0],Points[3][1],Points[1][0],Points[1][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 13:
                    // edges.push_back(Points[0]);
                    // edges.push_back(Points[1]);
                    edges.push_back({Points[0][0],Points[0][1],Points[1][0],Points[1][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 14:
                    // edges.push_back(Points[3]);
                    // edges.push_back(Points[0]);
                    edges.push_back({Points[3][0],Points[3][1],Points[0][0],Points[0][1]});
                    std::cout << "Case " << cases[i][j] << " Result: ";
                    std::cout << edges[edges.size()-1][0] << ",";
                    std::cout << edges[edges.size()-1][1] << ",";
                    std::cout << edges[edges.size()-1][2] << ",";
                    std::cout << edges[edges.size()-1][3] << "\n";
                    break;
                case 15:
                    break;
                default:
                    break;
            }
        }
    }    

    // Loop over all quads and add everything to our file
    myfile.open("marching_cubes_output/edgeOutput-manson.txt");
    // First, insert grid characteristics
    myfile << xmin << ",";
    myfile << xmax << ",";
    myfile << ymin << ",";
    myfile << ymax << ",";
    myfile << a_nx << "\n";
    for(int i = 0; i < edges.size(); i++) {
        for(int j = 0; j < edges[0].size(); j++) {
            myfile << edges[i][j] << ",";
        }
        myfile << "\n";
    }
    myfile.close();

    // Now, make the connections ordered
    std::vector<std::vector<double>> ret = {{edges[0][0],edges[0][1]}}; // Start at first point
    std::cout << "Edges Size" << edges.size() << "\n";
    for(int i = 1; i < edges.size()+1; i++) {
        std::vector<double> Pcurr = ret[i-1]; // Current Point
        // Loop over edges and find next Point
        for(int j = 0; j < edges.size(); j++) {
            std::vector<double> Pfirst = edges[j];
            if(Pfirst[0] == Pcurr[0] && Pfirst[1] == Pcurr[1]) {
                // std::cout << "j = " << j << "\n";
                ret.push_back({Pfirst[2],Pfirst[3]});
                break;
            }
        }
    }
    myfile.open("marching_cubes_output/edgeOutput-ret.txt");
    for(int i = 0; i < ret.size(); i++) {
        for(int j = 0; j < ret[0].size(); j++) {
            myfile << ret[i][j] << ",";
            std::cout << ret[i][j] << ",";
        }
        myfile << "\n";
        std::cout << "\n";
    }
    myfile.close();

    // At this point we have the Ordered set of points, now let us get the tangents 
    std::vector<std::vector<double>> tangents = {};
    std::cout << "ret Size = " << ret.size() << "\n";
    for(int i = 0; i < ret.size(); i++) {
        int iNext = (i+1);
        int iPrev = (i-1);
        while(iNext >= ret.size()){
            iNext -= ret.size();
        }
        while(iPrev < 0) {
            iPrev += ret.size();
        }
        // Tangent points from previous point to next point 
        std::vector<double> Pprev = ret[iPrev];
        std::vector<double> Pcurr = ret[i];
        std::vector<double> Pnext = ret[iNext];

        std::vector<double> tCurr = {Pnext[0]-Pprev[0],Pnext[1]-Pprev[1]};
        //Normalize
        double tCurrMag = std::sqrt(tCurr[0]*tCurr[0]+tCurr[1]*tCurr[1]);
        tCurr[0] = tCurr[0]/tCurrMag;
        tCurr[1] = tCurr[1]/tCurrMag;
        tangents.push_back(tCurr);
    }
    std::cout << "Tangent Size = " << tangents.size(); 
    // Print Ret and Tangents
    std::cout << std::setprecision(3);
    for(int i =0; i < ret.size(); i++) {
        std::cout <<"(" <<ret[i][0] << ","<<ret[i][1] << ") at (" <<tangents[i][0] << ","<<tangents[i][1] << ")\n";
    }
    // Make Spline
    Spline s = Spline<double>::LocalRQuadInterp(ret,tangents);
    s.saveToVTK("marching_cubes_output/SplineFromMarchingCubes");
    s.printControlPoints();
    return 0;
}