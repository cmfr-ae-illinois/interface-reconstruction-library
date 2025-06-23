#include <array>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>

#include "examples/spline_optimization/basic_mesh.h"
#include "examples/spline_optimization/irl2d.h"
#include "examples/spline_optimization/data.h"
#include "irl/splines/Spline.h"
#include "irl/splines/MarchingSquares.h"
int main() {
    // Field Initialization
    const int a_nx = 50;
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
    std::cout << "imax = " << mesh.imax() << "\n";
    for(int i =0; i <= a_nx; i++) {
        std::cout << mesh.x(i) << ",";
    }
    std::cout << "\n";
    // Initialize
    const auto circle_center = IRL2D::Vec(0,0);
    const double circle_radius = 0.15;

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

    // Treshold
    double tresh = 0.5;
    
    MarchingSquares<double> testMS(&mesh,tresh,2);
    std::vector<std::vector<double>> ret = testMS.vertexPoints(VOFarray); // Note this does not have the repeating start point
    
    std::vector<std::vector<double>> tangents = {};
    // Here we calculate the tangents. This is done by assuming it will be periodic and estimating the tangent as the vector
    // Pointing from the previous point to the next point (central difference like)
    for(int i = 0; i < ret.size(); i++) {
        int iNext = (i+1);
        int iPrev = (i-1);
        while(iNext >= ret.size()) {
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

    std::vector<std::vector<double>> ret2 = {ret[0]};
    std::vector<std::vector<double>> tangents2 = {tangents[0]};
    double dx = mesh.dx();
    // std::cout <<"Here 120\n";
    for(int i = 1; i < ret.size(); i++) {
        std::vector<double> pcurr = ret[i];
        bool add = true;
        for(int j = 0; j < ret2.size(); j++) {
            std::vector<double> ptest = ret2[j];
            if(fabs(pcurr[0] - ptest[0]) < dx/2 && fabs(pcurr[1] - ptest[1]) < dx/2) {
                add = false;
                break;
            }
        }
        if(add) {
            ret2.push_back(pcurr);
        }
    }

    // Recalculate Tangents

    for(int i = 0; i < ret2.size(); i++) {
        int iNext = (i+1);
        int iPrev = (i-1);
        while(iNext >= ret2.size()) {
            iNext -= ret2.size();
        }
        while(iPrev < 0) {
            iPrev += ret2.size();
        }
        // Tangent points from previous point to next point 
        std::vector<double> Pprev = ret2[iPrev];
        std::vector<double> Pcurr = ret2[i];
        std::vector<double> Pnext = ret2[iNext];

        std::vector<double> tCurr = {Pnext[0]-Pprev[0],Pnext[1]-Pprev[1]};
        //Normalize
        double tCurrMag = std::sqrt(tCurr[0]*tCurr[0]+tCurr[1]*tCurr[1]);
        tCurr[0] = tCurr[0]/tCurrMag;
        tCurr[1] = tCurr[1]/tCurrMag;
        tangents2.push_back(tCurr);
    }
    // Now we have points and tangents. Make them periodic
    ret.push_back(ret[0]);
    tangents.push_back(tangents[0]);
    ret2.push_back(ret2[0]);
    tangents2.push_back(tangents2[0]);
    Spline s = Spline<double>::LocalRQuadInterp(ret2,tangents2);
    // s.meshSurfaceTensionForce(&mesh,1);
    s.saveToVTK("marching_cubes_output/FilteredSplineFromMarchingCubes_N="+std::to_string(a_nx));
    return 0;
}