// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/2d_advector/solver.h"

#include <stdio.h>
#include <cstdio>
#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iomanip>

#include "examples/2d_advector/basic_mesh.h"
#include "examples/2d_advector/data.h"
#include "examples/2d_advector/rotation_2d.h"

void setPhaseQuantities(const Data<IRL2D::Parabola>& a_interface,
                        Data<IRL2D::Moments>* a_liquid_moments,
                        Data<IRL2D::Moments>* a_gas_moments,
                        Data<double>* a_vfrac) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      const auto cell =
          IRL2D::RectangleFromBounds(IRL2D::Vec(mesh.x(i), mesh.y(j)),
                                     IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
      const auto cell_moments = IRL2D::ComputeMoments(cell);
      (*a_liquid_moments)(i, j) =
          IRL2D::ComputeMoments(cell, a_interface(i, j));
      (*a_gas_moments)(i, j) = cell_moments - (*a_liquid_moments)(i, j);
      (*a_vfrac)(i, j) = (*a_liquid_moments)(i, j).m0() / mesh.cell_volume();
    }
  }
  a_liquid_moments->updateBorder();
  a_gas_moments->updateBorder();
  a_vfrac->updateBorder();
}

void writeDiagnosticsHeader(void) {
  printf("%10s %20s %12s %20s %20s %20s %20s %20s %20s %20s %16s\n", "Iteration",
         "Time", "CFL", "liquidVFSum", "liquidVolSum", "ChangeLiquidVFSum",
         "ChangeLiquidVolSum", "AdvectionDuration", "ReconDuration",
         "OutputDuration", "InterfaceCells");
}

void writeOutDiagnostics(const int a_iteration, const double a_dt,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<IRL2D::Moments>& a_liquid_moments,
                         const Data<IRL2D::Parabola>& a_interface,
                         std::chrono::duration<double> a_VOF_duration,
                         std::chrono::duration<double> a_recon_duration,
                         std::chrono::duration<double> a_write_duration) {
  const BasicMesh& mesh = a_U.getMesh();
  static double initial_liquid_volume_fraction_sum;
  static double initial_liquid_volume_sum;
  // Calculate CFL
  double CFL = -DBL_MAX;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      CFL = std::fmax(CFL, std::fmax(a_U(i, j) * a_dt / mesh.dx(),
                                     a_V(i, j) * a_dt / mesh.dy()));
    }
  }

  // Calculate sum of volume fraction and sum of liquid volume
  double liquid_volume_fraction_sum = 0.0;
  double liquid_volume_sum = 0.0;
  int number_of_interface_cells = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          a_liquid_moments(i, j).m0() / mesh.cell_volume();
      liquid_volume_fraction_sum += liquid_volume_fraction;
      liquid_volume_sum += a_liquid_moments(i, j).m0();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        ++number_of_interface_cells;
      }
    }
  }
  // Save initial values to compare against.
  if (a_iteration == 0) {
    initial_liquid_volume_fraction_sum = liquid_volume_fraction_sum;
    initial_liquid_volume_sum = liquid_volume_sum;
  }
  printf(
      "%10d %20.4E %12.3F %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %16d"
      "\n",
      a_iteration, a_simulation_time, CFL, liquid_volume_fraction_sum,
      liquid_volume_sum,
      liquid_volume_fraction_sum - initial_liquid_volume_fraction_sum,
      liquid_volume_sum - initial_liquid_volume_sum, a_VOF_duration.count(),
      a_recon_duration.count(), a_write_duration.count(),
      number_of_interface_cells);
}

void printError(const BasicMesh& mesh,
                const Data<IRL2D::Moments>& liquid_moments,
                const Data<IRL2D::Moments>& starting_liquid_moments) {
  double linf_error_m0 = 0.0;
  double linf_error_m1 = 0.0;
  double linf_error_m2 = 0.0;
  double l1_error_m0 = 0.0;
  double l1_error_m1 = 0.0;
  double l1_error_m2 = 0.0;
  double l2_error_m0 = 0.0;
  double l2_error_m1 = 0.0;
  double l2_error_m2 = 0.0;
  double scale_m0 = 1.0 / std::pow(mesh.dx(), 2.0);
  double scale_m1 = 1.0 / std::pow(mesh.dx(), 3.0);
  double scale_m2 = 1.0 / std::pow(mesh.dx(), 4.0);
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const double liquid_volume_fraction =
          liquid_moments(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        auto mom_err = (liquid_moments(i, j) - starting_liquid_moments(i, j));
        linf_error_m0 = std::max(linf_error_m0, std::abs(mom_err.m0()));
        linf_error_m1 = std::max(linf_error_m1, std::abs(mom_err.m1()[0]));
        linf_error_m1 = std::max(linf_error_m1, std::abs(mom_err.m1()[1]));
        linf_error_m2 = std::max(linf_error_m2, std::abs(mom_err.m2()[0][0]));
        linf_error_m2 = std::max(linf_error_m2, std::abs(mom_err.m2()[1][0]));
        linf_error_m2 = std::max(linf_error_m2, std::abs(mom_err.m2()[1][1]));
        linf_error_m2 = std::max(linf_error_m2, std::abs(mom_err.m2()[0][1]));
        l1_error_m0 += std::abs(mom_err.m0());
        l1_error_m1 += std::abs(mom_err.m1()[0]);
        l1_error_m1 += std::abs(mom_err.m1()[1]);
        l1_error_m2 += std::abs(mom_err.m2()[0][0]);
        l1_error_m2 += std::abs(mom_err.m2()[1][0]);
        l1_error_m2 += std::abs(mom_err.m2()[1][1]);
        l1_error_m2 += std::abs(mom_err.m2()[0][1]);
        l2_error_m0 += mom_err.m0() * mom_err.m0();
        l2_error_m1 += mom_err.m1()[0] * mom_err.m1()[0];
        l2_error_m1 += mom_err.m1()[1] * mom_err.m1()[1];
        l2_error_m2 += mom_err.m2()[0][0] * mom_err.m2()[0][0];
        l2_error_m2 += mom_err.m2()[1][0] * mom_err.m2()[1][0];
        l2_error_m2 += mom_err.m2()[1][1] * mom_err.m2()[1][1];
        l2_error_m2 += mom_err.m2()[0][1] * mom_err.m2()[0][1];
      }
    }
  }
  l1_error_m0 /= (static_cast<double>(mesh.getNx() * mesh.getNy()));
  l1_error_m1 /= (static_cast<double>(mesh.getNx() * mesh.getNy()));
  l1_error_m2 /= (static_cast<double>(mesh.getNx() * mesh.getNy()));
  l2_error_m0 /= (static_cast<double>(mesh.getNx() * mesh.getNy()));
  l2_error_m1 /= (static_cast<double>(mesh.getNx() * mesh.getNy()));
  l2_error_m2 /= (static_cast<double>(mesh.getNx() * mesh.getNy()));
  linf_error_m0 *= scale_m0;
  linf_error_m1 *= scale_m1;
  linf_error_m2 *= scale_m2;
  l1_error_m0 *= scale_m0;
  l1_error_m1 *= scale_m1;
  l1_error_m2 *= scale_m2;
  l2_error_m0 = std::sqrt(l2_error_m0) * scale_m0;
  l2_error_m1 = std::sqrt(l2_error_m1) * scale_m1;
  l2_error_m2 = std::sqrt(l2_error_m2) * scale_m2;
  std::cout << std::scientific << std::setprecision(3)
            << "Linf M0 = " << linf_error_m0 << std::endl;
  std::cout << "Linf M1 = " << linf_error_m1 << std::endl;
  std::cout << "Linf M2 = " << linf_error_m2 << std::endl;
  std::cout << "L1   M0 = " << l1_error_m0 << std::endl;
  std::cout << "L1   M1 = " << l1_error_m1 << std::endl;
  std::cout << "L1   M2 = " << l1_error_m2 << std::endl;
  std::cout << "L2   M0 = " << l2_error_m0 << std::endl;
  std::cout << "L2   M1 = " << l2_error_m1 << std::endl;
  std::cout << "L2   M2 = " << l2_error_m2 << std::endl;
}

// output data to csv
void writeToCSV(const BasicMesh& mesh, 
                const Data<IRL2D::Parabola>& a_interface,
                const Data<IRL2D::Moments>& a_liquid_moments,
                const std::string& filepath){

  // finding mixed cells
  Data<int> band(&mesh);
  for (int i = mesh.imin(); i < mesh.imax(); ++i){
    for (int j = mesh.jmin(); j < mesh.jmax(); ++j){
      const double liquid_volume_fraction =
          (a_liquid_moments)(i, j).m0() / mesh.cell_volume();
      if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
          liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
        band(i, j) = 1;
      }
    }
  }

  // csv file
  // stores indices, end points, normal and tangent components
  std::ofstream csvfile(filepath);
  csvfile << "i,j,ax,ay,bx,by,tx,ty,nx,ny,vf\n";
  double ax, ay, bx, by, tx, ty, nx, ny, vf;
  // storing interface endpoints for all line segments
  for (int i = mesh.imin(); i < mesh.imax(); ++i){
    for (int j = mesh.jmin(); j < mesh.jmax(); ++j){
      if (band(i,j) == 1){
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
          IRL2D::Vec(mesh.x(i), mesh.y(j)), 
          IRL2D::Vec(mesh.x(i+1), mesh.y(j+1))
        );
        IRL2D::BezierList clipped_plic = IRL2D::ParabolaClip(cell, a_interface(i,j), true);
        ax = clipped_plic[0].first[0]; ay = clipped_plic[0].first[1];
        bx = clipped_plic[1].first[0]; by = clipped_plic[1].first[1];
        tx = a_interface(i,j).frame()[0][0]; ty = a_interface(i,j).frame()[0][1];
        nx = a_interface(i,j).frame()[1][0]; ny = a_interface(i,j).frame()[1][1];
        vf = a_liquid_moments(i,j).m0() / IRL2D::ComputeArea(cell);
        csvfile << i << "," << j << ',' << ax << "," << ay << ","
                << bx << "," << by << "," << tx << "," << ty << ","
                << nx << "," << ny << "," << vf << "\n";
      }
    }
  }
}

void readCSV(const std::string& filepath,
             std::vector<IRL2D::InterfaceEndPoints>& data){
  
  std::ifstream file(filepath);

  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filepath << std::endl;
  }

  std::string line;

  // skip header
  std::getline(file, line);

  while (std::getline(file, line)){
    std::stringstream ss(line);   // separates by ,
    std::string info;
    IRL2D::InterfaceEndPoints iep;

    // reading in each value and storing it in struct
    std::getline(ss, info, ','); iep.xIndex = std::stoi(info);
    std::getline(ss, info, ','); iep.yIndex = std::stoi(info);
    std::getline(ss, info, ','); iep.ax = std::stod(info);
    std::getline(ss, info, ','); iep.ay = std::stod(info);
    std::getline(ss, info, ','); iep.bx = std::stod(info);
    std::getline(ss, info, ','); iep.by = std::stod(info);
    std::getline(ss, info, ','); iep.tx = std::stod(info);
    std::getline(ss, info, ','); iep.ty = std::stod(info);
    std::getline(ss, info, ','); iep.nx = std::stod(info);
    std::getline(ss, info, ','); iep.ny = std::stod(info);
    std::getline(ss, info, ','); iep.vf = std::stod(info);
    iep.mixed = true;
    data.push_back(iep);
  }
  file.close();
}