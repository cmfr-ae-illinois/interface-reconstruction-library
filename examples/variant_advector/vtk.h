// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_SIMPLE_VARIANT_ADVECTOR_VTK_H_
#define EXAMPLES_SIMPLE_VARIANT_ADVECTOR_VTK_H_

#include <string>
#include <vector>

#include "irl/cylinder_reconstruction/cylinder_parametrized_surface.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/paraboloid_reconstruction/paraboloid_parametrized_surface.h"
#include "irl/surface_mesher/triangulated_surface.h"

#include "examples/variant_advector/basic_mesh.h"
#include "examples/variant_advector/data.h"

struct InterfaceScalarField {
  std::string name;
  Data<double> polygon_scalar_data;
  Data<double> paraboloid_scalar_data;
  std::vector<double> flattened_polygon_scalar_data;
  std::vector<double> flattened_paraboloid_scalar_data;

  InterfaceScalarField() = default;

  InterfaceScalarField(const std::string& n, const BasicMesh* mesh)
      : name(n), polygon_scalar_data(mesh), paraboloid_scalar_data(mesh) {}
};

class VTKOutput {
  struct DataIO {
    DataIO(const std::string& a_name, const Data<double>& a_data)
        : name(a_name), pointer(&a_data) {}

    std::string name;
    const Data<double>* pointer;
  };

 public:
  VTKOutput(const std::string& a_directory, const std::string& a_file_name_base,
            const BasicMesh& a_mesh);

  void addData(const std::string& a_name, const Data<double>& a_data);

  void writeVTKFile(const double a_time);

  void writeVTKInterface(
      const double a_time, std::vector<IRL::Polygon>& a_polygons,
      std::vector<IRL::ParaboloidParametrizedSurfaceOutput>& a_paraboloids,
      std::vector<IRL::CylinderParametrizedSurfaceOutput>& a_cylinders,
      const bool a_print_info = false);

  void writeVTKInterfaceWithScalar(
      const double a_time, std::vector<IRL::Polygon>& a_polygons,
      std::vector<IRL::ParaboloidParametrizedSurfaceOutput>& a_paraboloids,
      std::vector<IRL::CylinderParametrizedSurfaceOutput>& a_cylinders,
      std::vector<InterfaceScalarField>* a_scalar_fields,
      const bool a_print_info = false);

  void writeParametrizedInterface(
      const double a_time,
      std::vector<IRL::ParaboloidParametrizedSurfaceOutput>& a_surface);

 private:
  std::string directory_m;
  std::string file_name_base_m;
  std::size_t data_files_written_m;
  std::size_t interface_files_written_m;
  const BasicMesh* mesh_m;
  std::vector<DataIO> data_to_write_m;
};

void writeCellsVTK(const std::string& path,
                   const std::vector<std::pair<IRL::Pt, IRL::Pt>>& cells);

bool writeScatterVTK(const std::vector<IRL::Pt>& points,
                     const std::string& filename);

void writePolygonVTK(const std::string& filename, const IRL::Polygon& poly);

void writePolygonsVTK(const std::string& filename,
                      const std::vector<IRL::Polygon>& polygons);

void writeVectorsVTK(const std::string& file, const std::vector<IRL::Pt>& pts,
                     const std::vector<IRL::Normal>& vecs);

bool writeLinesVTK(const std::vector<std::pair<IRL::Pt, IRL::Pt>>& lines,
                   const std::string& filename);

bool writePlanePatchVTK(const Eigen::Vector3d& origin,
                        const Eigen::Vector3d& normal, double h,
                        const std::string& filename);

void writeParaboloidVTK(const std::string& filename,
                        const IRL::Paraboloid& paraboloid, const int& Nx,
                        const int& Ny, const double& xmin, const double& xmax,
                        const double& ymin, const double& ymax);

#endif  // EXAMPLES_SIMPLE_VARIANT_ADVECTOR_VTK_H_
