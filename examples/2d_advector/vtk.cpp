
#include "examples/2d_advector/vtk.h"

#include <mpi.h>
#include <stdio.h>
#include <sys/stat.h>

VTKOutput::VTKOutput(const std::string& a_directory,
                     const std::string& a_file_name_base,
                     const BasicMesh& a_mesh)
    : directory_m(a_directory),
      file_name_base_m(a_file_name_base),
      data_files_written_m(0),
      interface_files_written_m(0),
      mesh_m(&a_mesh),
      data_to_write_m() {
  const int dir_err = mkdir(directory_m.c_str(), 0777);
}

void VTKOutput::addData(const std::string& a_name, const Data<double>& a_data) {
  data_to_write_m.push_back(DataIO(a_name, a_data));
}

void VTKOutput::writeVTKFile(const double a_time) {
  const auto file_name = directory_m + "/" + file_name_base_m + "_" +
                         std::to_string(data_files_written_m) + ".vtr";

  FILE* file;
  file = fopen(file_name.c_str(), "w");

  fprintf(file, "<VTKFile type=\"RectilinearGrid\">\n");
  fprintf(file, "<RectilinearGrid WholeExtent=\"%d %d %d %d 0 1\">\n",
          mesh_m->imin(), mesh_m->imax() + 1, mesh_m->jmin(),
          mesh_m->jmax() + 1);
  fprintf(file, "<Piece Extent=\"%d %d %d %d 0 1\">\n", mesh_m->imin(),
          mesh_m->imax() + 1, mesh_m->jmin(), mesh_m->jmax() + 1);

  fprintf(file, "<Coordinates>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = mesh_m->imin(); i <= mesh_m->imax() + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(mesh_m->x(i)));
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = mesh_m->jmin(); i <= mesh_m->jmax() + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(mesh_m->y(i)));
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  fprintf(file, "%15.8E ", static_cast<double>(0.0));
  fprintf(file, "%15.8E ", static_cast<double>(0.0));
  fprintf(file, "\n</DataArray>\n");

  fprintf(file, "</Coordinates>\n");

  fprintf(file, "<PointData>\n</PointData>\n");

  fprintf(file, "<CellData Scalars=\"");
  for (auto& data : data_to_write_m) {
    fprintf(file, "%s ", data.name.c_str());
  }
  fprintf(file, "\" >\n");
  for (auto& data : data_to_write_m) {
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\" "
            "format=\"ascii\">\n",
            data.name.c_str());
    for (int j = mesh_m->jmin(); j <= mesh_m->jmax(); ++j) {
      for (int i = mesh_m->imin(); i <= mesh_m->imax(); ++i) {
        fprintf(file, "%15.8E ", static_cast<double>((*data.pointer)(i, j)));
      }
    }
    fprintf(file, "\n</DataArray>\n");
  }
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n</RectilinearGrid>\n</VTKFile>\n");
  fclose(file);
  ++data_files_written_m;
}

void VTKOutput::writeVTKInterface(const double a_time,
                                  const Data<IRL2D::Parabola> a_interface) {
  const double tol = 5. * std::numeric_limits<double>::epsilon();
  const auto surface_file_name = directory_m + "/" + file_name_base_m +
                                 "_interface_" +
                                 std::to_string(interface_files_written_m);

  std::vector<IRL2D::BezierList> interfaces;
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const auto cell =
          IRL2D::RectangleFromBounds(IRL2D::Vec(mesh.x(i), mesh.y(j)),
                                     IRL2D::Vec(mesh.x(i + 1), mesh.y(j + 1)));
      const auto vfrac = IRL2D::ComputeVFrac(cell, a_interface(i, j));
      if (vfrac > 1.0e-8 && vfrac < 1.0 - 1.0e-8) {
        const auto arcs = IRL2D::ParabolaClip(cell, a_interface(i, j), true);
        for (int k = 0; k < arcs.size(); k += 2) {
          interfaces.push_back(IRL2D::BezierList{arcs[k], arcs[k + 1]});
        }
      }
    }
  }
  IRL2D::ToVTK(interfaces, surface_file_name);

  ++interface_files_written_m;
}