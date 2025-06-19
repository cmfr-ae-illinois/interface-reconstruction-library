
#include "examples/variant_advector/vtk.h"

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
  fprintf(file, "<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n",
          mesh_m->imin(), mesh_m->imax() + 1, mesh_m->jmin(),
          mesh_m->jmax() + 1, mesh_m->kmin(), mesh_m->kmax() + 1);
  fprintf(file, "<Piece Extent=\"%d %d %d %d %d %d\">\n", mesh_m->imin(),
          mesh_m->imax() + 1, mesh_m->jmin(), mesh_m->jmax() + 1,
          mesh_m->kmin(), mesh_m->kmax() + 1);

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
  for (int i = mesh_m->kmin(); i <= mesh_m->kmax() + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(mesh_m->z(i)));
  }
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
    for (int k = mesh_m->kmin(); k <= mesh_m->kmax(); ++k) {
      for (int j = mesh_m->jmin(); j <= mesh_m->jmax(); ++j) {
        for (int i = mesh_m->imin(); i <= mesh_m->imax(); ++i) {
          fprintf(file, "%15.8E ",
                  static_cast<double>((*data.pointer)(i, j, k)));
        }
      }
    }
    fprintf(file, "\n</DataArray>\n");
  }
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n</RectilinearGrid>\n</VTKFile>\n");
  fclose(file);
  ++data_files_written_m;
}

void VTKOutput::writeParametrizedInterface(
    const double a_time,
    std::vector<IRL::ParaboloidParametrizedSurfaceOutput>& a_surface) {
  const auto surface_file_name =
      directory_m + "/" + file_name_base_m + "_interface_" +
      std::to_string(interface_files_written_m) + ".irl";
  FILE* file;

  int dummy1, dummy2;
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status status;

  int number_of_surfaces = a_surface.size(), total_surfaces = 0;
  MPI_Allreduce(&number_of_surfaces, &total_surfaces, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

  if (rank == 0) {
    file = fopen(surface_file_name.c_str(), "w");
    fprintf(file, "Number of surface patches: %i\n", total_surfaces);
    fclose(file);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  for (int r = 0; r < size; r++) {
    if (rank == r) {
      if (rank > 0) {
        MPI_Recv(&dummy1, 1, MPI_INT, r - 1, 1234 + r - 1, MPI_COMM_WORLD,
                 &status);
      }
      file = fopen(surface_file_name.c_str(), "a");
      for (std::size_t i = 0; i < a_surface.size(); ++i) {
        auto paraboloid = a_surface[i].getParaboloid();
        auto datum = paraboloid.getDatum();
        auto ref_frame = paraboloid.getReferenceFrame();
        auto aligned_paraboloid = paraboloid.getAlignedParaboloid();
        auto arc_list = a_surface[i].getArcs();
        fprintf(file, "Number of arcs: %i\n",
                static_cast<int>(arc_list.size()));
        fprintf(file,
                "Reference frame: ( %+.16e %+.16e %+.16e ) ( %+.16e %+.16e "
                "%+.16e ) ( %+.16e %+.16e %+.16e )\n",
                static_cast<double>(ref_frame[0][0]),
                static_cast<double>(ref_frame[0][1]),
                static_cast<double>(ref_frame[0][2]),
                static_cast<double>(ref_frame[1][0]),
                static_cast<double>(ref_frame[1][1]),
                static_cast<double>(ref_frame[1][2]),
                static_cast<double>(ref_frame[2][0]),
                static_cast<double>(ref_frame[2][1]),
                static_cast<double>(ref_frame[2][2]));
        fprintf(file, "Datum: ( %+.16e %+.16e %+.16e )\n",
                static_cast<double>(datum[0]), static_cast<double>(datum[1]),
                static_cast<double>(datum[2]));
        fprintf(file, "Coefficients: ( %+.16e %+.16e )\n",
                static_cast<double>(aligned_paraboloid.a()),
                static_cast<double>(aligned_paraboloid.b()));
        for (std::size_t j = 0; j < arc_list.size(); ++j) {
          const auto arc = arc_list[j];
          fprintf(file,
                  "Arc %i: %i %i %+.16e ( %+.16e %+.16e %+.16e ) ( %+.16e "
                  "%+.16e %+.16e ) ( %+.16e %+.16e %+.16e )\n",
                  j, arc.start_point_id(), arc.end_point_id(), arc.weight(),
                  arc.start_point()[0], arc.start_point()[1],
                  arc.start_point()[2], arc.control_point()[0],
                  arc.control_point()[1], arc.control_point()[2],
                  arc.end_point()[0], arc.end_point()[1], arc.end_point()[2]);
        }
      }

      fclose(file);
      if (size > 1 && rank < size - 1) {
        MPI_Send(&dummy2, 1, MPI_INT, r + 1, 1234 + r, MPI_COMM_WORLD);
      }
    }
  }
  ++interface_files_written_m;
}

void VTKOutput::writeVTKInterface(
    const double a_time, std::vector<IRL::Polygon>& a_polygons,
    std::vector<IRL::ParaboloidParametrizedSurfaceOutput>& a_paraboloids,
    std::vector<IRL::CylinderParametrizedSurfaceOutput>& a_cylinders,
    const bool a_print_info) {
  const auto surface_file_name = directory_m + "/" + file_name_base_m +
                                 "_interface_" +
                                 std::to_string(interface_files_written_m);
  IRL::MixedPolygonBezierSurface surface_output;
  for (int i = 0; i < a_paraboloids.size(); ++i) {
    surface_output.addSurface(
        a_paraboloids[i].getQuadraticBezierTriangleApprox());
  }
  for (int i = 0; i < a_cylinders.size(); ++i) {
    // surface_output.addSurface(
    //   a_cylinders[i].getQuadraticBezierTriangleApprox());
  }
  surface_output.addPolygons(a_polygons);
  surface_output.write(surface_file_name);

  ++interface_files_written_m;
}
