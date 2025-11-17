
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

void writeCellsVTK(const std::string& path,
                   const std::vector<std::pair<IRL::Pt, IRL::Pt>>& cells) {
  const int pts_per_cell = 8;
  const int vtk_hex_type = 12;  // VTK_HEXAHEDRON

  const std::size_t nCells = cells.size();
  const std::size_t nPoints = nCells * pts_per_cell;

  std::ofstream out(path);
  if (!out) return;

  out << "# vtk DataFile Version 4.2\n";
  out << "Axis-aligned hex cells\n";
  out << "ASCII\n";
  out << "DATASET UNSTRUCTURED_GRID\n";
  out << "POINTS " << nPoints << " float\n";

  // 0:(x0,y0,z0) 1:(x1,y0,z0) 2:(x1,y1,z0) 3:(x0,y1,z0)
  // 4:(x0,y0,z1) 5:(x1,y0,z1) 6:(x1,y1,z1) 7:(x0,y1,z1)
  std::vector<std::array<float, 3>> points;
  points.reserve(nPoints);
  for (const auto& b : cells) {
    float x0 = b.first[0], x1 = b.second[0], y0 = b.first[1], y1 = b.second[1],
          z0 = b.first[2], z1 = b.second[2];
    points.push_back({x0, y0, z0});
    points.push_back({x1, y0, z0});
    points.push_back({x1, y1, z0});
    points.push_back({x0, y1, z0});
    points.push_back({x0, y0, z1});
    points.push_back({x1, y0, z1});
    points.push_back({x1, y1, z1});
    points.push_back({x0, y1, z1});
  }
  for (auto& p : points) out << p[0] << " " << p[1] << " " << p[2] << "\n";

  out << "CELLS " << nCells << " " << nCells * (1 + pts_per_cell) << "\n";
  for (std::size_t c = 0; c < nCells; ++c) {
    std::size_t base = c * pts_per_cell;
    out << pts_per_cell;
    for (int k = 0; k < pts_per_cell; ++k) out << " " << (base + k);
    out << "\n";
  }

  out << "CELL_TYPES " << nCells << "\n";
  for (std::size_t c = 0; c < nCells; ++c) out << vtk_hex_type << "\n";
  out.close();
}

bool writeScatterVTK(const std::vector<IRL::Pt>& points,
                     const std::string& filename) {
  if (points.empty()) {
    std::cerr << "No points to write.\n";
    return false;
  }

  std::ofstream os(filename);
  if (!os) {
    std::cerr << "Could not open file: " << filename << "\n";
    return false;
  }

  os << "# vtk DataFile Version 3.0\n";
  os << "3D scatter points\n";
  os << "ASCII\n";
  os << "DATASET POLYDATA\n";

  os << "POINTS " << points.size() << " float\n";
  for (auto& p : points) {
    os << static_cast<double>(p[0]) << " " << static_cast<double>(p[1]) << " "
       << static_cast<double>(p[2]) << "\n";
  }

  const std::size_t n = points.size();
  os << "VERTICES " << n << " " << 2 * n << "\n";
  for (std::size_t i = 0; i < n; ++i) {
    os << 1 << " " << i << "\n";
  }

  return true;
}

void writePolygonVTK(const std::string& filename,
                     const std::vector<IRL::Pt>& poly) {
  std::ofstream out(filename);
  out << "# vtk DataFile Version 4.2\n";
  out << "Plane patch inside cell\n";
  out << "ASCII\n";
  out << "DATASET POLYDATA\n";
  out << "POINTS " << poly.size() << " float\n";
  for (auto& p : poly) {
    out << p.x() << " " << p.y() << " " << p.z() << "\n";
  }

  out << "POLYGONS 1 " << poly.size() + 1 << "\n";
  out << poly.size();
  for (size_t i = 0; i < poly.size(); ++i) {
    out << " " << i;
  }
  out << "\n";
}

void writePolygonsVTK(const std::string& filename,
                      const std::vector<std::vector<IRL::Pt>>& polygons) {
  std::vector<std::vector<IRL::Pt>> polys;
  polys.reserve(polygons.size());
  for (const auto& poly : polygons)
    if (poly.size() >= 3) polys.push_back(poly);

  // Count totals
  std::size_t totalPoints = 0;
  std::size_t totalPolyRecs = 0;  // sum of (1 + n_i)
  for (const auto& poly : polys) {
    totalPoints += poly.size();
    totalPolyRecs += 1 + poly.size();
  }

  std::ofstream out(filename);
  if (!out) return;

  out << "# vtk DataFile Version 4.2\n";
  out << "Multiple plane patches\n";
  out << "ASCII\n";
  out << "DATASET POLYDATA\n";

  // POINTS block
  out << "POINTS " << totalPoints << " float\n";
  for (const auto& poly : polys) {
    for (const auto& p : poly) {
      out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
  }

  // POLYGONS block
  out << "POLYGONS " << polys.size() << " " << totalPolyRecs << "\n";
  std::size_t offset = 0;
  for (const auto& poly : polys) {
    out << poly.size();
    for (std::size_t k = 0; k < poly.size(); ++k) out << " " << (offset + k);
    out << "\n";
    offset += poly.size();
  }
}

void writeVectorsVTK(const std::string& file, const std::vector<IRL::Pt>& pts,
                     const std::vector<IRL::Normal>& vecs) {
  std::ofstream out(file);
  out << "# vtk DataFile Version 4.2\nVectors\nASCII\n";
  out << "DATASET POLYDATA\n";
  out << "POINTS " << pts.size() << " float\n";
  for (auto& p : pts) out << p.x() << " " << p.y() << " " << p.z() << "\n";

  out << "POINT_DATA " << vecs.size() << "\n";
  out << "VECTORS myVectors float\n";
  for (auto& v : vecs) out << v[0] << " " << v[1] << " " << v[2] << "\n";
}

bool writeLinesVTK(const std::vector<std::pair<IRL::Pt, IRL::Pt>>& lines,
                   const std::string& filename) {
  std::ofstream os(filename);
  if (!os) {
    std::cerr << "Could not open file: " << filename << "\n";
    return false;
  }

  // Header
  os << "# vtk DataFile Version 3.0\n";
  os << "Line segments\n";
  os << "ASCII\n";
  os << "DATASET POLYDATA\n";

  // Collect points
  std::size_t npoints = 2 * lines.size();
  os << "POINTS " << npoints << " float\n";
  for (auto& seg : lines) {
    os << seg.first[0] << " " << seg.first[1] << " " << seg.first[2] << "\n";
    os << seg.second[0] << " " << seg.second[1] << " " << seg.second[2] << "\n";
  }

  // Define lines (cells)
  std::size_t ncells = lines.size();
  os << "LINES " << ncells << " " << ncells * 3 << "\n";
  for (std::size_t i = 0; i < lines.size(); ++i) {
    os << "2 " << 2 * i << " " << 2 * i + 1 << "\n";
  }

  return true;
}

bool writePlanePatchVTK(const Eigen::Vector3d& origin,
                        const Eigen::Vector3d& normal, double h,
                        const std::string& filename) {
  Eigen::Vector3d n = normal.normalized();
  Eigen::Vector3d a = (std::fabs(n.x()) < 0.9) ? Eigen::Vector3d::UnitX()
                                               : Eigen::Vector3d::UnitY();
  // Orthonormal basis (u,v) spanning the plane
  Eigen::Vector3d u = n.cross(a).normalized();
  Eigen::Vector3d v = n.cross(u);

  // Square corners (counter-clockwise)
  Eigen::Vector3d p0 = origin - 1.25 * h * u - 0.75 * h * v;
  Eigen::Vector3d p1 = origin + 1.25 * h * u - 0.75 * h * v;
  Eigen::Vector3d p2 = origin + 1.25 * h * u + 0.5 * h * v;
  Eigen::Vector3d p3 = origin - 1.25 * h * u + 0.5 * h * v;

  std::ofstream os(filename);
  if (!os) return false;

  os << "# vtk DataFile Version 3.0\n";
  os << "Plane patch\n";
  os << "ASCII\n";
  os << "DATASET POLYDATA\n";

  os << "POINTS 4 float\n";
  auto emit = [&](const Eigen::Vector3d& p) {
    os << static_cast<float>(p.x()) << " " << static_cast<float>(p.y()) << " "
       << static_cast<float>(p.z()) << "\n";
  };
  emit(p0);
  emit(p1);
  emit(p2);
  emit(p3);

  os << "POLYGONS 1 5\n";
  os << "4 0 1 2 3\n";

  os << "CELL_DATA 1\n";
  os << "NORMALS plane_normals float\n";
  os << static_cast<float>(n.x()) << " " << static_cast<float>(n.y()) << " "
     << static_cast<float>(n.z()) << "\n";

  return true;
}
