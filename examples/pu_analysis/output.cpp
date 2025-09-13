#include <iostream>

#include "examples/pu_analysis/functions.h"

#include <fstream>
#include <iomanip>
#include <string>

// volume fraction data
void write_vtr(const std::string& filepath, const Data<double>& vf,
               const BasicMesh& mesh) {
  std::ofstream out(filepath);
  if (!out.is_open()) {
    std::cerr << "Failed to open " << filepath << " for writing.\n";
    return;
  }
  out << std::fixed << std::setprecision(6);

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" "
         "byte_order=\"LittleEndian\">\n";
  out << "  <RectilinearGrid WholeExtent=\"0 " << mesh.imax() + 1 << " 0 "
      << mesh.jmax() + 1 << " 0 0\">\n";
  out << "    <Piece Extent=\"0 " << mesh.imax() + 1 << " 0 " << mesh.jmax() + 1
      << " 0 0\">\n";

  out << "      <CellData Scalars=\"volumeFraction\">\n";
  out << "        <DataArray type=\"Float32\" Name=\"volumeFraction\" "
         "format=\"ascii\">\n";
  for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
    for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
      out << "          " << vf(i, j) << "\n";
    }
  }
  out << "        </DataArray>\n";
  out << "      </CellData>\n";

  out << "      <Coordinates>\n";
  out << "        <DataArray type=\"Float32\" Name=\"X\" format=\"ascii\">\n";
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    out << "          " << mesh.x(i) << "\n";
  }
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float32\" Name=\"Y\" format=\"ascii\">\n";
  for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
    out << "          " << mesh.y(j) << "\n";
  }
  out << "        </DataArray>\n";

  out << "        <DataArray type=\"Float32\" Name=\"Z\" format=\"ascii\">\n";
  out << "          0.0\n";
  out << "        </DataArray>\n";
  out << "      </Coordinates>\n";

  out << "    </Piece>\n";
  out << "  </RectilinearGrid>\n";
  out << "</VTKFile>\n";
  out.close();
}

// nested grid
void write_vtu(const std::string& filepath,
               const std::vector<const Cell*>& cells) {
  std::ofstream out(filepath);
  if (!out.is_open()) {
    std::cerr << "Failed to open " << filepath << " for writing.\n";
    return;
  }

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n";
  out << "  <UnstructuredGrid>\n";
  out << "    <Piece NumberOfPoints=\"" << 4 * cells.size()
      << "\" NumberOfCells=\"" << cells.size() << "\">\n";

  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         "format=\"ascii\">\n";
  for (const Cell* c : cells) {
    double x = c->x, y = c->y, dx = c->dx;
    out << x << " " << y << " 0\n";
    out << x + dx << " " << y << " 0\n";
    out << x + dx << " " << y + dx << " 0\n";
    out << x << " " << y + dx << " 0\n";
  }
  out << "        </DataArray>\n";
  out << "      </Points>\n";

  out << "      <Cells>\n";

  // Connectivity
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
         "format=\"ascii\">\n";
  for (int i = 0; i < cells.size(); ++i) {
    int base = 4 * i;
    out << base << " " << base + 1 << " " << base + 2 << " " << base + 3
        << "\n";
  }
  out << "        </DataArray>\n";

  // Offsets
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" "
         "format=\"ascii\">\n";
  for (int i = 1; i <= cells.size(); ++i) out << 4 * i << "\n";
  out << "        </DataArray>\n";

  // Types (VTK_QUAD = 9)
  out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (int i = 0; i < cells.size(); ++i) out << "9\n";
  out << "        </DataArray>\n";

  out << "      </Cells>\n";

  out << "    </Piece>\n";
  out << "  </UnstructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}

// scatter point vtk file
bool writeScatterVTK(const std::vector<IRL2D::Vec>& points,
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
  os << "2D scatter points\n";
  os << "ASCII\n";
  os << "DATASET POLYDATA\n";

  os << "POINTS " << points.size() << " float\n";
  for (auto& p : points) {
    os << static_cast<float>(p[0]) << " " << static_cast<float>(p[1]) << " "
       << 0.0f << "\n";
  }

  const std::size_t n = points.size();
  os << "VERTICES " << n << " " << 2 * n << "\n";
  for (std::size_t i = 0; i < n; ++i) {
    os << 1 << " " << i << "\n";
  }

  std::cout << "Wrote " << filename << " with " << n << " points.\n";
  return true;
}

// for pu implicit surface
void writeVTKStructuredPoints_PointData(const std::string& filename, int nx_pts,
                                        int ny_pts, double xmin, double xmax,
                                        double ymin, double ymax,
                                        const std::vector<double>& Fpoints) {
  assert((int)Fpoints.size() == nx_pts * ny_pts);

  const double dx = (nx_pts > 1) ? (xmax - xmin) / (nx_pts - 1) : 1.0;
  const double dy = (ny_pts > 1) ? (ymax - ymin) / (ny_pts - 1) : 1.0;

  std::ofstream out(filename);
  out << "# vtk DataFile Version 3.0\n";
  out << "PoU scalar field (point-centered)\n";
  out << "ASCII\n";
  out << "DATASET STRUCTURED_POINTS\n";
  out << "DIMENSIONS " << nx_pts << " " << ny_pts << " 1\n";
  out << "ORIGIN " << std::setprecision(16) << xmin << " " << ymin << " 0\n";
  out << "SPACING " << dx << " " << dy << " 1\n";
  out << "POINT_DATA " << nx_pts * ny_pts << "\n";
  out << "SCALARS F double\n";
  out << "LOOKUP_TABLE default\n";
  for (int j = 0; j < ny_pts; ++j) {
    for (int i = 0; i < nx_pts; ++i) {
      out << std::setprecision(16) << Fpoints[j * nx_pts + i] << "\n";
    }
  }
  out.close();
}

void dumpPoUFieldVTK(const std::string& filename,
                     const Data<IRL2D::Parabola>& interface,
                     const Data<IRL2D::Moments>& liquid_moments,
                     const BasicMesh& coarse_mesh, const BasicMesh& fine_mesh,
                     IRL2D::Vec minPt, IRL2D::Vec maxPt) {
  std::vector<IRL2D::Vec> centroids, datums, normals, tangents;
  std::vector<double> a_list;

  for (int i = coarse_mesh.imin(); i <= coarse_mesh.imax(); ++i) {
    for (int j = coarse_mesh.jmin(); j <= coarse_mesh.jmax(); ++j) {
      const double lvf =
          (liquid_moments)(i, j).m0() / coarse_mesh.cell_volume();
      if (lvf >= IRL::global_constants::VF_LOW &&
          lvf <= IRL::global_constants::VF_HIGH) {
        IRL2D::BezierList cell = IRL2D::RectangleFromBounds(
            IRL2D::Vec(coarse_mesh.x(i), coarse_mesh.y(j)),
            IRL2D::Vec(coarse_mesh.x(i + 1), coarse_mesh.y(j + 1)));
        IRL2D::BezierList clipped_interface =
            IRL2D::ParabolaClip(cell, interface(i, j), true);
        IRL2D::Vec start = clipped_interface[0].first;
        IRL2D::Vec end = clipped_interface[1].first;
        IRL2D::Vec centroid = IRL2D::centroidParabolaSegmentAnalytic(
            interface(i, j).datum(), interface(i, j).frame()[0],
            interface(i, j).frame()[1], interface(i, j).coeff(), start, end);
        centroids.push_back(centroid);
        datums.push_back(interface(i, j).datum());
        normals.push_back(interface(i, j).frame()[1]);
        tangents.push_back(interface(i, j).frame()[0]);
        a_list.push_back(interface(i, j).coeff());
      }
    }
  }

  IRL2D::ImplicitSurfaceParabola IS(datums, tangents, normals, a_list,
                                    centroids, 2.5 * coarse_mesh.dx());
  std::vector<double> Fpoints;
  for (int j = fine_mesh.jmin(); j <= (fine_mesh.jmax() + 1); j++) {
    for (int i = fine_mesh.imin(); i <= (fine_mesh.imax() + 1); i++) {
      double Fval = IS.F(IRL2D::Vec(fine_mesh.x(i), fine_mesh.y(j)));
      if (std::isnan(Fval)) {
        Fval = 10000.;
      }
      double mag = std::sqrt(fine_mesh.xm(i) * fine_mesh.xm(i) +
                             fine_mesh.ym(j) * fine_mesh.ym(j));
      // if (mag < 0.13) Fval = 10000.;
      Fpoints.push_back(Fval);
    }
  }
  const int nx_pts = fine_mesh.getNx() + 1;
  const int ny_pts = nx_pts;

  // write vtk
  writeVTKStructuredPoints_PointData(filename, nx_pts, ny_pts, minPt[0],
                                     maxPt[0], minPt[1], maxPt[1], Fpoints);
}
