#include <chrono>
#include <ctime>
#include <iostream>
#include <string>

#include "irl/generic_cutting/generic_cutting.h"

template <class Data>
void WriteField(const int ncells, const std::vector<double>& coords,
                const Data& field, const std::string& field_name) {
  const auto file_name = "./" + field_name + ".vtr";

  FILE* file;
  file = fopen(file_name.c_str(), "w");

  fprintf(file, "<VTKFile type=\"RectilinearGrid\">\n");
  fprintf(file, "<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0,
          ncells, 0, ncells, 0, ncells);
  fprintf(file, "<Piece Extent=\"%d %d %d %d %d %d\">\n", 0, ncells, 0, ncells,
          0, ncells);

  fprintf(file, "<Coordinates>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < ncells + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(coords[i]));
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < ncells + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(coords[i]));
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < ncells + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(coords[i]));
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file, "</Coordinates>\n");

  fprintf(file, "<PointData>\n</PointData>\n");

  fprintf(file, "<CellData Scalars=\"");
  fprintf(file, "%s ", field_name.c_str());
  fprintf(file, "\" >\n");
  fprintf(file,
          "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n",
          field_name.c_str());
  for (int k = 0; k < ncells; ++k) {
    for (int j = 0; j < ncells; ++j) {
      for (int i = 0; i < ncells; ++i) {
        fprintf(file, "%15.8E ", static_cast<double>(field[i][j][k]));
      }
    }
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n</RectilinearGrid>\n</VTKFile>\n");
  fclose(file);
}

void WriteSurface(const std::vector<IRL::ParametrizedSurfaceOutput>& surfaces,
                  const std::string& file_subname) {
  const auto file_name = "./" + file_subname + ".vtu";
  std::vector<IRL::TriangulatedSurfaceOutput> triangulated_surface;
  triangulated_surface.resize(surfaces.size());
  for (std::size_t i = 0; i < surfaces.size(); ++i) {
    triangulated_surface[i] = surfaces[i].triangulate(1.0e-2, 5);
  }

  int number_of_vertices = 0;
  std::vector<int> offset(triangulated_surface.size() + 1, 0);
  for (int i = 0; i < triangulated_surface.size(); ++i) {
    const auto& vlist = triangulated_surface[i].getVertexList();
    number_of_vertices += vlist.size();
    offset[i + 1] = offset[i] + vlist.size();
  }

  int number_of_triangles = 0;
  for (int i = 0; i < triangulated_surface.size(); ++i) {
    const auto& tlist = triangulated_surface[i].getTriangleList();
    number_of_triangles += tlist.size();
  }

  FILE* file;
  file = fopen(file_name.c_str(), "w");
  fprintf(file, "<VTKFile type=\"UnstructuredGrid\">\n");
  fprintf(file, "<UnstructuredGrid>\n");
  fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          number_of_vertices, number_of_triangles);
  fprintf(file, "<Points>\n");
  fprintf(file, "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n");
  for (std::size_t i = 0; i < triangulated_surface.size(); ++i) {
    const auto& vlist = triangulated_surface[i].getVertexList();
    for (const auto& vertex : vlist) {
      fprintf(file, "%15.8E %15.8E %15.8E ", static_cast<double>(vertex[0]),
              static_cast<double>(vertex[1]), static_cast<double>(vertex[2]));
    }
  }
  fprintf(file, "</DataArray>\n</Points>\n");

  fprintf(file, "<Cells>\n");
  fprintf(file,
          "<DataArray type=\"Int32\" Name=\"connectivity\" "
          "format=\"ascii\">\n");
  for (int i = 0; i < triangulated_surface.size(); ++i) {
    const auto& tlist = triangulated_surface[i].getTriangleList();
    const auto off = offset[i];
    for (const auto& triangle : tlist) {
      const auto& index_mapping = triangle.getIndexMapping();
      fprintf(file, "%d %d %d ", static_cast<int>(off + index_mapping[0]),
              static_cast<int>(off + index_mapping[1]),
              static_cast<int>(off + index_mapping[2]));
    }
  }
  fprintf(file, "</DataArray>\n");

  fprintf(file,
          "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for (std::size_t i = 0; i < number_of_triangles; ++i) {
    fprintf(file, "%d ", static_cast<int>(3 * (i + 1)));
  }
  fprintf(file, "</DataArray>\n");

  fprintf(file, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (std::size_t i = 0; i < number_of_triangles; ++i) {
    fprintf(file, "5 ");
  }
  fprintf(file, "</DataArray>\n");

  fprintf(file, "</Cells>\n");
  fprintf(file, "<CellData>\n");
  fprintf(file, "<DataArray type=\"Int32\" Name=\"ID\" format=\"ascii\">\n");
  for (std::size_t i = 0; i < triangulated_surface.size(); ++i) {
    const auto& tlist = triangulated_surface[i].getTriangleList();
    for (std::size_t j = 0; j < tlist.size(); ++j) {
      fprintf(file, "%d ", static_cast<int>(i));
    }
  }
  fprintf(file, "</DataArray>\n");
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
  fclose(file);
}

int main(int argc, char* argv[]) {
  // Create 5x5x5 array of volume fractions and surfaces
  std::array<std::array<std::array<double, 5>, 5>, 5> vfrac;

  // Defining cell coordinates
  auto coords = std::vector<double>(6);
  for (int i = 0; i < 6; i++) {
    coords[i] = -2.5 + static_cast<double>(i);
  }
  const double cell_volume = 1.0;

  // Create random number generator and seed it with entropy
  std::random_device rd;
  std::mt19937_64 eng(static_cast<int>(std::time(0)));

  // Bounds of paraboloid parameters
  std::uniform_real_distribution<double> random_rotation(-0.5 * M_PI,
                                                         0.5 * M_PI);
  std::uniform_real_distribution<double> random_coeff(-5.0, 5.0);
  std::uniform_real_distribution<double> random_translation(-0.5, 0.5);

  // Create reference frame
  const auto frame = IRL::ReferenceFrame(
      IRL::Normal(1, 0, 0), IRL::Normal(0, 1, 0), IRL::Normal(0, 0, 1));
  // Create random rotation angles, datum and coefficiens
  double angles[3] = {random_rotation(eng), random_rotation(eng), 0.0};
  IRL::Pt datum(random_translation(eng), random_translation(eng),
                random_translation(eng));
  double coeffs[2] = {random_coeff(eng), random_coeff(eng)};
  // Generate rotation matrices
  IRL::UnitQuaternion x_rot(angles[0], frame[0]);
  IRL::UnitQuaternion y_rot(angles[1], frame[1]);
  IRL::UnitQuaternion z_rot(angles[2], frame[2]);
  // Create random paraboloid
  const auto paraboloid = IRL::Paraboloid(datum, x_rot * y_rot * z_rot * frame,
                                          coeffs[0], coeffs[1]);

  // Initialize field
  std::vector<IRL::ParametrizedSurfaceOutput> surfaces;
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      for (int k = 0; k < 5; k++) {
        // Create cell
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(coords[i], coords[j], coords[k]),
            IRL::Pt(coords[i + 1], coords[j + 1], coords[k + 1]));
        // Intersect cell with paraboloid -- return volume and surface
        auto volume_fraction = IRL::getVolumeFraction(cell, paraboloid);
        std::cout << "VFRAC(" << i << ", " << j << ", " << k
                  << ") = " << volume_fraction << std::endl;
        auto volume_and_surface = IRL::getVolumeMoments<
            IRL::AddSurfaceOutput<IRL::Volume, IRL::ParametrizedSurfaceOutput>>(
            cell, paraboloid);
        // Store surface and volume (fraction)
        surfaces.push_back(volume_and_surface.getSurface());
        vfrac[i][j][k] = volume_and_surface.getMoments().volume() / cell_volume;
      }
    }
  }

  // Print field to file
  WriteField(5, coords, vfrac, "vfrac");
  WriteSurface(surfaces, "surface");

  return 0;
}