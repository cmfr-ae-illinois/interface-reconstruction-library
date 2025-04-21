#include <string>
#include <cstdlib>
#include <iostream>
#include <chrono>

template <class Data>
void WriteField(const int ncells, const std::vector<double>& coords,
                const Data& field, const std::string& field_name) {
  const auto file_name = "./" + field_name + ".vtr";
  std::cout<<"in write field"<<std::endl;
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
  std::cout<<"loop start";
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
    triangulated_surface[i] = surfaces[i].triangulate(1.0e-1, 5);
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