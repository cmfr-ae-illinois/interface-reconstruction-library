#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <random>

#include "irl/ml_classification/data_gen.h"

//#include <eigen3> do not need to include? error when including, why?

//#include "irl/generic_cutting/generic_cutting.h"
/*
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

// Function to generate a random vector within given bounds
Eigen::Vector3d generateRandomPoint(double min_bound, double max_bound, std::default_random_engine& eng) {
    std::uniform_real_distribution<double> dist(min_bound, max_bound);
    return Eigen::Vector3d(dist(eng), dist(eng), dist(eng));  // Random point in 3D space
}
*/

int main(int argc, char* argv[]) {
  /*
  // Define bounds for the random points (you can change these values)
  double min_bound = -1.0;  // Minimum bound for each coordinate
  double max_bound = 1.0;   // Maximum bound for each coordinate

  // Initialize random number generator
  std::random_device rd;
  std::default_random_engine eng(rd());

  // Generate random points for axis_point1 and axis_point2
  Eigen::Vector3d axis_point1 = generateRandomPoint(min_bound, max_bound, eng);
  Eigen::Vector3d axis_point2 = generateRandomPoint(min_bound, max_bound, eng);

  Eigen::Vector3d axis_direction = axis_point2-axis_point1;
  double radius = 0.7;
  IRL::Pt datum(0.1, -0.3, 0.2);

  // Create nxnxn vector of volume fractions and surfaces
  int stencil_size = 9;
  std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));
  
  // Defining cell coordinates
  auto coords = std::vector<double>(stencil_size+1);
  for (int i = 0; i < stencil_size+1; i++) {
    coords[i] = -0.5*stencil_size + static_cast<double>(i);
  }
  const double cell_volume = 1.0;

  // Bounds of paraboloid parameters
  //std::uniform_real_distribution<double> random_rotation(-0.5 * M_PI,
  //                                                       0.5 * M_PI);
  //std::uniform_real_distribution<double> random_coeff(-5.0, 5.0);
  //std::uniform_real_distribution<double> random_translation(-0.5, 0.5);

  // Create reference frame
  //const auto frame = IRL::ReferenceFrame(
  //    IRL::Normal(1, 0, 0), IRL::Normal(0, 1, 0), IRL::Normal(0, 0, 1));
  // Create random rotation angles, datum and coefficiens
  //double angles[3] = {random_rotation(eng), random_rotation(eng), 0.0};
  //IRL::Pt datum(random_translation(eng), random_translation(eng),
                //random_translation(eng));
  //double coeffs[2] = {random_coeff(eng), random_coeff(eng)};
  // Generate rotation matrices
  //IRL::UnitQuaternion x_rot(angles[0], frame[0]);
  //IRL::UnitQuaternion y_rot(angles[1], frame[1]);
  //IRL::UnitQuaternion z_rot(angles[2], frame[2]);
  // Create random paraboloid
  //const auto paraboloid = IRL::Paraboloid(datum, x_rot * y_rot * z_rot * frame,
                                          //coeffs[0], coeffs[1]);

  // Initialize field
  std::vector<IRL::ParametrizedSurfaceOutput> surfaces;
  for (int i = 0; i < stencil_size; i++) {
    for (int j = 0; j < stencil_size; j++) {
      for (int k = 0; k < stencil_size; k++) {
        // Create cell
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(coords[i], coords[j], coords[k]),
            IRL::Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

        // Create paraboloid for this cell in the steps below
        // Find center point of cell
        Eigen::Vector3d cell_center((coords[i+1]+coords[i])/2.0,
                                    (coords[j+1]+coords[j])/2.0,
                                    (coords[k+1]+coords[k])/2.0);
        Eigen::Vector3d axis_point1_to_cell_center = cell_center - axis_point1;
        Eigen::Vector3d axis_point2_to_cell_center = cell_center - axis_point2;
        // Calculate the projection of the cell center onto the cylinder axis
        double projection_factor = axis_point1_to_cell_center.dot(axis_direction) / axis_direction.squaredNorm();
        Eigen::Vector3d closest_point_on_axis = axis_point1 + projection_factor * axis_direction;
        // Find datum of paraboloid
        Eigen::Vector3d datum_paraboloid_eVec = closest_point_on_axis + radius * (cell_center - closest_point_on_axis).normalized();
        IRL::Pt datum_paraboloid(datum_paraboloid_eVec.x(), datum_paraboloid_eVec.y(), datum_paraboloid_eVec.z());

        // Create frame
        Eigen::Vector3d paraboloid_x = axis_direction.normalized();
        Eigen::Vector3d paraboloid_z = cell_center-closest_point_on_axis;
        paraboloid_z.normalize();
        Eigen::Vector3d paraboloid_y = paraboloid_z.cross(paraboloid_x);
        paraboloid_y.normalize();
        //frame[1] = IRL::crossProduct(frame[2],frame[0]); frame2 is normal, frame0 is axis
        //frame[1].normalize();

        const auto frame = IRL::ReferenceFrame(IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                                              IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                                              IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));
        // Create paraboloid
        const auto paraboloid = IRL::Paraboloid(datum_paraboloid, frame, 0, 1/(2*radius));

        // Intersect cell with paraboloid -- return volume and surface
        auto volume_fraction = IRL::getVolumeFraction(cell, paraboloid);
        std::cout << "VFRAC(" << i << ", " << j << ", " << k
                  << ") = " << volume_fraction << std::endl;
        std::cout << "vfracs generated" << std::endl;
        auto volume_and_surface = IRL::getVolumeMoments<
            IRL::AddSurfaceOutput<IRL::Volume, IRL::ParametrizedSurfaceOutput>>(
            cell, paraboloid);
        // Store surface and volume (fraction)
        
        vfrac[i][j][k] = volume_and_surface.getMoments().volume() / cell_volume;
        //if vfrac[i][j][k]<0.9999999 {
          surfaces.push_back(volume_and_surface.getSurface());
        //}
      }
    }
  }
  // Print field to file
  WriteField(stencil_size, coords, vfrac, "vfrac");
  WriteSurface(surfaces, "surface");
  */
  int stencil_size = 3;
  IRL::Data_gen data_gen;
  //std::vector<double> state = data_gen.generate_Paraboloid(stencil_size, 1, 0.05, true);
  std::vector<double> state = data_gen.generate_Sheet(stencil_size, true);
  //std::vector<double> state = data_gen.generate_Cylinder(stencil_size, true);
  //std::vector<double> state = data_gen.generate_Sphere(stencil_size, true);
  return 0;
}