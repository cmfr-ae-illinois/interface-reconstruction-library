#include "irl/geometry/general/new_pt_calculation_functors.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/rotations.h"

#include <sys/time.h>
#include <cmath>
#include <random>
#include <chrono>

#include <omp.h>
#include <getopt.h>

#include "irl/data_structures/small_vector.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"

#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection.h"
#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection_amr.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron_quadratic.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron_quadratic.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
// #include
// "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/moments/general_moments.h"

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

#include "irl/planar_reconstruction/planar_separator.h"

#include "tools/src/cylinder_testing/geometry.h"

using namespace IRL;

int first_vertex_on_surface;
static int polyhedron = 0;
static int old_polyhedron = 0;
int re_scale;

static struct option long_options[] =
    {
        {"number_of_configuration", required_argument, 0, 'n'},
        {"AMR_level",               required_argument, 0, 'l'},
        {"elliptic",                no_argument,       0, 'e'},
        {"hyperbolic",              no_argument,       0, 'h'},
        {"both",                    no_argument,       0, 'b'},
        {"seed",                    required_argument, 0, 's'},
        {"file_name",               required_argument, 0, 'o'},
        {"graded_gen",              no_argument,       0, 'g'},
        {"help",                    no_argument,       0, 'H'},
        {"first_vertex_on_surface", no_argument, &first_vertex_on_surface, 1},
        {"no_rescale",              no_argument, &re_scale, 0},
        {"force_rescale",           no_argument, &re_scale, 1},
        {"cube",                    no_argument, &polyhedron, 0},
        {"old_cube",                no_argument, &old_polyhedron, 1},
        {"tet",                     no_argument, &polyhedron, 1},
        {"dodecahedron",            no_argument, &polyhedron, 2},
        {"cube_hole",               no_argument, &polyhedron, 3},
        {"cube_hole_convex",        no_argument, &polyhedron, 4},
        {"bunny",                   no_argument, &polyhedron, 5},
        {0, 0, 0, 0}
    };

static const std::array<double, 9> gratted_b({0.9, 1.0, 16.0 / 9.0, 2.0, 9.0 / 4.0, 4, -0.75, -1.0, -1.25});
static const std::array<double, 5> gratted_r({0.0625, 0.25, 0.5, 9.0 / 16.0, 1.0});
static const std::array<double, 5> gratted_rotation({-M_PI, -M_PI / 2.0, 0.0, M_PI / 2.0, M_PI});
static const std::array<double, 5> gratted_translation({- 0.5, - 0.25, 0.0, 0.25, 0.5});

std::array<double, 8> get_config(int i) {
  int index = i;
  int b_index = index % 9;
  index /= 9;
  int r_index = index % 5;
  index /= 5;
  int rotation_index_1 = index % 5;
  index /= 5;
  int rotation_index_2 = index % 5;
  index /= 5;
  int rotation_index_3 = index % 5;
  index /= 5;
  int translation_index_1 = index % 5;
  index /= 5;
  int translation_index_2 = index % 5;
  index /= 5;
  int translation_index_3 = index % 5;
  return {gratted_b[b_index], gratted_r[r_index], gratted_rotation[rotation_index_1], 
      gratted_rotation[rotation_index_2], gratted_rotation[rotation_index_3], gratted_translation[translation_index_1],
      gratted_translation[translation_index_2], gratted_translation[translation_index_3]};
}

void print_help(std::string prg_name) {
  std::cout << "Usage: " << prg_name << " [options]\n\n";
  std::cout << "Options:\n";
  std::cout << "  -H, --help                             Show this help message\n";
  std::cout << "  -n, --number_of_configuration <int>    Number of configurations (default: 1000)\n";
  std::cout << "  -l, --AMR_level <int>                  AMR level (default: 15)\n";
  std::cout << "  -s, --seed <int>                       Random seed (default: current time)\n";
  std::cout << "  -o, --file_name <string>               Output file name (default: result.csv)\n";
  std::cout << "  -e, --elliptic                         Use elliptic bounds (min_b=0, max_b=10)\n";
  std::cout << "  -h, --hyperbolic                       Use hyperbolic bounds (min_b=-10, max_b=0)\n";
  std::cout << "  -b, --both                             Use full bounds (min_b=-10, max_b=10)\n";
  std::cout << "  -g, --graded_gen                       Use graded generation mode\n";
  std::cout << "      --first_vertex_on_surface          Start with first vertex on surface\n";
  std::cout << "      --no_rescale                       Disable rescaling\n";
  std::cout << "      --force_rescale                    Force rescaling\n";
  std::cout << "      --cube                             Use cube polyhedron (default)\n";
  std::cout << "      --old_cube                         Use old cube (center is 0.5 0.5 0.5 and is not rescaled)\n";
  std::cout << "      --tet                              Use tetrahedron\n";
  std::cout << "      --dodecahedron                     Use dodecahedron\n";
  std::cout << "      --cube_hole                        Use cube with hole\n";
  std::cout << "      --cube_hole_convex                 Use convex cube with hole\n";
  std::cout << "      --bunny                            Use bunny shape\n";
}


int main(int argc, char* argv[]) {
  int Ntests = 1000;
  int seed = std::chrono::system_clock::now().time_since_epoch().count();
  float min_b = -3.0;
  float max_b = 5.0;
  float min_r = 0.0;
  float max_r = 5.0;
  int long_id = 0;
  int AMR_levels = 15;
  re_scale = 1;
  first_vertex_on_surface = 0;
  bool graded_gen = false;
  std::string filename = "result.csv";

  while(1) {
      int result = getopt_long(argc, argv, "l:n:s:ehbo:gH", long_options, &long_id);
      if (result == -1)
      {
          break;
      }
      switch (result)
      {
      case 'l':
          AMR_levels = std::atoi(optarg);
          break;
      case 'n':
          Ntests = std::atoi(optarg);
          break;
      case 's':
          seed = std::atoi(optarg);
          break;
      case 'o':
          filename = optarg;
          break;
      case 'e':
          min_b = 0.0;
          max_b = 10.0;
          break;
      case 'h':
          min_b = -10.0;
          max_b = 0.0;
          break;
      case 'b':
          min_b = -10.0;
          max_b = 10.0;
          break;
      
      case 'g':
          graded_gen = true;
          Ntests = 9*5*5*5*5*5*5*5;
          re_scale = 0;
          break;

      case 'H': // help
          print_help(argv[0]);
          exit(0);
      case '?':
          std::cerr << "Invalid argument detected.\n\n";
          print_help(argv[0]);
          exit(1);
      
      case 0:
          break;
      
      default:
          std::cout << "error while reading arguments : " << result << std::endl;
          return 1;
      }
  }

  int geom = re_scale ? polyhedron : - polyhedron;

  #pragma omp parallel
  {
    #pragma omp single
    {
      int polytype = old_polyhedron ? -1 : polyhedron;
      std::cout << "running with " << omp_get_max_threads() << " threads (seed is :" << seed << ") polygon type : " << polytype << std::endl;
      std::ofstream file(filename);
  
      // Write the CSV header
      file << "b,r,"
          << "rotation_1,rotation_2,rotation_3,"
          << "datum_x,datum_y,datum_z,"
          << "is_on_surface,"
          << "geometry,"
          << "Volume,"
          << "Centroid_x,Centroid_y,Centroid_z"
          << std::endl;
      file.close();
    }

    #pragma omp barrier
    int thread_id = omp_get_thread_num(); // Get thread ID

    // Create random number generator and seed it with entropy
    std::random_device rd;
    std::mt19937_64 eng(seed + thread_id);  // rd());

    // Bounds of paraboloid parameters
    std::uniform_real_distribution<double> random_rotation(-M_PI,
                                                            M_PI);
    std::uniform_real_distribution<double> random_coeffs_b(min_b, max_b);
    std::uniform_real_distribution<double> random_coeffs_r(0.0, 1.2);
    std::uniform_real_distribution<double> random_translation(-0.5, 0.5);

    auto start = std::chrono::high_resolution_clock::now();

    int count = 0;

    #pragma omp for nowait schedule(dynamic)
    for (int i = 0; i < Ntests; i++) {

      double b, r, rot_1, rot_2, rot_3, tra_1, tra_2, tra_3;

      if (graded_gen) {
        auto config = get_config(i);
        b = config[0];
        r = config[1];
        rot_1 = config[2];
        rot_2 = config[3];
        rot_3 = config[4];
        tra_1 = config[5];
        tra_2 = config[6];
        tra_3 = config[7];
      }
      else {
        b = random_coeffs_b(eng);
        r = random_coeffs_r(eng);
        rot_1 = random_rotation(eng);
        rot_2 = random_rotation(eng);
        rot_3 = random_rotation(eng);
        tra_1 = random_translation(eng);
        tra_2 = random_translation(eng);
        tra_3 = random_translation(eng);
      }

      auto aligned_cylinder =
          AlignedCylinder({b, r});

      double angles[3] = {rot_1, rot_2, rot_3};
      Pt datum(tra_1, tra_2, tra_3);

      ReferenceFrame frame(Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0),
                          Normal(0.0, 0.0, 1.0));

      UnitQuaternion x_rotation(angles[0], frame[0]);
      UnitQuaternion y_rotation(angles[1], frame[1]);
      UnitQuaternion z_rotation(angles[2], frame[2]);
      frame = x_rotation * y_rotation * z_rotation * frame;

      HalfEdgePolyhedronQuadratic<Pt> half_edge;

      if (old_polyhedron) {
        // Defining the square
        RectangularCuboid cube = RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0),
                                                Pt(1.0, 1.0, 1.0));



        int pt_above = 0, pt_below = 0;
        for (auto& vertex : cube) {
          Pt tmp_pt = vertex - datum;
          for (UnsignedIndex_t d = 0; d < 3; ++d) {
            vertex[d] = frame[d] * tmp_pt;
          }
        }
        
        if (first_vertex_on_surface)
        {
          Pt clip_translation = Pt(0.0, 0.0, sqrt(r))-cube[0];
          for (auto& vertex : cube) {
            Pt tmp_pt = vertex + clip_translation;
            for (int d = 0; d < 3; d++) {
                vertex[d] = tmp_pt[d];
            }
          }
        }

        cube.setHalfEdgeVersion(&half_edge);
      } else {
        auto geometri_and_connectivity_amr = getGeometry(polyhedron, re_scale);
        auto geometri_amr = geometri_and_connectivity_amr.first;
        auto connectivity_amr = geometri_and_connectivity_amr.second;
        for (auto& vertex : geometri_amr) {
            Pt tmp_pt = vertex - datum;
            for (UnsignedIndex_t d = 0; d < 3; ++d) {
                vertex[d] = frame[d] * tmp_pt;
            }
        }
        if (first_vertex_on_surface)
        {
          Pt clip_translation = Pt(0.0, 0.0, sqrt(r))-geometri_amr[0];
          for (auto& vertex : geometri_amr) {
            Pt tmp_pt = vertex + clip_translation;
            vertex = tmp_pt;
          }
        }
        
        geometri_amr.setHalfEdgeVersion(&half_edge);

      }
      auto seg_half_edge = half_edge.generateSegmentedPolyhedron();

      // Calculate volume of unclipped dodecahedron using AMR
      VolumeMoments amr_volume_moments =
          intersectPolyhedronWithCylinderAMR<VolumeMoments>(
              &seg_half_edge, &half_edge, aligned_cylinder, AMR_levels);

      #pragma omp critical
      {
        std::ofstream file(filename, std::ios::app);
        file << std::setprecision(20);
        file << b << "," << r << ",";
        file << angles[0] << "," << angles[1] << ","<< angles[2] << ",";
        file << datum[0] << "," << datum[1] << "," << datum[2] << ",";
        file << first_vertex_on_surface << ",";
        file << geom << ",";
        file << amr_volume_moments.volume() << ",";
        auto amr_Centroid = amr_volume_moments.centroid();
        file << amr_Centroid[0] << "," << amr_Centroid[1] << "," << amr_Centroid[2] << std::endl;
        file.close();
      }

      count++;
    }
    auto end = std::chrono::high_resolution_clock::now();



    auto duration = std::chrono::duration<double>(end - start).count();

    auto s = fmod(duration, 60.0);
    auto m = fmod((duration - s) / 60.0, 60.0);
    auto h = fmod(((duration - s) / 60.0 - m) / 60.0, 24.0);
    auto j = (((duration - s) / 60.0 - m) / 60.0 - h) / 24.0;

    #pragma omp critical
    {
      std::cout << "Results written to " << filename << ", time taken to do " << count << "cases : " <<
      j << ":" << h << ":" << m << ":" << s << std::endl;
    }
  }
  return 0;
}