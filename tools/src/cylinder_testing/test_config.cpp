#include "irl/geometry/general/new_pt_calculation_functors.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/rotations.h"

#include <sys/time.h>
#include <cmath>
#include <random>
#include <chrono>

#include <omp.h>
#include <stdlib.h>
#include <getopt.h>

// #define VALDEBUG_S
// #define DEBUG_CYL_IRL

// SURFACE_DEBUG add log about the cylinder output surface, including:
// when new arcs are added and the computation of the closed loops
// #define VALDEBUG2
// #define NUDGE_REGION

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

static int no_surface_output = 0;
int first_vertex_on_surface;
static int old_format = 0;
static int polyhedron = 0;
static int old_polyhedron = 0;
int re_scale;

void print_time(std::chrono::_V2::system_clock::time_point start, 
                std::chrono::_V2::system_clock::time_point end) {

    auto duration = std::chrono::duration<double>(end - start).count();

    auto s = fmod(duration, 60.0);
    auto m = fmod((duration - s) / 60.0, 60.0);
    auto h = fmod(((duration - s) / 60.0 - m) / 60.0, 24.0);
    auto j = (((duration - s) / 60.0 - m) / 60.0 - h) / 24.0;

    std::cout <<
      j << ":" << h << ":" << m << ":" << s;
}

static struct option long_options[] =
    {
        {"no_surface_output",   no_argument, &no_surface_output, 1},
        {"first_vertex_on_surface",no_argument, &first_vertex_on_surface, 1},
        {"old_format",          no_argument, &old_format, 1},
        {"cube",                no_argument, &polyhedron, 0},
        {"old_cube",            no_argument, &old_polyhedron, 1},
        {"tet",                 no_argument, &polyhedron, 1},
        {"dodecahedron",        no_argument, &polyhedron, 2},
        {"cube_hole",           no_argument, &polyhedron, 3},
        {"cube_hole_convex",    no_argument, &polyhedron, 4},
        {"bunny",               no_argument, &polyhedron, 5},
        {"no_rescale",              no_argument, &re_scale, 0},
        {"force_rescale",           no_argument, &re_scale, 1},
        {"AMR_level",           required_argument,       0, 'l'},
        {"max_AMR_level",       required_argument,       0, 'l'},
        {"triangluation_Level", required_argument,       0, 'L'},
        {"file",                required_argument,       0, 'f'},
        {"random_position",     optional_argument,       0, 'r'},
        {"help",                    no_argument,       0,   'H'},
        {0, 0, 0, 0}
    };

void print_help(const char* program_name) {
    std::cout << "Usage:\n";
    std::cout << "  " << program_name << " [OPTIONS] -- b r datum_x datum_y datum_z rotation_x rotation_y rotation_z [first_vertex_on_surface]\n";
    std::cout << "  " << program_name << " [OPTIONS] -- b,r,rotation_x,rotation_y,rotation_z,datum_x,datum_y,datum_z[,first_vertex_on_surface]\n";
    std::cout << "  " << program_name << " [OPTIONS] --file <path>\n";
    std::cout << "  " << program_name << " [OPTIONS] -r[<seed>]\n\n";

    std::cout << "Options:\n";
    std::cout << "  -H, --help                                         Show this help message\n";
    std::cout << "  -l, --AMR_level <int> / --max_AMR_level <int>      Set the maximum AMR level (default: 10)\n";
    std::cout << "  -f, --file <string>                                Input file path\n";
    std::cout << "  -r, --random_position[=<seed>]                     Generate random coordinates; optional seed\n";
    std::cout << "  -L, --triangluation_Level                          Set the max distance used for triangulating the output surface (default : 0.1)\n";
    std::cout << "      --no_surface_output                            Disable surface output\n";
    std::cout << "      --first_vertex_on_surface                      Place first vertex on surface\n";
    std::cout << "      --no_rescale                                   Disable rescaling\n";
    std::cout << "      --force_rescale                                Force rescaling\n";
    std::cout << "      --old_format                                   Use old input format (remove rotation_z from input)\n";
    std::cout << "      --cube                                         Use cube polyhedron\n";
    std::cout << "      --old_cube                                     Use old cube\n";
    std::cout << "      --tet                                          Use tetrahedron\n";
    std::cout << "      --dodecahedron                                 Use dodecahedron\n";
    std::cout << "      --cube_hole                                    Use cube with a hole\n";
    std::cout << "      --cube_hole_convex                             Use convex cube with a hole\n";
    std::cout << "      --bunny                                        Use bunny mesh\n";
}

template <class HalfEdgeType, class SegmentedHalfEdge, class CuttingType>
void procede_case(Cylinder a_cylinder, CuttingType irl_cutting, SegmentedHalfEdge a_segmented_half_edge, HalfEdgeType a_half_edge, int max_amr_level, float tri_level) {
    using VolumeAndSuface =
        AddSurfaceOutput<VolumeMoments, CylinderParametrizedSurfaceOutput>;

    if (no_surface_output == 0) {
        std::ofstream file("./full_geometry_test_config.vtu");
        file << a_segmented_half_edge;
        file.close();
    }

    auto start = std::chrono::high_resolution_clock::now();
    auto temp_surface_and_moments =
        getVolumeMoments<VolumeAndSuface>(irl_cutting, a_cylinder);
    auto end = std::chrono::high_resolution_clock::now();
    auto temp_moments =
        getVolumeMoments<VolumeMoments>(irl_cutting, a_cylinder);
    auto start_amr = std::chrono::high_resolution_clock::now();
    std::string output_file = no_surface_output ? IRL::no_amr_output : "_test_config_cube";
    #ifdef DEBUG_CYL_IRL
    std::cout << "computing AMR" << std::endl;
    #endif
    auto amr_moments =
        intersectPolyhedronWithCylinderAMR<VolumeMoments>(
            &(a_segmented_half_edge), &(a_half_edge), a_cylinder.getAlignedCylinder(), max_amr_level, output_file);
    auto end_amr = std::chrono::high_resolution_clock::now();
    std::cout << "the computed volume (with surface) is :" << std::setprecision(20)
              << temp_surface_and_moments.getMoments().volume().volume()
              << std::endl;
    std::cout << "the computed volume (w/o  surface) is :" << std::setprecision(20)
              << temp_moments.volume().volume()
              << std::endl;
    std::cout << "the amr      volume                is :"
              << amr_moments.volume().volume()
              << std::endl;
    std::cout << "the computed centroid (with surface) is :"
              << temp_surface_and_moments.getMoments().centroid() << std::endl;
    std::cout << "the computed centroid (w/o  surface) is :"
              << temp_moments.centroid() << std::endl;
    std::cout << "the amr      centroid                is :"
              << amr_moments.centroid() << std::endl;

    auto amr_Centroid = amr_moments.centroid().getPt();
    auto Centroid = temp_moments.centroid().getPt();
    auto CentroidS = temp_surface_and_moments.getMoments().centroid().getPt();

    auto error_centroid = Centroid - amr_Centroid;
    auto error_centroidS = CentroidS - amr_Centroid;
    auto& centroid = temp_surface_and_moments.getMoments().centroid().getPt();
    centroid /= temp_surface_and_moments.getMoments().volume().volume();
    auto& centroid_2 = temp_moments.centroid().getPt();
    centroid_2 /= temp_moments.volume().volume();
    auto& amr_centroid = amr_moments.centroid().getPt();
    amr_centroid /= amr_moments.volume().volume();
    std::cout << "the computed center of mass (with surface) is :" << centroid
              << std::endl;
    std::cout << "the computed center of mass (w/o  surface) is :" << centroid_2
              << std::endl;
    std::cout << "the amr      center of mass                is :" << amr_centroid
              << std::endl << std::endl;

    std::cout << "time taken to compute volume : ";
    print_time(start, end);
    std::cout << "\ntime taken to compute volume with amr : ";
    print_time(start_amr, end_amr);
    std::cout << std::endl;

    auto amr_volum = amr_moments.volume().volume();
    auto volume = temp_moments.volume().volume();
    auto volumeS = temp_surface_and_moments.getMoments().volume().volume();
    auto volume_error = abs(volume-amr_volum);
    auto volume_errorS = abs(volumeS-amr_volum);

    std::cout << "errors :\n" << volume_error
              << "," << error_centroid[0]
              << "," << error_centroid[1]
              << "," << error_centroid[2] << std::endl;
    std::cout << "errors (surface ):\n" << volume_errorS
              << "," << error_centroidS[0]
              << "," << error_centroidS[1]
              << "," << error_centroidS[2] << std::endl;

    if (no_surface_output == 0) {
        auto temp_param_surface = temp_surface_and_moments.getSurface();
        auto temp_tri_surface = temp_param_surface.triangulate(tri_level);
        temp_tri_surface.write("test_config_surface");
    }
}

struct DataEntry {
    double b, r, rotation_1, rotation_2, rotation_3, datum_x, datum_y, datum_z, is_on_surface, geometry;
};

void read_file(char* filepath, int max_amr_level, float tri_level) {
    std::ifstream file(filepath);
    std::vector<DataEntry> datas;
    std::string line;
    int totalEntries = 0;
    
    if (!file.is_open()) {
        std::cerr << "Error opening file!\n";
        return;
    }
    
    std::getline(file, line); // Skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        DataEntry entry;
        char comma;
        ss >> entry.b >> comma >> entry.r >> comma >> entry.rotation_1 >> comma
        >> entry.rotation_2 >> comma >> entry.rotation_3 >> comma >> entry.datum_x >> comma >> entry.datum_y >> comma
        >> entry.datum_z >> comma >> entry.is_on_surface >> comma >> entry.geometry;
        
        datas.push_back(entry);
        totalEntries++;
    }
    for (auto data : datas) {

        if (data.geometry < 0) {
            re_scale = 0;
            polyhedron = -data.geometry;
        } else {
            re_scale = 1;
            polyhedron = data.geometry;
        }
        std::cout << "configuration :\n";

        std::cout << std::setprecision(20);

        std::cout << "b : " << data.b << std::endl
                << "r : " << data.r << std::endl
                << "datum_x : " << data.datum_x << std::endl
                << "datum_y : " << data.datum_y << std::endl
                << "datum_z : " << data.datum_z << std::endl
                << "rot_x : " << data.rotation_1 << std::endl
                << "rot_y : " << data.rotation_2 << std::endl
                << "rot_z : " << data.rotation_3 << std::endl
                << "is_on_surface : " << data.is_on_surface << std::endl
                << "is rescaled : " << re_scale << std::endl;

        AlignedCylinder aligned_cylinder = AlignedCylinder({data.b, data.r});

        ReferenceFrame frame(Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0),
                            Normal(0.0, 0.0, 1.0));
        ReferenceFrame frame_amr(Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0),
                            Normal(0.0, 0.0, 1.0));

        double angles[3] = {data.rotation_1, data.rotation_2, 0.0};
        Pt datum(data.datum_x, data.datum_y, data.datum_z);
        Pt cyl_datum(0.0, 0.0, 0.0);

        UnitQuaternion x_rotation(angles[0], frame_amr[0]);
        UnitQuaternion y_rotation(angles[1], frame_amr[1]);
        UnitQuaternion z_rotation(angles[2], frame_amr[2]);
        frame_amr = x_rotation * y_rotation * z_rotation * frame_amr;

        if (old_polyhedron) {
            RectangularCuboid cube = 
                RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0),
                                                    Pt(1.0, 1.0, 1.0));
            RectangularCuboid amr_cube = 
                RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0),
                                                    Pt(1.0, 1.0, 1.0));
            for (auto& vertex : cube) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame_amr[d] * tmp_pt;
                }
            }
            for (auto& vertex : amr_cube) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame_amr[d] * tmp_pt;
                }
            }
            if (data.is_on_surface)
            {
              Pt clip_translation = Pt(0.0, 0.0, sqrt(data.r))-cube[0];
              for (auto& vertex : cube) {
                Pt tmp_pt = vertex + clip_translation;
                vertex = tmp_pt;
              }
              for (auto& vertex : amr_cube) {
                Pt tmp_pt = vertex + clip_translation;
                vertex = tmp_pt;
              }
            }

            HalfEdgePolyhedronQuadratic<Pt> half_edge;
            amr_cube.setHalfEdgeVersion(&half_edge);
            auto seg_half_edge = half_edge.generateSegmentedPolyhedron();
    
    
            Cylinder cylinder(cyl_datum, frame, aligned_cylinder.b(), aligned_cylinder.r());
    
            procede_case(cylinder, cube, seg_half_edge, half_edge, max_amr_level, tri_level);
        } else {
            auto geometri_and_connectivity = getGeometry(polyhedron, re_scale);
            auto geometri = geometri_and_connectivity.first;
            auto connectivity = geometri_and_connectivity.second;
            auto geometri_and_connectivity_amr = getGeometry(polyhedron, re_scale);
            auto geometri_amr = geometri_and_connectivity_amr.first;
            auto connectivity_amr = geometri_and_connectivity_amr.second;
            for (auto& vertex : geometri) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame_amr[d] * tmp_pt;
                }
            }
            for (auto& vertex : geometri_amr) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame_amr[d] * tmp_pt;
                }
            }
            if (data.is_on_surface)
            {
              Pt clip_translation = Pt(0.0, 0.0, sqrt(data.r))-geometri[0];
              for (auto& vertex : geometri) {
                Pt tmp_pt = vertex + clip_translation;
                vertex = tmp_pt;
              }
              for (auto& vertex : geometri_amr) {
                Pt tmp_pt = vertex + clip_translation;
                vertex = tmp_pt;
              }
            }
            
            HalfEdgePolyhedronQuadratic<Pt> half_edge;
            geometri_amr.setHalfEdgeVersion(&half_edge);
            auto seg_half_edge = half_edge.generateSegmentedPolyhedron();
    
    
            Cylinder cylinder(cyl_datum, frame, aligned_cylinder.b(), aligned_cylinder.r());
    
            procede_case(cylinder, geometri, seg_half_edge, half_edge, max_amr_level, tri_level);
        }
      

    }
}

int main(int argc, char* argv[]) {
    int amr_level = 10;
    char* filepath = nullptr;
    int long_id = 0;
    first_vertex_on_surface = 0;
    bool random_coords = false;
    int seed = 0;
    re_scale = 1;
    float tri_level = 0.1;

    while(1) {
        int result = getopt_long(argc, argv, "l:L:f:r::H", long_options, &long_id);
        if (result == -1)
        {
            break;
        }
        switch (result)
        {
        case 'l':
            amr_level = std::atoi(optarg);
            break;
        case 'L':
            tri_level = std::atof(optarg);
            break;
        case 'f':
            filepath = optarg;
            break;
        case 'r':
            random_coords = true;
            if (optarg == nullptr) {
                seed = std::chrono::system_clock::now().time_since_epoch().count();
            } else {
                seed = std::atoi(optarg);
            }
            break;

        case 'H':
            print_help(argv[0]);
            return 0;

        case '?':
            std::cerr << "Invalid argument.\n\n";
            print_help(argv[0]);
            return 1;
        
        case 0:
            break;
        
        default:
            std::cout << "error while reading arguments : " << result << std::endl;
            return 1;
        }
    }

    if (filepath == nullptr && (optind != argc - 1 && optind != argc - 8) && !random_coords) {
        std::cerr << "Usage: " << argv[0] << " -- b r datum_x datum_y datum_z rotation_x rotation_y rotation_z\n";
        std::cerr << "or   : " << argv[0] << " -- b,r,rotation_x,rotation_y,rotation_z,datum_x,datum_y,datum_z[,first_vertex_on_surface[,geometry]]\n";
        std::cerr << "or   : " << argv[0] << " -file [file path]\n";
        std::cerr << "or   : " << argv[0] << " -r <seed>\n";
        return 1;
    }

    std::cout << "amr level : " << amr_level<< std::endl;

    if (filepath == nullptr) {
        double b, r;
        double rotation_x, rotation_y, rotation_z;
        double datum_x, datum_y, datum_z;

        if (random_coords) {
            std::cout << "using random parameters initialized with seed : " << seed << std::endl;
            std::random_device rd;
            std::mt19937_64 eng(seed);

            std::uniform_real_distribution<double> random_rotation(-0.5 * M_PI,
                                                                  0.5 * M_PI);
            std::uniform_real_distribution<double> random_coeffs_b(-5.0, 5.0);
            std::uniform_real_distribution<double> random_coeffs_r(0, 2);
            std::uniform_real_distribution<double> random_translation(-0.5, 0.5);
            b = random_coeffs_b(eng);
            r = random_coeffs_r(eng);
            rotation_x = random_rotation(eng);
            rotation_y = random_rotation(eng);
            rotation_z = random_rotation(eng);
            datum_x = random_translation(eng);
            datum_y = random_translation(eng);
            datum_z = random_translation(eng);

            std::cout << "config : " ;
            std::cout << std::setprecision(20);
            std::cout << b << "," << r << ",";
            std::cout << rotation_x << "," << rotation_y << ","<< rotation_z << ",";
            std::cout << datum_x << "," << datum_y << "," << datum_z << ",";
            std::cout << first_vertex_on_surface << std::endl;
        } else {
            if (optind == argc - 8) {
                sscanf(argv[optind],"%lf",&b);
                sscanf(argv[optind+1],"%lf",&r);

                sscanf(argv[optind+2],"%lf",&datum_x);
                sscanf(argv[optind+3],"%lf",&datum_y);
                sscanf(argv[optind+4],"%lf",&datum_z);

                sscanf(argv[optind+5],"%lf",&rotation_x);
                sscanf(argv[optind+6],"%lf",&rotation_y);
                sscanf(argv[optind+7],"%lf",&rotation_z);
            }

            if (optind == argc - 1) {
                std::stringstream ss(argv[optind]);
                char comma;
                double dbl_fvos;
                double geom;
                ss >> b >> comma >> r >> comma >> rotation_x >> comma
                >> rotation_y >> comma
                >> rotation_z >> comma >> datum_x >> comma >> datum_y >> comma
                >> datum_z >> comma >> dbl_fvos >> comma >> geom;
                first_vertex_on_surface = dbl_fvos;
                if (geom < 0) {
                    re_scale = 0;
                    polyhedron = -geom;
                } else {
                    polyhedron = geom;
                }
            }
        }
        std::cout << "configuration :\n";

        std::cout << std::setprecision(20);

        std::cout << "b : " << b << std::endl
                << "r : " << r << std::endl
                << "datum_x : " << datum_x << std::endl
                << "datum_y : " << datum_y << std::endl
                << "datum_z : " << datum_z << std::endl
                << "rot_x : " << rotation_x << std::endl
                << "rot_y : " << rotation_y << std::endl
                << "rot_z : " << rotation_z << std::endl
                << "is the first vertex on surface ? : " << first_vertex_on_surface << std::endl
                << "geometry use :" << polyhedron << std::endl
                << "rescale ? " << re_scale << std::endl;

        AlignedCylinder aligned_cylinder = AlignedCylinder({b, r});

        ReferenceFrame frame(Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0),
                            Normal(0.0, 0.0, 1.0));
        ReferenceFrame frame_amr(Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0),
                            Normal(0.0, 0.0, 1.0));

        double angles[3] = {rotation_x, rotation_y, rotation_z};
        Pt datum(datum_x, datum_y, datum_z);
        Pt cyl_datum(0.0, 0.0, 0.0);

        UnitQuaternion x_rotation(angles[0], frame_amr[0]);
        UnitQuaternion y_rotation(angles[1], frame_amr[1]);
        UnitQuaternion z_rotation(angles[2], frame_amr[2]);
        frame_amr = x_rotation * y_rotation * z_rotation * frame_amr;
        if (old_polyhedron) {
            RectangularCuboid cube = 
                RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0),
                                                    Pt(1.0, 1.0, 1.0));
            RectangularCuboid amr_cube = 
                RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0),
                                                    Pt(1.0, 1.0, 1.0));
            for (auto& vertex : cube) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame_amr[d] * tmp_pt;
                }
            }
            for (auto& vertex : amr_cube) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame_amr[d] * tmp_pt;
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
              for (auto& vertex : amr_cube) {
                Pt tmp_pt = vertex + clip_translation;
                for (int d = 0; d < 3; d++) {
                    vertex[d] = tmp_pt[d];
                }
              }
            }
    
            HalfEdgePolyhedronQuadratic<Pt> half_edge;
            amr_cube.setHalfEdgeVersion(&half_edge);
            auto seg_half_edge = half_edge.generateSegmentedPolyhedron();
    
    
            Cylinder cylinder(cyl_datum, frame, aligned_cylinder.b(), aligned_cylinder.r());
    
            procede_case(cylinder, cube, seg_half_edge, half_edge, amr_level, tri_level);
        } else {
            auto geometri_and_connectivity = getGeometry(polyhedron, re_scale);
            auto geometri = geometri_and_connectivity.first;
            auto connectivity = geometri_and_connectivity.second;
            auto geometri_and_connectivity_amr = getGeometry(polyhedron, re_scale);
            auto geometri_amr = geometri_and_connectivity_amr.first;
            auto connectivity_amr = geometri_and_connectivity_amr.second;
            for (auto& vertex : geometri) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame_amr[d] * tmp_pt;
                }
            }
            for (auto& vertex : geometri_amr) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame_amr[d] * tmp_pt;
                }
            }
          
            if (first_vertex_on_surface)
            {
              Pt clip_translation = Pt(0.0, 0.0, sqrt(r))-geometri[0];
              for (auto& vertex : geometri) {
                Pt tmp_pt = vertex + clip_translation;
                for (int d = 0; d < 3; d++) {
                    vertex[d] = tmp_pt[d];
                }
              }
              for (auto& vertex : geometri_amr) {
                Pt tmp_pt = vertex + clip_translation;
                for (int d = 0; d < 3; d++) {
                    vertex[d] = tmp_pt[d];
                }
              }
            }
            #ifdef DEBUG_CYL_IRL
            std::cout << "this is the points befor scaling" << std::endl;
            for (auto& vertex : geometri) {
                std::cout << vertex << std::endl;
            }
            #endif
    
            HalfEdgePolyhedronQuadratic<Pt> half_edge;
            geometri_amr.setHalfEdgeVersion(&half_edge);
            auto seg_half_edge = half_edge.generateSegmentedPolyhedron();
    
    
            Cylinder cylinder(cyl_datum, frame, aligned_cylinder.b(), aligned_cylinder.r());
    
            procede_case(cylinder, geometri, seg_half_edge, half_edge, amr_level, tri_level);

        }

    }
    else {
        read_file(filepath, amr_level, tri_level);
    }

    return 0;
}