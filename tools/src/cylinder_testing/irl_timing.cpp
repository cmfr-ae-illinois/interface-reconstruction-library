
#include "irl/geometry/general/new_pt_calculation_functors.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/rotations.h"

#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>
#include <getopt.h>

#include "irl/data_structures/small_vector.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"

#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron_quadratic.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron_quadratic.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/moments/general_moments.h"

#include "irl/cylinder_reconstruction/cylinder.h"
#include "tools/src/cylinder_testing/geometry.h"

using namespace IRL;

int polyhedron;
int old_polyhedron;
int re_scale;

static struct option long_options[] =
    {
        {"cube",                    no_argument, &polyhedron, 0},
        {"old_cube",                no_argument, &old_polyhedron, 1},
        {"tet",                     no_argument, &polyhedron, 1},
        {"dodecahedron",            no_argument, &polyhedron, 2},
        {"cube_hole",               no_argument, &polyhedron, 3},
        {"cube_hole_convex",        no_argument, &polyhedron, 4},
        {"bunny",                   no_argument, &polyhedron, 5},
        {"no_rescale",              no_argument, &re_scale, 0},
        {"force_rescale",           no_argument, &re_scale, 1},
        {"input_file",              required_argument, 0,   'f'},
        {"help",                    no_argument,       0,   'H'},
        {0, 0, 0, 0}
    };

void print_help(std::string prg_name) {
    std::cout << "Usage: " << prg_name << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -H, --help                             Show this help message\n";
    std::cout << "  -f, --input_file <string>              Path to input file (default: result.csv)\n";
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

void displayProgress(int current, int total) {
    int barWidth = 50;
    float progress = (float)current / total;
    int pos = barWidth * progress;
    
    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "]" << current << "/" << total << std::endl;
}

void processEntries(const std::string& filename) {
    std::ifstream datafilecount(filename);
    if (!datafilecount) {
        std::cerr << "Error opening " << filename << std::endl;
        return ;
    }
    datafilecount.unsetf(std::ios_base::skipws);
    std::cout << "counting number of cases" << std::endl;
    std::string line;
    int total = std::count(
        std::istream_iterator<char>(datafilecount),
        std::istream_iterator<char>(), 
        '\n') - 1;
    std::cout << "there is " << total << " cases" << std::endl;
    datafilecount.close();
    std::ifstream datafile(filename);
    datafile.unsetf(std::ios_base::skipws);
    int wrong_result = 0;
    double min_volume_error = __DBL_MAX__;
    double max_volume_error = __DBL_MIN__;
    double min_centroid_x_error = __DBL_MAX__;
    double max_centroid_x_error = __DBL_MIN__;
    double min_centroid_y_error = __DBL_MAX__;
    double max_centroid_y_error = __DBL_MIN__;
    double min_centroid_z_error = __DBL_MAX__;
    double max_centroid_z_error = __DBL_MIN__;
    double avg_volume_error = 0.0;
    double avg_centroid_x_error = 0.0;
    double avg_centroid_y_error = 0.0;
    double avg_centroid_z_error = 0.0;
    double threshold = 1e-13;
    int last_update = 0;
    int update_step = total/50;

    displayProgress(0, total);

    std::string filename1 = "current_configuration.txt";
    std::ofstream file(filename1);
    std::string filename2 = "errors.txt";
    std::ofstream file2(filename2);
    std::string time_data = "time_data.csv";
    std::ofstream time_file(time_data);

    if (!file){
        std::cerr << "Error opening " << filename1 << std::endl;
        return;
    }
    file << std::setprecision(20);
    if (!file2){
        std::cerr << "Error opening " << filename2 << std::endl;
        return;
    }
    file2 << std::setprecision(20);
    if (!time_file){
        std::cerr << "Error opening " << time_data << std::endl;
        return;
    }
    time_file << "time_Volume,time_Moment" << std::endl;
    time_file << std::setprecision(20);
    std::chrono::duration<double, std::micro> time_Volume(0);
    std::chrono::duration<double, std::micro> time_Moment(0);
    int i = 0;

    std::getline(datafile, line);
    while (std::getline(datafile, line)) {
        i++;
        // Example operation: compute sum of all values in the row
        std::stringstream ss(line);
        char comma;
        double b, r, rotation_1, rotation_2, rotation_3, datum_x, datum_y, datum_z, is_on_surface, VolumeData, Centroid_x, Centroid_y, Centroid_z;
        
        ss >> b >> comma >> r >> comma >> rotation_1 >> comma
           >> rotation_2 >> comma >> rotation_3 >> comma >> datum_x >> comma >> datum_y >> comma
           >> datum_z >> comma >> is_on_surface >> comma >> VolumeData >> comma >> Centroid_x >> comma
           >> Centroid_y >> comma >> Centroid_z;

        file << b << "," << r << "," << rotation_1 << "," << rotation_2 << "," << rotation_3 << "," << datum_x << "," << datum_y << "," << datum_z << "," << is_on_surface << "," << VolumeData << "," << Centroid_x << "," << Centroid_y << "," << Centroid_z << std::endl;

        IRL::Pt datum(datum_x, datum_y, datum_z);
        IRL::ReferenceFrame frame(IRL::Normal(1.0, 0.0, 0.0), IRL::Normal(0.0, 1.0, 0.0),
                            IRL::Normal(0.0, 0.0, 1.0));

        double angles[3] = {rotation_1, rotation_2, rotation_3};
        IRL::UnitQuaternion x_rotation(angles[0], frame[0]);
        IRL::UnitQuaternion y_rotation(angles[1], frame[1]);
        IRL::UnitQuaternion z_rotation(angles[2], frame[2]);
        frame = x_rotation * y_rotation * z_rotation * frame;

        Volume result_moments;
        VolumeMoments result_momentsS;

        if (old_polyhedron) {
            IRL::RectangularCuboid cube = IRL::RectangularCuboid::fromBoundingPts(
                IRL::Pt(0.0, 0.0, 0.0), IRL::Pt(1.0, 1.0, 1.0));
    
            for (auto& vertex : cube) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame[d] * tmp_pt;
                }
            }
          
            if (is_on_surface)
            {
              Pt clip_translation = Pt(0.0, 0.0, sqrt(r))-cube[0];
              for (auto& vertex : cube) {
                Pt tmp_pt = vertex + clip_translation;
                vertex = tmp_pt;
              }
            }
    
            IRL::Pt datum0(0.0, 0.0, 0.0);
            IRL::ReferenceFrame frame0(IRL::Normal(1.0, 0.0, 0.0), IRL::Normal(0.0, 1.0, 0.0),
                                IRL::Normal(0.0, 0.0, 1.0));
    
            IRL::Cylinder cylinder(datum0, frame0, b, r);
            asm volatile ("" ::: "memory");
            auto startV = std::chrono::steady_clock::now();
            asm volatile ("" ::: "memory");
            result_moments =
                getVolumeMoments<Volume>(cube, cylinder);
            asm volatile ("" ::: "memory");
            auto endV = std::chrono::steady_clock::now();
            asm volatile ("" ::: "memory");
            time_Volume += std::chrono::duration<double, std::micro>(endV - startV);

            asm volatile ("" ::: "memory");
            auto startM = std::chrono::steady_clock::now();
            asm volatile ("" ::: "memory");
            result_momentsS =
                getVolumeMoments<VolumeMoments>(cube, cylinder);
                asm volatile ("" ::: "memory");
            auto endM = std::chrono::steady_clock::now();
            asm volatile ("" ::: "memory");

            time_Moment += std::chrono::duration<double, std::micro>(endM - startM);

            time_file << std::chrono::duration<double, std::micro>(endV - startV).count() << 
                  "," << std::chrono::duration<double, std::micro>(endM - startM).count() << std::endl;

        } else {

            auto geometri_and_connectivityV = getGeometry(polyhedron, re_scale);
            auto geometriV = geometri_and_connectivityV.first;
            auto connectivityV = geometri_and_connectivityV.second;
            auto geometri_and_connectivityM = getGeometry(polyhedron, re_scale);
            auto geometriM = geometri_and_connectivityM.first;
            auto connectivityM = geometri_and_connectivityM.second;
            for (auto& vertex : geometriV) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame[d] * tmp_pt;
                }
            }
            for (auto& vertex : geometriM) {
                Pt tmp_pt = vertex - datum;
                for (UnsignedIndex_t d = 0; d < 3; ++d) {
                    vertex[d] = frame[d] * tmp_pt;
                }
            }
          
            if (is_on_surface)
            {
              Pt clip_translationV = Pt(0.0, 0.0, sqrt(r))-geometriV[0];
              for (auto& vertex : geometriV) {
                Pt tmp_pt = vertex + clip_translationV;
                for (int d = 0; d < 3; d++) {
                    vertex[d] = tmp_pt[d];
                }
              }
              Pt clip_translationM = Pt(0.0, 0.0, sqrt(r))-geometriM[0];
              for (auto& vertex : geometriM) {
                Pt tmp_pt = vertex + clip_translationM;
                for (int d = 0; d < 3; d++) {
                    vertex[d] = tmp_pt[d];
                }
              }
            }
    
            IRL::Pt datum0(0.0, 0.0, 0.0);
            IRL::ReferenceFrame frame0(IRL::Normal(1.0, 0.0, 0.0), IRL::Normal(0.0, 1.0, 0.0),
                                IRL::Normal(0.0, 0.0, 1.0));
    
            IRL::Cylinder cylinder(datum0, frame0, b, r);
            asm volatile ("" ::: "memory");
            auto startV = std::chrono::steady_clock::now();
            asm volatile ("" ::: "memory");
            result_moments =
                getVolumeMoments<Volume>(geometriV, cylinder);
            asm volatile ("" ::: "memory");
            auto endV = std::chrono::steady_clock::now();
            asm volatile ("" ::: "memory");
            time_Volume += std::chrono::duration<double, std::micro>(endV - startV);

            asm volatile ("" ::: "memory");
            auto startM = std::chrono::steady_clock::now();
            asm volatile ("" ::: "memory");
            result_momentsS =
                getVolumeMoments<VolumeMoments>(geometriM, cylinder);
                asm volatile ("" ::: "memory");
            auto endM = std::chrono::steady_clock::now();
            asm volatile ("" ::: "memory");

            time_Moment += std::chrono::duration<double, std::micro>(endM - startM);

            time_file << std::chrono::duration<double, std::micro>(endV - startV).count() << 
                  "," << std::chrono::duration<double, std::micro>(endM - startM).count() << std::endl;
        }


        auto volume = result_moments;
        auto volumeS = result_momentsS.volume();
        auto volume_error = abs(volume-VolumeData);
        auto volume_errorS = abs(volumeS-VolumeData);

        min_volume_error = minimum(min_volume_error, volume_error);
        min_volume_error = minimum(min_volume_error, volume_errorS);
        max_volume_error = maximum(max_volume_error, volume_error);
        max_volume_error = maximum(max_volume_error, volume_errorS);

        avg_volume_error += volume_error / (2.0 * total);
        avg_volume_error += volume_errorS / (2.0 * total);

        // auto Centroid = result_moments.centroid().getPt();
        auto CentroidS = result_momentsS.centroid().getPt();

        // auto Centroid_x_error = abs(Centroid[0] - Centroid_x);
        // auto Centroid_y_error = abs(Centroid[1] - Centroid_y);
        // auto Centroid_z_error = abs(Centroid[2] - Centroid_z);
        auto Centroid_x_errorS = abs(CentroidS[0] - Centroid_x);
        auto Centroid_y_errorS = abs(CentroidS[1] - Centroid_y);
        auto Centroid_z_errorS = abs(CentroidS[2] - Centroid_z);
        
        // min_centroid_x_error = minimum(min_centroid_x_error, Centroid_x_error);
        // max_centroid_x_error = maximum(max_centroid_x_error, Centroid_x_error);
        // min_centroid_y_error = minimum(min_centroid_y_error, Centroid_y_error);
        // max_centroid_y_error = maximum(max_centroid_y_error, Centroid_y_error);
        // min_centroid_z_error = minimum(min_centroid_z_error, Centroid_z_error);
        // max_centroid_z_error = maximum(max_centroid_z_error, Centroid_z_error);
        min_centroid_x_error = minimum(min_centroid_x_error, Centroid_x_errorS);
        max_centroid_x_error = maximum(max_centroid_x_error, Centroid_x_errorS);
        min_centroid_y_error = minimum(min_centroid_y_error, Centroid_y_errorS);
        max_centroid_y_error = maximum(max_centroid_y_error, Centroid_y_errorS);
        min_centroid_z_error = minimum(min_centroid_z_error, Centroid_z_errorS);
        max_centroid_z_error = maximum(max_centroid_z_error, Centroid_z_errorS);
        
        // avg_centroid_x_error += Centroid_x_error / (2.0 * total);
        // avg_centroid_y_error += Centroid_y_error / (2.0 * total);
        // avg_centroid_z_error += Centroid_z_error / (2.0 * total);
        avg_centroid_x_error += Centroid_x_errorS / (2.0 * total);
        avg_centroid_y_error += Centroid_y_errorS / (2.0 * total);
        avg_centroid_z_error += Centroid_z_errorS / (2.0 * total);

        if (volume_error > threshold || //  || Centroid_x_error > threshold || Centroid_y_error > threshold || Centroid_z_error > threshold
            volume_errorS > threshold || Centroid_x_errorS > threshold || Centroid_y_errorS > threshold || Centroid_z_errorS > threshold) {

            file2 << "The moment were wrong for this configuration :\n" << b << "," << r << "," << rotation_1 << "," << rotation_2 << "," << rotation_3 << "," << datum_x << "," << datum_y << "," << datum_z << "," << is_on_surface << "\n";
            file2 << 
            "actual   results :\n" << volume << "\n" << // "," << Centroid[0] << "," << Centroid[1] << "," << Centroid[2] << "\n" <<
            "actual   results :\n" << volumeS << "," << CentroidS[0] << "," << CentroidS[1] << "," << CentroidS[2] << "\n"
            << "expected results :\n" << VolumeData << "," << Centroid_x << "," << Centroid_y << "," << Centroid_z << "\n"
            << "errors :\n" << volume_error // << "," << Centroid_x_error << "," << Centroid_y_error << "," << Centroid_z_error << "\n"
            << "errors :\n" << volume_errorS << "," << Centroid_x_errorS << "," << Centroid_y_errorS << "," << Centroid_z_errorS << "\n" <<
            "\n\n" << std::endl;
            wrong_result++;
        }

        if (i >= last_update + update_step) {
            // Display progress bar
            displayProgress(i, total);
            last_update = i;
        }
    }
    displayProgress(total, total);
    std::cout << std::endl;

    datafile.close();

    file.close();
    file2.close();
    time_file.close();

    std::cout << std::endl; // Move to new line after completion
    
    std::cout << "done processing all the cases.\nnumber of wrong result : " << wrong_result << std::endl;
    std::cout << std::setprecision(20);
    std::cout << "average error:\n";
    std::cout << "    Volume     : " << avg_volume_error << "\n";
    std::cout << "    Centroid x : " << avg_centroid_x_error << "\n";
    std::cout << "    Centroid y : " << avg_centroid_y_error << "\n";
    std::cout << "    Centroid z : " << avg_centroid_z_error << "\n";

    std::cout << "min:\n";
    std::cout << "    Volume     : " << min_volume_error << "\n";
    std::cout << "    Centroid x : " << min_centroid_x_error << "\n";
    std::cout << "    Centroid y : " << min_centroid_y_error << "\n";
    std::cout << "    Centroid z : " << min_centroid_z_error << "\n";

    std::cout << "max:\n";
    std::cout << "    Volume     : " << max_volume_error << "\n";
    std::cout << "    Centroid x : " << max_centroid_x_error << "\n";
    std::cout << "    Centroid y : " << max_centroid_y_error << "\n";
    std::cout << "    Centroid z : " << max_centroid_z_error << "\n";

    std::cout << "\n\ntime taken to do " << (total) << " cases : \n"
        << "    Only Volume :" << time_Volume.count() / 1e6 << " s or an average of " << time_Volume.count() / (total) << " µs per case\n"
        << "    with Moment :" << time_Moment.count() / 1e6 << " s or an average of " << time_Moment.count() / (total) << " µs per case" << std::endl;
    
    std::cout << "sanity check : i = " << i <<", total = " << total << std::endl;
}

int main(int argc, char* argv[]) {
    std::string file_path = "result.csv";
    int long_id = 0;
    re_scale = 1;
    polyhedron = 0;
    old_polyhedron = 0;
    
    while(1) {
        int result = getopt_long(argc, argv, "f:H", long_options, &long_id);
        if (result == -1)
        {
            break;
        }
        switch (result)
        {
        case 'f':
            file_path = optarg;
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
    
    processEntries(file_path);
    
    return 0;
}
