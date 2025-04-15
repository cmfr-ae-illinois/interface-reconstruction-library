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

#include "irl/data_structures/small_vector.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"

#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron_quadratic.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron_quadratic.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"

#ifndef GETGEOMETRY
#define GETGEOMETRY


using namespace IRL;

std::pair<IRL::StoredGeneralPolyhedron<IRL::Pt, 20>, IRL::PolyhedronConnectivity*> getGeometry(
    const int a_geometry, bool re_scale) {
    std::vector<IRL::Pt> vertices;
    IRL::PolyhedronConnectivity* connectivity = nullptr;
    if (a_geometry == 1) {
        if (re_scale){
            vertices = std::vector<IRL::Pt>{
                {IRL::Pt(1.0 / std::sqrt(3.0), 0.0, -1.0 / std::sqrt(6.0)),
                    IRL::Pt(-std::sqrt(3.0) / 6.0, 0.5, -1.0 / std::sqrt(6.0)),
                    IRL::Pt(-std::sqrt(3.0) / 6.0, -0.5, -1.0 / std::sqrt(6.0)),
                    IRL::Pt(0.0, 0.0, 1.0 / std::sqrt(6.0))}};
            std::array<std::array<IRL::UnsignedIndex_t, 3>, 4> face_mapping{
                {{0, 1, 3}, {1, 2, 3}, {1, 0, 2}, {0, 3, 2}}};
            connectivity = new IRL::PolyhedronConnectivity(face_mapping);
        } else {
            vertices = std::vector<IRL::Pt>{
                {IRL::Pt(0.5, 0.0, -0.5),
                    IRL::Pt(-0.5, 0.5, -0.5),
                    IRL::Pt(-0.5, -0.5, -0.5),
                    IRL::Pt(0.0, 0.0, 0.5)}};
            std::array<std::array<IRL::UnsignedIndex_t, 3>, 4> face_mapping{
                {{0, 1, 3}, {1, 2, 3}, {1, 0, 2}, {0, 3, 2}}};
            connectivity = new IRL::PolyhedronConnectivity(face_mapping);
        }
    } else if (a_geometry == 0) {
        vertices = std::vector<IRL::Pt>{
            {IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(-0.5, 0.5, -0.5),
                IRL::Pt(-0.5, 0.5, 0.5), IRL::Pt(-0.5, -0.5, 0.5),
                IRL::Pt(0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, -0.5),
                IRL::Pt(0.5, 0.5, 0.5), IRL::Pt(0.5, -0.5, 0.5)}};
        std::array<std::array<IRL::UnsignedIndex_t, 4>, 6> face_mapping{
            {{3, 2, 1, 0},
                {0, 1, 5, 4},
                {4, 5, 6, 7},
                {2, 3, 7, 6},
                {1, 2, 6, 5},
                {3, 0, 4, 7}}};
        connectivity = new IRL::PolyhedronConnectivity(face_mapping);
    } else if (a_geometry == 2) {
        const double tau = (sqrt(5.0) + 1.0) / 2.0;
        std::array<IRL::Pt, 12> M{{IRL::Pt(0.0, tau, 1.0), IRL::Pt(0.0, -tau, 1.0),
                                    IRL::Pt(0.0, tau, -1.0),
                                    IRL::Pt(0.0, -tau, -1.0), IRL::Pt(1.0, 0.0, tau),
                                    IRL::Pt(-1.0, 0.0, tau), IRL::Pt(1.0, 0.0, -tau),
                                    IRL::Pt(-1.0, 0.0, -tau), IRL::Pt(tau, 1.0, 0.0),
                                    IRL::Pt(-tau, 1.0, 0.0), IRL::Pt(tau, -1.0, 0.0),
                                    IRL::Pt(-tau, -1.0, 0.0)}};
        vertices = std::vector<IRL::Pt>{{(1.0 / 3.0) * (M[0] + M[8] + M[2]),
                                            (1.0 / 3.0) * (M[0] + M[4] + M[8]),
                                            (1.0 / 3.0) * (M[0] + M[5] + M[4]),
                                            (1.0 / 3.0) * (M[0] + M[9] + M[5]),
                                            (1.0 / 3.0) * (M[0] + M[2] + M[9]),
                                            (1.0 / 3.0) * (M[2] + M[8] + M[6]),
                                            (1.0 / 3.0) * (M[8] + M[10] + M[6]),
                                            (1.0 / 3.0) * (M[8] + M[4] + M[10]),
                                            (1.0 / 3.0) * (M[4] + M[1] + M[10]),
                                            (1.0 / 3.0) * (M[4] + M[5] + M[1]),
                                            (1.0 / 3.0) * (M[5] + M[11] + M[1]),
                                            (1.0 / 3.0) * (M[5] + M[9] + M[11]),
                                            (1.0 / 3.0) * (M[9] + M[7] + M[11]),
                                            (1.0 / 3.0) * (M[9] + M[2] + M[7]),
                                            (1.0 / 3.0) * (M[2] + M[6] + M[7]),
                                            (1.0 / 3.0) * (M[3] + M[10] + M[1]),
                                            (1.0 / 3.0) * (M[3] + M[1] + M[11]),
                                            (1.0 / 3.0) * (M[3] + M[11] + M[7]),
                                            (1.0 / 3.0) * (M[3] + M[7] + M[6]),
                                            (1.0 / 3.0) * (M[3] + M[6] + M[10])}};
        const double scale = 0.5;
        for (auto& pt : vertices) {
            pt *= scale;
        }
        std::array<std::array<IRL::UnsignedIndex_t, 5>, 12> face_mapping{
            {{5, 4, 3, 2, 1},
                {1, 2, 8, 7, 6},
                {2, 3, 10, 9, 8},
                {3, 4, 12, 11, 10},
                {4, 5, 14, 13, 12},
                {5, 1, 6, 15, 14},
                {16, 17, 18, 19, 20},
                {16, 20, 7, 8, 9},
                {9, 10, 11, 17, 16},
                {11, 12, 13, 18, 17},
                {6, 7, 20, 19, 15},
                {13, 14, 15, 19, 18}}};
        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 5; ++j) {
            --face_mapping[i][j];
            }
        }

        connectivity = new IRL::PolyhedronConnectivity(face_mapping);
    } else if (a_geometry == 3) {
        vertices = std::vector<IRL::Pt>{
            {IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(-0.5, 0.5, -0.5),
                IRL::Pt(-0.5, 0.5, 0.5), IRL::Pt(-0.5, -0.5, 0.5),
                IRL::Pt(0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, -0.5),
                IRL::Pt(0.5, 0.5, 0.5), IRL::Pt(0.5, -0.5, 0.5),
                IRL::Pt(-0.25, 0.5, -0.25), IRL::Pt(-0.25, 0.5, 0.25),
                IRL::Pt(0.25, 0.5, 0.25), IRL::Pt(0.25, 0.5, -0.25),
                IRL::Pt(-0.25, -0.5, -0.25), IRL::Pt(-0.25, -0.5, 0.25),
                IRL::Pt(0.25, -0.5, 0.25), IRL::Pt(0.25, -0.5, -0.25)}};

        std::vector<std::vector<IRL::UnsignedIndex_t>> face_mapping(12);
        face_mapping[0] = std::vector<IRL::UnsignedIndex_t>{{3, 2, 1, 0}};
        face_mapping[1] = std::vector<IRL::UnsignedIndex_t>{{0, 1, 5, 4}};
        face_mapping[2] = std::vector<IRL::UnsignedIndex_t>{{4, 5, 6, 7}};
        face_mapping[3] = std::vector<IRL::UnsignedIndex_t>{{6, 2, 3, 7}};
        face_mapping[4] = std::vector<IRL::UnsignedIndex_t>{{1, 2, 6, 10, 9, 8}};
        face_mapping[5] = std::vector<IRL::UnsignedIndex_t>{{6, 5, 1, 8, 11, 10}};
        face_mapping[6] = std::vector<IRL::UnsignedIndex_t>{{0, 4, 7, 14, 15, 12}};
        face_mapping[7] = std::vector<IRL::UnsignedIndex_t>{{7, 3, 0, 12, 13, 14}};
        face_mapping[8] = std::vector<IRL::UnsignedIndex_t>{{8, 9, 13, 12}};
        face_mapping[9] = std::vector<IRL::UnsignedIndex_t>{{9, 10, 14, 13}};
        face_mapping[10] = std::vector<IRL::UnsignedIndex_t>{{10, 11, 15, 14}};
        face_mapping[11] = std::vector<IRL::UnsignedIndex_t>{{11, 8, 12, 15}};

        connectivity = new IRL::PolyhedronConnectivity(face_mapping);

    } else if (a_geometry == 4) {
        vertices = std::vector<IRL::Pt>{
            {IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(-0.5, 0.5, -0.5),
                IRL::Pt(-0.5, 0.5, 0.5), IRL::Pt(-0.5, -0.5, 0.5),
                IRL::Pt(0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, -0.5),
                IRL::Pt(0.5, 0.5, 0.5), IRL::Pt(0.5, -0.5, 0.5),
                IRL::Pt(-0.25, 0.5, -0.25), IRL::Pt(-0.25, 0.5, 0.25),
                IRL::Pt(0.25, 0.5, 0.25), IRL::Pt(0.25, 0.5, -0.25),
                IRL::Pt(-0.25, -0.5, -0.25), IRL::Pt(-0.25, -0.5, 0.25),
                IRL::Pt(0.25, -0.5, 0.25), IRL::Pt(0.25, -0.5, -0.25)}};

        std::vector<std::vector<IRL::UnsignedIndex_t>> face_mapping(16);
        face_mapping[0] = std::vector<IRL::UnsignedIndex_t>{{3, 2, 1, 0}};
        face_mapping[1] = std::vector<IRL::UnsignedIndex_t>{{0, 1, 5, 4}};
        face_mapping[2] = std::vector<IRL::UnsignedIndex_t>{{4, 5, 6, 7}};
        face_mapping[3] = std::vector<IRL::UnsignedIndex_t>{{6, 2, 3, 7}};
        face_mapping[4] = std::vector<IRL::UnsignedIndex_t>{{1, 2, 9, 8}};
        face_mapping[5] = std::vector<IRL::UnsignedIndex_t>{{2, 6, 10, 9}};
        face_mapping[6] = std::vector<IRL::UnsignedIndex_t>{{6, 5, 11, 10}};
        face_mapping[7] = std::vector<IRL::UnsignedIndex_t>{{5, 1, 8, 11}};
        face_mapping[8] = std::vector<IRL::UnsignedIndex_t>{{0, 4, 15, 12}};
        face_mapping[9] = std::vector<IRL::UnsignedIndex_t>{{4, 7, 14, 15}};
        face_mapping[10] = std::vector<IRL::UnsignedIndex_t>{{7, 3, 13, 14}};
        face_mapping[11] = std::vector<IRL::UnsignedIndex_t>{{3, 0, 12, 13}};
        face_mapping[12] = std::vector<IRL::UnsignedIndex_t>{{8, 9, 13, 12}};
        face_mapping[13] = std::vector<IRL::UnsignedIndex_t>{{9, 10, 14, 13}};
        face_mapping[14] = std::vector<IRL::UnsignedIndex_t>{{10, 11, 15, 14}};
        face_mapping[15] = std::vector<IRL::UnsignedIndex_t>{{11, 8, 12, 15}};

        connectivity = new IRL::PolyhedronConnectivity(face_mapping);
    } else if (a_geometry == 5) {
        std::ifstream myfile("../vtk_data/bunny.vtk");
        std::string line;
        char space_char = ' ';
        std::vector<std::array<IRL::UnsignedIndex_t, 3>> face_mapping{};
        if (myfile.is_open()) {
            std::vector<std::string> words{};
            bool read_points = false;
            bool read_connectivity = false;
            while (getline(myfile, line)) {
            std::stringstream line_stream(line);
            std::string split_line;
            while (std::getline(line_stream, split_line, space_char)) {
                words.push_back(split_line);
            }

            if (read_points) {
                for (unsigned i = 0; i < words.size() / 3; i++) {
                vertices.push_back(IRL::Pt(std::stod(words[3 * i + 0]),
                                            std::stod(words[3 * i + 1]),
                                            std::stod(words[3 * i + 2])));
                }
            }
            if (read_connectivity) {
                for (unsigned i = 0; i < words.size() / 3; i++) {
                face_mapping.push_back(std::array<IRL::UnsignedIndex_t, 3>{
                    static_cast<IRL::UnsignedIndex_t>(std::stoi(words[3 * i + 0])),
                    static_cast<IRL::UnsignedIndex_t>(std::stoi(words[3 * i + 1])),
                    static_cast<IRL::UnsignedIndex_t>(
                        std::stoi(words[3 * i + 2]))});
                }
            }
            if (words.size() == 0 || words[0] == "METADATA") {
                read_points = false;
                read_connectivity = false;
            } else if (words[0] == "POINTS") {
                read_points = true;
                read_connectivity = false;
            } else if (words[0] == "CONNECTIVITY") {
                read_points = false;
                read_connectivity = true;
            } else if (words[0] == "CELL_DATA") {
                read_points = false;
                read_connectivity = false;
            }
            words.clear();
            }
            myfile.close();
        }
        std::cout << "The bunny has " << vertices.size() << " vertices!"
                    << std::endl;
        std::cout << "The bunny has " << face_mapping.size() << " triangles!"
                    << std::endl;
        connectivity = new IRL::PolyhedronConnectivity(face_mapping);
    } else {
        std::cout << "Unkown geometry type \"" << a_geometry << "\"" << std::endl;
        std::exit(-1);
    }

    auto polyhedron = IRL::StoredGeneralPolyhedron<IRL::Pt, 20>(vertices, connectivity);

    if (re_scale) {
        // Normalize volume
        const double normalization_factor =
            1.0 / std::pow(polyhedron.calculateVolume(), 1.0 / 3.0);
        for (auto& vertex : polyhedron) {
            for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
                vertex[d] *= normalization_factor;
            }
        }
        // Place centroid at origin
        const auto centroid = polyhedron.calculateCentroid();
        for (auto& vertex : polyhedron) {
            for (IRL::UnsignedIndex_t d = 0; d < 3; ++d) {
                vertex[d] -= centroid[d];
            }
        }
    }

    return std::make_pair(polyhedron, connectivity);
}

#endif