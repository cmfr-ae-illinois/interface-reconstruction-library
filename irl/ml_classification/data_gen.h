#pragma once

#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/ml_classification/vtk_out.h"
#include "irl/ml_classification/inertia_calc.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

namespace IRL {
    class Data_gen {

        public:
        std::mt19937_64 eng;
        double machineZero = 10e-14;
        //int stencil_size = 3;
        //std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));


        Data_gen(){
            // Initialize random number generator with a seed from current time
            std::random_device rd;
            eng = std::mt19937_64(rd());  // Use random_device for seeding
            std::cout << "I'm a data generator!" << std::endl;
        }

        Eigen::Vector3d computeCentroidFromFirstMoment(const Eigen::Vector3d& firstMoment, double volume){
            if (volume == 0.0) {
                return Eigen::Vector3d::Zero();  // safe fallback
            }
            return firstMoment / volume;
        }

        // Print a 3D field of Eigen::Vector3d centroids, used if visualize=true for data generation
        void printCentroids(const std::vector<std::vector<std::vector<Eigen::Vector3d>>>& centroid)
        {
            std::cout << "Centroid:" << std::endl;

            for (int i = 0; i < centroid.size(); ++i) {
                for (int j = 0; j < centroid[i].size(); ++j) {
                    for (int k = 0; k < centroid[i][j].size(); ++k) {
                        const auto& c = centroid[i][j][k];
                        std::cout << "[" << i << "," << j << "," << k << "] "
                                << std::fixed << std::setprecision(3)
                                << c.x() << ", " << c.y() << ", " << c.z() << "   ";
                    }
                    std::cout << std::endl;
                }
            }
        }

        Eigen::Vector3d generateRandomDirection(std::mt19937_64& eng) {
            std::uniform_real_distribution<double> dist01(0.0, 1.0);
            double u = dist01(eng);              // uniform in [0,1]
            double v = dist01(eng);              // uniform in [0,1]

            double theta = 2.0 * M_PI * u;       // azimuth angle uniform in [0, 2π)
            double z = 2.0 * v - 1.0;            // uniform in [-1,1]
            double r = std::sqrt(1.0 - z*z);

            return Eigen::Vector3d(r * std::cos(theta), 
                                r * std::sin(theta), 
                                z);
        }

        std::vector<double> generateParaboloid(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double coeff_stddev = 0.1, bool visualize = false) 
        {
            // repeat until center cell is cut by surface
            while (true) {
                // make centroid, only used for visualization
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                // Defining cell coordinates
                auto coords = std::vector<double>(stencil_size+1);
                for (int i = 0; i < stencil_size+1; i++) {
                    coords[i] = -0.5*static_cast<double>(stencil_size) + static_cast<double>(i);
                }
                const double cell_volume = 1.0;

                // Random paraboloid parameters
                Eigen::Vector3d datumVec = generateRandomPoint(
                    -0.5*static_cast<double>(stencil_size),
                    0.5*static_cast<double>(stencil_size), eng);

                Eigen::Vector3d direction = generateRandomDirection(eng);

                // Build orthonormal frame aligned with direction
                Eigen::Vector3d helper = generateRandomDirection(eng); // safe random helper
                Eigen::Vector3d paraboloid_x = direction.cross(helper);
                if (paraboloid_x.norm() < 1e-12) {
                    // fallback: pick fixed axis if helper was parallel
                    paraboloid_x = direction.cross(Eigen::Vector3d(1,0,0));
                }
                paraboloid_x.normalize();
                Eigen::Vector3d paraboloid_y = direction.cross(paraboloid_x);
                paraboloid_y.normalize();
                Eigen::Vector3d paraboloid_z = direction;

                const auto frame = IRL::ReferenceFrame(
                    IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                    IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                    IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));

                IRL::Pt datum(datumVec.x(), datumVec.y(), datumVec.z());

                std::normal_distribution<double> random_coeff(0.0, coeff_stddev);
                double coeff1 = random_coeff(eng);
                double coeff2 = random_coeff(eng);

                const auto paraboloid = Paraboloid(datum, frame, coeff1, coeff2);

                // Initialize field
                std::vector<ParaboloidParametrizedSurfaceOutput> surfaces;
                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;
                for (int i = 0; i < stencil_size; i++) {
                    for (int j = 0; j < stencil_size; j++) {
                        for (int k = 0; k < stencil_size; k++) {
                            auto cell = RectangularCuboid::fromBoundingPts(
                                Pt(coords[i], coords[j], coords[k]),
                                Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                            // auto volume_and_surface = getVolumeMoments<
                             //   AddSurfaceOutput<VolumeMoments, ParametrizedSurfaceOutput>>(
                              //  cell, paraboloid);

                            auto volume_and_surface = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid);

                            auto surface = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid).getSurface();

                            //surfaces.push_back(volume_and_surface.getSurface());
                            surfaces.push_back(surface);
                            vfrac[i][j][k] = volume_and_surface.getMoments().volume() / cell_volume;
                            firstMoment[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                    volume_and_surface.getMoments().centroid().y(),
                                                    volume_and_surface.getMoments().centroid().z();

                            if (visualize) {
                                centroid[i][j][k] = computeCentroidFromFirstMoment(
                                    firstMoment[i][j][k], vfrac[i][j][k] * cell_volume);
                            }
                        }
                    }
                }

                // --- check central cell ---
                int mid = stencil_size / 2;
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0-machineZero) {
                    // Accept this sample
                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    // Flatten the 3D vector vfrac into a 1D vector
                    std::vector<double> flattened_vfrac;
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {
                                flattened_vfrac.push_back(vfrac[i][j][k]);
                            }
                        }
                    }
                    return flattened_vfrac;
                }

                // else: reject and try again (loop restarts)
            }
        }

        double sample_truncated_normal(double mean, double stddev,
                               double lower, double upper)
        {
            std::normal_distribution<double> normal(mean, stddev);
            double x;
            do {
                x = normal(eng);
            } while (x < lower || x > upper);
            return x;
        }


        std::vector<double> generateSheet(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double coeff_stddev = 0.1, double max_thickness = 0.5, double thickness_stddev = 0.0, bool visualize = false) 
        {
            while (true) { // keep trying until center cell has surface crossing
                // make centroid, only used for visualization
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                // Defining cell coordinates
                auto coords = std::vector<double>(stencil_size+1);
                for (int i = 0; i < stencil_size+1; i++) {
                    coords[i] = -0.5*static_cast<double>(stencil_size) + static_cast<double>(i);
                }
                const double cell_volume = 1.0;

                // Random datum anywhere in stencil
                Eigen::Vector3d datum = generateRandomPoint(
                    -0.5*static_cast<double>(stencil_size),
                    0.5*static_cast<double>(stencil_size), eng);

                // Random unbiased direction
                Eigen::Vector3d direction = generateRandomDirection(eng);

                // Random sheet thickness
                double thickness = max_thickness;
                if (thickness_stddev > 0.0) {
                    thickness = sample_truncated_normal(0, thickness_stddev, machineZero, max_thickness);
                }else{
                    std::uniform_real_distribution<double> random_thickness(machineZero, max_thickness);
                    thickness = random_thickness(eng);
                }

                // Two paraboloid datums offset along the direction
                Eigen::Vector3d datum_paraboloid1_eVec = datum - direction.normalized() * (thickness/2.0);
                Eigen::Vector3d datum_paraboloid2_eVec = datum + direction.normalized() * (thickness/2.0);

                // Build orthonormal frame aligned with direction
                Eigen::Vector3d helper = generateRandomDirection(eng);
                Eigen::Vector3d paraboloid_x = direction.cross(helper);
                if (paraboloid_x.norm() < 1e-12) {
                    paraboloid_x = direction.cross(Eigen::Vector3d(1,0,0));
                }
                paraboloid_x.normalize();
                Eigen::Vector3d paraboloid_y = direction.cross(paraboloid_x);
                paraboloid_y.normalize();
                Eigen::Vector3d paraboloid_z = direction;

                const auto frame = IRL::ReferenceFrame(
                    IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                    IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                    IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));

                IRL::Pt datum_paraboloid1(datum_paraboloid1_eVec.x(), datum_paraboloid1_eVec.y(), datum_paraboloid1_eVec.z());
                IRL::Pt datum_paraboloid2(datum_paraboloid2_eVec.x(), datum_paraboloid2_eVec.y(), datum_paraboloid2_eVec.z());

                // Random coefficients
                std::normal_distribution<double> random_coeff(0.0, coeff_stddev);
                double coeff1 = random_coeff(eng);
                double coeff2 = random_coeff(eng);

                const auto paraboloid1 = Paraboloid(datum_paraboloid1, frame, coeff1, coeff2);
                const auto paraboloid2 = Paraboloid(datum_paraboloid2, frame, coeff1, coeff2);

                // Initialize field
                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;
                std::vector<ParaboloidParametrizedSurfaceOutput> surfaces;
                for (int i = 0; i < stencil_size; i++) {
                    for (int j = 0; j < stencil_size; j++) {
                        for (int k = 0; k < stencil_size; k++) {
                            auto cell = RectangularCuboid::fromBoundingPts(
                                Pt(coords[i], coords[j], coords[k]),
                                Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                            /*
                            auto volume_and_surface1 = getVolumeMoments<
                                AddSurfaceOutput<VolumeMoments, ParametrizedSurfaceOutput>>(
                                cell, paraboloid1);
                            auto volume_and_surface2 = getVolumeMoments<
                                AddSurfaceOutput<VolumeMoments, ParametrizedSurfaceOutput>>(
                                cell, paraboloid2);
                            */

                            auto volume_and_surface1 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid1);

                            auto surface1 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid1).getSurface();

                            auto volume_and_surface2 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid2);

                            auto surface2 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid2).getSurface();

                            //surfaces.push_back(volume_and_surface1.getSurface());
                            //surfaces.push_back(volume_and_surface2.getSurface());
                            surfaces.push_back(surface1);
                            surfaces.push_back(surface2);

                            double V1 = volume_and_surface1.getMoments().volume();
                            double V2 = volume_and_surface2.getMoments().volume();
                            double Vdiff = V2 - V1;

                            if (Vdiff < 0.0) Vdiff = 0.0;

                            vfrac[i][j][k] = Vdiff / cell_volume;

                            Eigen::Vector3d M1(volume_and_surface1.getMoments().centroid().x(),
                                            volume_and_surface1.getMoments().centroid().y(),
                                            volume_and_surface1.getMoments().centroid().z());

                            Eigen::Vector3d M2(volume_and_surface2.getMoments().centroid().x(),
                                            volume_and_surface2.getMoments().centroid().y(),
                                            volume_and_surface2.getMoments().centroid().z());

                            firstMoment[i][j][k] = M2 - M1;

                            if (visualize) {
                                centroid[i][j][k] = computeCentroidFromFirstMoment(
                                    firstMoment[i][j][k], vfrac[i][j][k] * cell_volume);
                            }
                        }
                    }
                }

                // --- check central cell ---
                int mid = stencil_size / 2;
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0-machineZero) {
                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    std::vector<double> flattened_vfrac;
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {
                                flattened_vfrac.push_back(vfrac[i][j][k]);
                            }
                        }
                    }
                    return flattened_vfrac; // accept this sample
                }

                // else: reject and regenerate
            }
        }

        /*
        std::vector<double> generateCutSheet(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double coeff_stddev = 0.1, double max_thickness = 0.5, double thickness_stddev = 0.0, bool visualize = false,
            int refinement_factr = 3,
            double Rmin = 5.0, double Rmax = 100.0, double Rmean = 5.0, double Rstd = 3.0)
        {
            while (true) { // keep trying until center cell has surface crossing
                // make centroid, only used for visualization
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                // Defining cell coordinates
                auto coords = std::vector<double>(stencil_size+1);
                for (int i = 0; i < stencil_size+1; i++) {
                    coords[i] = -0.5*static_cast<double>(stencil_size) + static_cast<double>(i);
                }
                const double cell_volume = 1.0;

                // Random datum anywhere in stencil
                Eigen::Vector3d datum = generateRandomPoint(
                    -0.5*static_cast<double>(stencil_size),
                    0.5*static_cast<double>(stencil_size), eng);

                // Random unbiased direction
                Eigen::Vector3d direction = generateRandomDirection(eng);

                // Random sheet thickness
                double thickness = max_thickness;
                if (thickness_stddev > 0.0) {
                    thickness = sample_truncated_normal(0, thickness_stddev, machineZero, max_thickness);
                }else{
                    std::uniform_real_distribution<double> random_thickness(machineZero, max_thickness);
                    thickness = random_thickness(eng);
                }

                // Two paraboloid datums offset along the direction
                Eigen::Vector3d datum_paraboloid1_eVec = datum - direction.normalized() * (thickness/2.0);
                Eigen::Vector3d datum_paraboloid2_eVec = datum + direction.normalized() * (thickness/2.0);

                // Build orthonormal frame aligned with direction
                Eigen::Vector3d helper = generateRandomDirection(eng);
                Eigen::Vector3d paraboloid_x = direction.cross(helper);
                if (paraboloid_x.norm() < 1e-12) {
                    paraboloid_x = direction.cross(Eigen::Vector3d(1,0,0));
                }
                paraboloid_x.normalize();
                Eigen::Vector3d paraboloid_y = direction.cross(paraboloid_x);
                paraboloid_y.normalize();
                Eigen::Vector3d paraboloid_z = direction;

                const auto frame = IRL::ReferenceFrame(
                    IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                    IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                    IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));

                IRL::Pt datum_paraboloid1(datum_paraboloid1_eVec.x(), datum_paraboloid1_eVec.y(), datum_paraboloid1_eVec.z());
                IRL::Pt datum_paraboloid2(datum_paraboloid2_eVec.x(), datum_paraboloid2_eVec.y(), datum_paraboloid2_eVec.z());

                // Random coefficients
                std::normal_distribution<double> random_coeff(0.0, coeff_stddev);
                double coeff1 = random_coeff(eng);
                double coeff2 = random_coeff(eng);

                const auto paraboloid1 = Paraboloid(datum_paraboloid1, frame, coeff1, coeff2);
                const auto paraboloid2 = Paraboloid(datum_paraboloid2, frame, coeff1, coeff2);
                
                // Sheet cutting
                std::uniform_real_distribution<double> uni01(0.0, 1.0);
                bool inward = (uni01(eng) < 0.5);
                double cuttingRadius = sample_truncated_normal(Rmean, Rstd, Rmin, Rmax);
                
                std::uniform_real_distribution<double> uniCu(cu_min, cu_max);
                std::uniform_real_distribution<double> uniCv(cv_min, cv_max);
                double cuttingRadiusCenterU = uniCu(eng);
                double cuttingRadiusCenterV = uniCv(eng);

                // Initialize field
                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;
                std::vector<ParaboloidParametrizedSurfaceOutput> surfaces;
                for (int i = 0; i < stencil_size; i++) {
                    for (int j = 0; j < stencil_size; j++) {
                        for (int k = 0; k < stencil_size; k++) {
                            auto cell = RectangularCuboid::fromBoundingPts(
                                Pt(coords[i], coords[j], coords[k]),
                                Pt(coords[i + 1], coords[j + 1], coords[k + 1]));



                            auto volume_and_surface1 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid1);

                            auto surface1 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid1).getSurface();

                            auto volume_and_surface2 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid2);

                            auto surface2 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid2).getSurface();

                            //surfaces.push_back(volume_and_surface1.getSurface());
                            //surfaces.push_back(volume_and_surface2.getSurface());
                            surfaces.push_back(surface1);
                            surfaces.push_back(surface2);

                            double V1 = volume_and_surface1.getMoments().volume();
                            double V2 = volume_and_surface2.getMoments().volume();
                            double Vdiff = V2 - V1;

                            if (Vdiff < 0.0) Vdiff = 0.0;

                            // Apply cutting radius -> hard mask


                            vfrac[i][j][k] = Vdiff / cell_volume;

                            Eigen::Vector3d M1(volume_and_surface1.getMoments().centroid().x(),
                                            volume_and_surface1.getMoments().centroid().y(),
                                            volume_and_surface1.getMoments().centroid().z());

                            Eigen::Vector3d M2(volume_and_surface2.getMoments().centroid().x(),
                                            volume_and_surface2.getMoments().centroid().y(),
                                            volume_and_surface2.getMoments().centroid().z());

                            firstMoment[i][j][k] = M2 - M1;

                            if (visualize) {
                                centroid[i][j][k] = computeCentroidFromFirstMoment(
                                    firstMoment[i][j][k], vfrac[i][j][k] * cell_volume);
                            }
                        }
                    }
                }

                // --- check central cell ---
                int mid = stencil_size / 2;
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0-machineZero) {
                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    std::vector<double> flattened_vfrac;
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {
                                flattened_vfrac.push_back(vfrac[i][j][k]);
                            }
                        }
                    }
                    return flattened_vfrac; // accept this sample
                }

                // else: reject and regenerate
            }
        }
        */

        // Function to generate a random vector within given bounds
        Eigen::Vector3d generateRandomPoint(double min_bound, double max_bound, std::mt19937_64& eng) {
            std::uniform_real_distribution<double> dist(min_bound, max_bound);
            return Eigen::Vector3d(dist(eng), dist(eng), dist(eng));  // Random point in 3D space
        }

        void compressStencilRefinedToCoarse(
            const std::vector<std::vector<std::vector<double>>>& volumes_refined,
            const std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoments_refined,
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size,
            int refinement_factor,
            double compressed_cell_volume,
            bool visualize,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>* centroid = nullptr)
        {
            double refinement_factor_double = static_cast<double>(refinement_factor);

            for (int i = 0; i < stencil_size; i++) {
                for (int j = 0; j < stencil_size; j++) {
                    for (int k = 0; k < stencil_size; k++) {
                        double vol_sum = 0.0;
                        Eigen::Vector3d firstMoment_sum = Eigen::Vector3d::Zero();

                        for (int m = i * refinement_factor; m < (i + 1) * refinement_factor; m++) {
                            for (int n = j * refinement_factor; n < (j + 1) * refinement_factor; n++) {
                                for (int o = k * refinement_factor; o < (k + 1) * refinement_factor; o++) {
                                    double v = volumes_refined[m][n][o];
                                    vol_sum += v;
                                    firstMoment_sum += firstMoments_refined[m][n][o];
                                }
                            }
                        }


                        if (vol_sum > machineZero) {
                            firstMoment[i][j][k] = firstMoment_sum;
                            vfrac[i][j][k] = vol_sum / compressed_cell_volume;
                        } else {
                            firstMoment[i][j][k] = Eigen::Vector3d::Zero();
                            vfrac[i][j][k] = 0.0;
                        }

                        if (visualize && centroid) {
                            (*centroid)[i][j][k] = computeCentroidFromFirstMoment(
                                firstMoment[i][j][k],
                                vol_sum
                            );
                        }
                    }
                }
            }
        }

        /*
        std::vector<double> generateCutSheet(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double coeff_stddev = 0.1, double max_thickness = 0.5, double thickness_stddev = 0.0, bool visualize = false,
            double Rmin = 5.0, double Rmax = 100.0, double Rmean = 5.0, double Rstd = 3.0)
        {
            while (true) { // keep trying until center cell has surface crossing
                // make centroid, only used for visualization
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                // for visualization option
                std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;

                const double cell_volume = 1.0;
                int refinement_factor = 3;
                double refinement_factor_double = static_cast<double>(refinement_factor);

                // Random datum anywhere in stencil
                Eigen::Vector3d datum = generateRandomPoint(
                    -0.5*static_cast<double>(stencil_size),
                    0.5*static_cast<double>(stencil_size), eng);

                // Random unbiased direction
                Eigen::Vector3d direction = generateRandomDirection(eng);

                // Random sheet thickness
                double thickness = max_thickness;
                if (thickness_stddev > 0.0) {
                    thickness = sample_truncated_normal(0, thickness_stddev, machineZero, max_thickness);
                }else{
                    std::uniform_real_distribution<double> random_thickness(machineZero, max_thickness);
                    thickness = random_thickness(eng);
                }

                // Two paraboloid datums offset along the direction
                Eigen::Vector3d datum_paraboloid1_eVec = datum - direction.normalized() * (thickness/2.0);
                Eigen::Vector3d datum_paraboloid2_eVec = datum + direction.normalized() * (thickness/2.0);

                // Refined mesh
                int refined_stencil_size = refinement_factor*stencil_size;
                std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                    std::vector<std::vector<double>>(refined_stencil_size,
                    std::vector<double>(refined_stencil_size)));

                std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                    std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

                // Refined cell coordinates
                auto coords = std::vector<double>(refined_stencil_size+1);
                for (int i = 0; i <= refined_stencil_size; i++) {
                    coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                }

                                // Build orthonormal frame aligned with direction
                Eigen::Vector3d helper = generateRandomDirection(eng);
                Eigen::Vector3d paraboloid_x = direction.cross(helper);
                if (paraboloid_x.norm() < 1e-12) {
                    paraboloid_x = direction.cross(Eigen::Vector3d(1,0,0));
                }
                paraboloid_x.normalize();
                Eigen::Vector3d paraboloid_y = direction.cross(paraboloid_x);
                paraboloid_y.normalize();
                Eigen::Vector3d paraboloid_z = direction;

                const auto frame = IRL::ReferenceFrame(
                    IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                    IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                    IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));

                IRL::Pt datum_paraboloid1(datum_paraboloid1_eVec.x(), datum_paraboloid1_eVec.y(), datum_paraboloid1_eVec.z());
                IRL::Pt datum_paraboloid2(datum_paraboloid2_eVec.x(), datum_paraboloid2_eVec.y(), datum_paraboloid2_eVec.z());

                // Random coefficients
                std::normal_distribution<double> random_coeff(0.0, coeff_stddev);
                double coeff1 = random_coeff(eng);
                double coeff2 = random_coeff(eng);

                const auto paraboloid1 = Paraboloid(datum_paraboloid1, frame, coeff1, coeff2);
                const auto paraboloid2 = Paraboloid(datum_paraboloid2, frame, coeff1, coeff2);

                // Initialize field
                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;
                std::vector<ParaboloidParametrizedSurfaceOutput> surfaces;

                // Fill refined stencil
                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;
                for (int i = 0; i < refined_stencil_size; i++) {
                    for (int j = 0; j < refined_stencil_size; j++) {
                        for (int k = 0; k < refined_stencil_size; k++) {
                            auto cell = IRL::RectangularCuboid::fromBoundingPts(
                                IRL::Pt(coords[i], coords[j], coords[k]),
                                IRL::Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                            // Find center of cell
                            Eigen::Vector3d cell_center((coords[i+1]+coords[i])/2.0,
                                                        (coords[j+1]+coords[j])/2.0,
                                                        (coords[k+1]+coords[k])/2.0);

                            // Projection of cell center onto cylinder axis
                            Eigen::Vector3d axis_to_cell = cell_center - axis_origin;
                            double projection_factor = axis_to_cell.dot(axis_direction) / axis_direction.squaredNorm();
                            Eigen::Vector3d closest_point_on_axis = axis_origin + projection_factor * axis_direction;

                            // Datum on paraboloid representation of cylinder
                            Eigen::Vector3d datum_paraboloid_eVec =
                                closest_point_on_axis + radius * (cell_center - closest_point_on_axis).normalized();
                            IRL::Pt datum_paraboloid(datum_paraboloid_eVec.x(), datum_paraboloid_eVec.y(), datum_paraboloid_eVec.z());

                            // Frame aligned with axis
                            Eigen::Vector3d paraboloid_x = axis_direction.normalized();
                            Eigen::Vector3d paraboloid_z = (cell_center - closest_point_on_axis).normalized();
                            Eigen::Vector3d paraboloid_y = paraboloid_z.cross(paraboloid_x);
                            paraboloid_y.normalize();

                            const auto frame = IRL::ReferenceFrame(
                                IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                                IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                                IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));

                            // Cylinder = paraboloid with coeffs (0, 1/(2R))
                            const auto paraboloid = IRL::Paraboloid(datum_paraboloid, frame, 0, 1/(2*radius));

                            auto volume_and_surface = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid);

                            volumes_refined[i][j][k] = volume_and_surface.getMoments().volume();
                            firstMoments_refined[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                            volume_and_surface.getMoments().centroid().y(),
                                                            volume_and_surface.getMoments().centroid().z();

                            if (visualize) {
                                auto surface = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid).getSurface();
                                surfaces.push_back(surface);
                            }
                        }
                    }
                }

                // Compress refined → coarse stencil
                compressStencilRefinedToCoarse(
                    volumes_refined,
                    firstMoments_refined,
                    vfrac,
                    firstMoment,
                    stencil_size,
                    refinement_factor,
                    cell_volume,
                    visualize,
                    &centroid
                );

                // --- check central cell ---
                int mid = stencil_size / 2;
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0-machineZero) {
                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    std::vector<double> flattened_vfrac;
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {
                                flattened_vfrac.push_back(vfrac[i][j][k]);
                            }
                        }
                    }
                    return flattened_vfrac; // accept this sample
                }
                // else: reject and regenerate
            }   
        }
        */


        std::vector<double> generateCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double max_radius = 0.5, double radius_stddev = 0.0, bool visualize = false) 
        {
            while (true) { // keep trying until center cell has surface crossing
                // make centroid, only used for visualization
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                // for visualization option
                std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;

                const double cell_volume = 1.0;
                int refinement_factor = 3;
                double refinement_factor_double = static_cast<double>(refinement_factor);

                Eigen::Vector3d axis_origin = generateRandomPoint(
                    -0.5 * stencil_size, 0.5 * stencil_size, eng);
                Eigen::Vector3d axis_direction = generateRandomDirection(eng);

                // Random radius
                double radius = max_radius;
                if (radius_stddev > 0.0) {
                    radius = sample_truncated_normal(0, radius_stddev, machineZero, max_radius);
                }else{
                    std::uniform_real_distribution<double> random_thickness(machineZero, max_radius);
                    radius = random_thickness(eng);
                }
                
                // Compute distance from stencil center (0,0,0) to cylinder axis
                Eigen::Vector3d center(0.0, 0.0, 0.0);
                Eigen::Vector3d origin_to_center = axis_origin - center;
                double projection = origin_to_center.dot(axis_direction);
                Eigen::Vector3d closest_point = axis_origin - projection * axis_direction;
                double distance_to_axis = (closest_point - center).norm();

                // Quick reject: if axis is too far from center, skip this attempt
                if (std::abs(distance_to_axis - radius) > 0.8661) {
                    continue; // try again
                }
                

                // Refined mesh
                int refined_stencil_size = refinement_factor*stencil_size;
                std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                    std::vector<std::vector<double>>(refined_stencil_size,
                    std::vector<double>(refined_stencil_size)));

                std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                    std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

                // Refined cell coordinates
                auto coords = std::vector<double>(refined_stencil_size+1);
                for (int i = 0; i <= refined_stencil_size; i++) {
                    coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                }

                // Fill refined stencil
                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;
                for (int i = 0; i < refined_stencil_size; i++) {
                    for (int j = 0; j < refined_stencil_size; j++) {
                        for (int k = 0; k < refined_stencil_size; k++) {
                            auto cell = IRL::RectangularCuboid::fromBoundingPts(
                                IRL::Pt(coords[i], coords[j], coords[k]),
                                IRL::Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                            // Find center of cell
                            Eigen::Vector3d cell_center((coords[i+1]+coords[i])/2.0,
                                                        (coords[j+1]+coords[j])/2.0,
                                                        (coords[k+1]+coords[k])/2.0);

                            // Projection of cell center onto cylinder axis
                            Eigen::Vector3d axis_to_cell = cell_center - axis_origin;
                            double projection_factor = axis_to_cell.dot(axis_direction) / axis_direction.squaredNorm();
                            Eigen::Vector3d closest_point_on_axis = axis_origin + projection_factor * axis_direction;

                            // Datum on paraboloid representation of cylinder
                            Eigen::Vector3d datum_paraboloid_eVec =
                                closest_point_on_axis + radius * (cell_center - closest_point_on_axis).normalized();
                            IRL::Pt datum_paraboloid(datum_paraboloid_eVec.x(), datum_paraboloid_eVec.y(), datum_paraboloid_eVec.z());

                            // Frame aligned with axis
                            Eigen::Vector3d paraboloid_x = axis_direction.normalized();
                            Eigen::Vector3d paraboloid_z = (cell_center - closest_point_on_axis).normalized();
                            Eigen::Vector3d paraboloid_y = paraboloid_z.cross(paraboloid_x);
                            paraboloid_y.normalize();

                            const auto frame = IRL::ReferenceFrame(
                                IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                                IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                                IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));

                            // Cylinder = paraboloid with coeffs (0, 1/(2R))
                            const auto paraboloid = IRL::Paraboloid(datum_paraboloid, frame, 0, 1/(2*radius));

                            /*
                            auto volume_and_surface = IRL::getVolumeMoments<
                                IRL::AddSurfaceOutput<IRL::VolumeMoments, IRL::ParametrizedSurfaceOutput>>(
                                cell, paraboloid);
                            */
                            auto volume_and_surface = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid);

                            volumes_refined[i][j][k] = volume_and_surface.getMoments().volume();
                            firstMoments_refined[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                            volume_and_surface.getMoments().centroid().y(),
                                                            volume_and_surface.getMoments().centroid().z();

                            if (visualize) {
                                auto surface = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid).getSurface();
                                surfaces.push_back(surface);
                            }
                        }
                    }
                }

                // Compress refined → coarse stencil
                compressStencilRefinedToCoarse(
                    volumes_refined,
                    firstMoments_refined,
                    vfrac,
                    firstMoment,
                    stencil_size,
                    refinement_factor,
                    cell_volume,
                    visualize,
                    &centroid
                );

                // --- check central cell ---
                int mid = stencil_size / 2;
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0-machineZero) {
                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    std::vector<double> flattened_vfrac;
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {
                                flattened_vfrac.push_back(vfrac[i][j][k]);
                            }
                        }
                    }
                    return flattened_vfrac; // accept this sample
                }
                // else: reject and regenerate
            }   
        }

        std::vector<double> generateTruncatedCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double max_radius = 0.5, double radius_stddev = 0.0, bool visualize = false) 
        {
            while (true) { // keep trying until center cell has surface crossing
                // make centroid, only used for visualization
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                // for visualization option
                std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;

                const double cell_volume = 1.0;
                int refinement_factor = 3;
                double refinement_factor_double = static_cast<double>(refinement_factor);

                Eigen::Vector3d axis_origin = generateRandomPoint(
                    -0.5 * stencil_size, 0.5 * stencil_size, eng);
                Eigen::Vector3d axis_direction = generateRandomDirection(eng);

                // Random radius
                double radius = max_radius;
                if (radius_stddev > 0.0) {
                    radius = sample_truncated_normal(0, radius_stddev, machineZero, max_radius);
                }else{
                    std::uniform_real_distribution<double> random_thickness(machineZero, max_radius);
                    radius = random_thickness(eng);
                }

                
                // Compute distance from stencil center (0,0,0) to cylinder axis
                Eigen::Vector3d center(0.0, 0.0, 0.0);
                Eigen::Vector3d origin_to_center = axis_origin - center;
                double projection = origin_to_center.dot(axis_direction);

                // For truncated case: closest point must be in +direction (projection < 0 means center is behind)
                if (projection > 0.0) {
                    // The axis runs away from center; no chance of intersection
                    continue; // try again
                }

                // Distance from axis to center
                Eigen::Vector3d closest_point = axis_origin + projection * axis_direction;
                double distance_to_axis = (center - closest_point).norm();
                // Intersection feasibility check
                if (std::abs(distance_to_axis - radius) > 0.8661) {
                    continue; // surface cannot cross central cell
                }
                

                // Refined mesh
                int refined_stencil_size = refinement_factor*stencil_size;
                std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                    std::vector<std::vector<double>>(refined_stencil_size,
                    std::vector<double>(refined_stencil_size)));

                std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                    std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

                // Refined cell coordinates
                auto coords = std::vector<double>(refined_stencil_size+1);
                for (int i = 0; i <= refined_stencil_size; i++) {
                    coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                }

                // Fill refined stencil
                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;
                for (int i = 0; i < refined_stencil_size; i++) {
                    for (int j = 0; j < refined_stencil_size; j++) {
                        for (int k = 0; k < refined_stencil_size; k++) {
                            auto cell = IRL::RectangularCuboid::fromBoundingPts(
                                IRL::Pt(coords[i], coords[j], coords[k]),
                                IRL::Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                            // Find center of cell
                            Eigen::Vector3d cell_center((coords[i+1]+coords[i])/2.0,
                                                        (coords[j+1]+coords[j])/2.0,
                                                        (coords[k+1]+coords[k])/2.0);

                            // Projection of cell center onto cylinder axis
                            Eigen::Vector3d axis_to_cell = cell_center - axis_origin;
                            double projection_factor = axis_to_cell.dot(axis_direction) / axis_direction.squaredNorm();

                            // Truncate cylinder at axis origin
                            if (projection_factor < 0.0) {
                                volumes_refined[i][j][k] = 0.0;
                                firstMoments_refined[i][j][k] = Eigen::Vector3d::Zero();
                                continue; // skip this cell
                            }

                            Eigen::Vector3d closest_point_on_axis = axis_origin + projection_factor * axis_direction;

                            // Datum on paraboloid representation of cylinder
                            Eigen::Vector3d datum_paraboloid_eVec =
                                closest_point_on_axis + radius * (cell_center - closest_point_on_axis).normalized();
                            IRL::Pt datum_paraboloid(datum_paraboloid_eVec.x(), datum_paraboloid_eVec.y(), datum_paraboloid_eVec.z());

                            // Frame aligned with axis
                            Eigen::Vector3d paraboloid_x = axis_direction.normalized();
                            Eigen::Vector3d paraboloid_z = (cell_center - closest_point_on_axis).normalized();
                            Eigen::Vector3d paraboloid_y = paraboloid_z.cross(paraboloid_x);
                            paraboloid_y.normalize();

                            const auto frame = IRL::ReferenceFrame(
                                IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                                IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                                IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));

                            // Cylinder = paraboloid with coeffs (0, 1/(2R))
                            const auto paraboloid = IRL::Paraboloid(datum_paraboloid, frame, 0, 1/(2*radius));

                            /*
                            auto volume_and_surface = IRL::getVolumeMoments<
                                IRL::AddSurfaceOutput<IRL::VolumeMoments, IRL::ParametrizedSurfaceOutput>>(
                                cell, paraboloid);
                            */
                            auto volume_and_surface = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid);

                            volumes_refined[i][j][k] = volume_and_surface.getMoments().volume();
                            firstMoments_refined[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                            volume_and_surface.getMoments().centroid().y(),
                                                            volume_and_surface.getMoments().centroid().z();

                            if (visualize) {
                                auto surface = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid).getSurface();
                                surfaces.push_back(surface);
                            }
                        }
                    }
                }

                // Compress refined → coarse stencil
                compressStencilRefinedToCoarse(
                    volumes_refined,
                    firstMoments_refined,
                    vfrac,
                    firstMoment,
                    stencil_size,
                    refinement_factor,
                    cell_volume,
                    visualize,
                    &centroid
                );

                // --- check central cell ---
                int mid = stencil_size / 2;
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0-machineZero) {
                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    std::vector<double> flattened_vfrac;
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {
                                flattened_vfrac.push_back(vfrac[i][j][k]);
                            }
                        }
                    }
                    return flattened_vfrac; // accept this sample
                }
                // else: reject and regenerate
            }
                
        }

        std::vector<double> generateSphere (
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double max_radius = 0.5, double radius_stddev = 0.0, bool visualize = false) 
            {
                while (true) { // keep trying until center cell has surface crossing
                

                    // make centroid, only used for visualization
                    std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                        stencil_size,
                        std::vector<std::vector<Eigen::Vector3d>>(
                            stencil_size,
                            std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                        )
                    );

                    // for visualization option:
                    std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;

                    // define parameters
                    //std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));
                    const double cell_volume = 1.0;
                    int refinement_factor = 3;
                    double refinement_factor_double = static_cast<double>(refinement_factor);

                    // Random radius
                    double radius = max_radius;
                    if (radius_stddev > 0.0) {
                        radius = sample_truncated_normal(0, radius_stddev, machineZero, max_radius);
                    }else{
                        std::uniform_real_distribution<double> random_thickness(machineZero, max_radius);
                        radius = random_thickness(eng);
                    }

                    Eigen::Vector3d origin = generateRandomPoint(-0.5 - radius , 0.5 + radius, eng);

                    // Create refined cells
                    int refined_stencil_size = refinement_factor*stencil_size;
                    std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                        std::vector<std::vector<double>>(refined_stencil_size,
                        std::vector<double>(refined_stencil_size)));

                    std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                        std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                        std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

                    //std::vector<std::vector<std::vector<double>>> volumes(refined_stencil_size, std::vector<std::vector<double>>(refined_stencil_size, std::vector<double>(refined_stencil_size)));

                    // Defining cell coordinates
                    auto coords = std::vector<double>(refined_stencil_size+1);
                    for (int i = 0; i <= refined_stencil_size; i++) {
                        coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                    }

                    using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

                    for (int i = 0; i < refined_stencil_size; i++) {
                        for (int j = 0; j < refined_stencil_size; j++) {
                            for (int k = 0; k < refined_stencil_size; k++) {
                                // Create cell
                                auto cell = IRL::RectangularCuboid::fromBoundingPts(
                                    IRL::Pt(coords[i], coords[j], coords[k]),
                                    IRL::Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                                // Create paraboloid for this cell in the steps below
                                // Find center point of cell
                                Eigen::Vector3d cell_center((coords[i+1]+coords[i])/2.0,
                                            (coords[j+1]+coords[j])/2.0,
                                            (coords[k+1]+coords[k])/2.0);

                                // Find datum of paraboloid
                                Eigen::Vector3d origin_to_cell_center = cell_center - origin;
                                Eigen::Vector3d datum_paraboloid_eVec = origin + radius * (cell_center - origin).normalized();
                                IRL::Pt datum_paraboloid(datum_paraboloid_eVec.x(), datum_paraboloid_eVec.y(), datum_paraboloid_eVec.z());

                                // Create frame
                                Eigen::Vector3d paraboloid_z = cell_center-origin;
                                paraboloid_z.normalize();
                                Eigen::Vector3d helper (paraboloid_z.x(), paraboloid_z.y()+1, paraboloid_z.z());
                                Eigen::Vector3d paraboloid_x = paraboloid_z.cross(helper);
                                paraboloid_x.normalize();
                                Eigen::Vector3d paraboloid_y = paraboloid_z.cross(paraboloid_x);
                                paraboloid_y.normalize();
                                //frame[1] = IRL::crossProduct(frame[2],frame[0]); frame2 is normal, frame0 is axis
                                //frame[1].normalize();

                                const auto frame = IRL::ReferenceFrame(IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                                                    IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                                                    IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));
                                // Create paraboloid
                                const auto paraboloid = IRL::Paraboloid(datum_paraboloid, frame, 1/(2*radius), 1/(2*radius));

                                // Intersect cell with paraboloid -- return volume and surface
                                //auto volume_fraction = IRL::getVolumeFraction(cell, paraboloid);
                                /*
                                auto volume_and_surface = IRL::getVolumeMoments<
                                    IRL::AddSurfaceOutput<IRL::VolumeMoments, IRL::ParametrizedSurfaceOutput>>(
                                    cell, paraboloid);
                                */
                                auto volume_and_surface = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid);

                                
                                //auto allMoments = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid);
                                
                                volumes_refined[i][j][k] = volume_and_surface.getMoments().volume();
                                                //volume_and_surface.getMoments().centroid() / vol

                                Eigen::Vector3d m1(volume_and_surface.getMoments().centroid().x(),
                                                volume_and_surface.getMoments().centroid().y(),
                                                volume_and_surface.getMoments().centroid().z());

                                firstMoments_refined[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                                volume_and_surface.getMoments().centroid().y(),
                                                                volume_and_surface.getMoments().centroid().z();

                                if (visualize) {
                                    auto surface = getVolumeMoments<
                                        VolumeMomentsAndSurface>(cell, paraboloid).getSurface();

                                    surfaces.push_back(surface);
                                }
                            }
                        }
                    }

                    // Compress refined mesh into original stencil size
                    compressStencilRefinedToCoarse(
                        volumes_refined,
                        firstMoments_refined,
                        vfrac,
                        firstMoment,
                        stencil_size,
                        refinement_factor,
                        cell_volume,
                        visualize,
                        &centroid
                    );

                    // --- check central cell ---
                    int mid = stencil_size / 2;
                    double center_vfrac = vfrac[mid][mid][mid];
                    if (center_vfrac > machineZero && center_vfrac < 1.0-machineZero) {
                        if (visualize) {
                            WriteField(stencil_size, coords, vfrac, "vfrac");
                            WriteSurface(surfaces, "surface");
                            printCentroids(centroid);
                        }

                        std::vector<double> flattened_vfrac;
                        for (int i = 0; i < stencil_size; ++i) {
                            for (int j = 0; j < stencil_size; ++j) {
                                for (int k = 0; k < stencil_size; ++k) {
                                    flattened_vfrac.push_back(vfrac[i][j][k]);
                                }
                            }
                        }
                        return flattened_vfrac; // accept this sample
                    }
                    // else: reject and regenerate
                    
                }
        }


        std::vector<float> generateState(int datapoint_type, int stencil_size, int include_Moments = 0, bool include_truncated_cylinder = false,
                                        double paraboloid_coeff_stddev = 0.1,
                                        double sheet_coeff_stddev = 0.1, double max_sheet_thickness = 0.5, double sheet_thickness_stddev = 0.0,
                                        double max_cylinder_radius = 0.5, double cylinder_radius_stddev = 0.0, 
                                        double max_sphere_radius = 0.5, double sphere_radius_stddev = 0.0,
                                        bool visualize = false)
        {
            // Create a vector of volume fractions and first moments
            std::vector<std::vector<std::vector<double>>> vfrac(
                stencil_size, 
                std::vector<std::vector<double>>(
                    stencil_size, 
                    std::vector<double>(
                        stencil_size
                    )
                )
            );

            std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoment(
                stencil_size,
                std::vector<std::vector<Eigen::Vector3d>>(
                    stencil_size,
                    std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                )
            );

            switch (datapoint_type) {
                case 1:
                    if (include_truncated_cylinder) {
                        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
                        double p = prob_dist(eng);  // draw a random number in [0,1)

                        if (p < 0.2) {
                            // 20% chance → truncated cylinder
                            generateTruncatedCylinder(
                                vfrac, firstMoment, stencil_size,
                                max_cylinder_radius, cylinder_radius_stddev, visualize);
                        } else {
                            // 80% chance → normal cylinder
                            generateCylinder(
                                vfrac, firstMoment, stencil_size,
                                max_cylinder_radius, cylinder_radius_stddev, visualize);
                        }
                    } else {
                        generateCylinder(vfrac, firstMoment, stencil_size, max_cylinder_radius, cylinder_radius_stddev, visualize);
                    }
                    break;
                case 2:
                    generateSphere(vfrac, firstMoment, stencil_size, max_sphere_radius, sphere_radius_stddev, visualize);
                    break;
                case 3:
                    generateSheet(vfrac, firstMoment, stencil_size, sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev, visualize);
                    break;
                case 4:
                    std::cout<<"Cut-Sheet generation not implemented yet."<<std::endl;
                    //generate_Cut_Sheet(vfrac, firstMoment, stencil_size, sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev, visualize);
                    break;
                default:
                    generateParaboloid(vfrac, firstMoment, stencil_size, paraboloid_coeff_stddev, visualize);
                    break;
            }

            // Flatten the 3D vectors into one 1D vector
            std::vector<double> flattened_state;
            for (int i = 0; i < stencil_size; ++i) {
                for (int j = 0; j < stencil_size; ++j) {
                    for (int k = 0; k < stencil_size; ++k) {
                        if (include_Moments >= 0) {
                            flattened_state.push_back(vfrac[i][j][k]);
                        }
                        if (include_Moments >= 1) {
                            flattened_state.push_back(firstMoment[i][j][k].x());
                            flattened_state.push_back(firstMoment[i][j][k].y());
                            flattened_state.push_back(firstMoment[i][j][k].z());
                        }
                    }
                }
            }

            if (include_Moments >= 2) {
                Eigen::Matrix3d inertia = IRL::computeInertiaTensor(flattened_state, stencil_size, 1);

                flattened_state.push_back(inertia(0, 0)); // Ixx
                flattened_state.push_back(inertia(1, 1)); // Iyy
                flattened_state.push_back(inertia(2, 2)); // Izz
                flattened_state.push_back(inertia(0, 1)); // Ixy
                flattened_state.push_back(inertia(0, 2)); // Ixz
                flattened_state.push_back(inertia(1, 2)); // Iyz
            }

            if (include_Moments >= 3) {
                // if include_Moments == 3, vector should have only inertia eigenvalues
                Eigen::Matrix3d I = IRL::computeInertiaTensor(flattened_state, stencil_size, 0);

                // Get eigenvalues
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(I);
                Eigen::Vector3d evals = solver.eigenvalues();

                // Sort eigenvalues descending: I1 >= I2 >= I3
                std::sort(evals.data(), evals.data() + 3, std::greater<double>());
                double I1 = evals[0], I2 = evals[1], I3 = evals[2];

                if (include_Moments == 4) {
                    // if include_Moments == 4, use only the three eigenvalues
                    flattened_state.clear();
                }

                flattened_state.push_back(I1);
                flattened_state.push_back(I2);
                flattened_state.push_back(I3);
            }

            // make flattened_state a vector of floats
            std::vector<float> flattened_state_float(flattened_state.begin(), flattened_state.end());
            return flattened_state_float;
        }

        void generateData (std::vector<std::vector<float>>* statesV, std::vector<int>* labelsV, int no_datapoints, int stencil_size = 3, int no_datapoint_types_in = 4, int include_Moments = 0,
                                        double paraboloid_coeff_stddev = 0.1,
                                        double sheet_coeff_stddev = 0.1, double max_sheet_thickness = 0.5, double sheet_thickness_stddev = 0.0,
                                        double max_cylinder_radius = 0.5, double cylinder_radius_stddev = 0.0, bool include_truncated_cylinder = false,
                                        double max_sphere_radius = 0.5, double sphere_radius_stddev = 0.0)
        
        {
            std::cout << no_datapoints << std::endl;
            std::cout << include_truncated_cylinder << std::endl;
            // Initialize random number generator with a seed from current time
            std::srand(std::time(0));

            for (int i=0; i<no_datapoints; i++) {
                
                // Generate the data, init with a random number 0 or 1, 0=paraboloid, 1=cylinder, 2=sphere, 3=sheet

                int datapoint_type = std::rand() % no_datapoint_types_in;
                //std::cout << "Generating datapoint " << i << " of type " << datapoint_type << std::endl;

                labelsV->push_back(datapoint_type);
                statesV->push_back(generateState(datapoint_type, stencil_size, include_Moments, include_truncated_cylinder,
                                        paraboloid_coeff_stddev,
                                        sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev,
                                        max_cylinder_radius, cylinder_radius_stddev, 
                                        max_sphere_radius, sphere_radius_stddev));
                if (i % 1000 == 0) {
                    std::cout << "Generated " << i << " datapoints." << std::endl;
                }
            }
        }
    };
}