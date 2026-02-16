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

        // Moment helpers below

        inline std::vector<double> makeCenteredCoords(int stencil_size) {
            std::vector<double> coords(stencil_size + 1);
            for (int i = 0; i <= stencil_size; ++i) {
                coords[i] = -0.5 * static_cast<double>(stencil_size) + static_cast<double>(i);
            }
            return coords;
        }

        inline std::vector<double> flattenVfrac(const std::vector<std::vector<std::vector<double>>>& vfrac,
                                                int stencil_size) {
            std::vector<double> flattened;
            flattened.reserve(static_cast<size_t>(stencil_size) * stencil_size * stencil_size);
            for (int i = 0; i < stencil_size; ++i)
                for (int j = 0; j < stencil_size; ++j)
                    for (int k = 0; k < stencil_size; ++k)
                        flattened.push_back(vfrac[i][j][k]);
            return flattened;
        }

        inline bool centerCellIsCut(const std::vector<std::vector<std::vector<double>>>& vfrac,
                                    int stencil_size,
                                    double machineZero) {
            int mid = stencil_size / 2;
            const double a = vfrac[mid][mid][mid];
            return (a > machineZero && a < 1.0 - machineZero);
        }

        // Accumulate moments and fill vfrac/firstMoment using GeneralMoments
        template <class CellType, class ReconstructionType>
        inline void fillCellFromGeneralMoments2(const CellType& cell,
                                                const ReconstructionType& recon,
                                                double cell_volume,
                                                double& out_vfrac,
                                                Eigen::Vector3d& out_firstMoment,
                                                IRL::GeneralMoments3D<2>& io_totalMoments) {
            auto gm = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, recon);

            const double vol = gm.volume();
            out_vfrac = vol / cell_volume;

            out_firstMoment << gm[1], gm[2], gm[3];

            io_totalMoments += gm; // accumulate for stencil -> calculate 2nd moment later
        }

        // Fill vfrac/firstMoment using VolumeMoments
        template <class VolumeMomentsAndSurfaceType, class CellType, class ReconstructionType>
        inline void fillCellFromVolumeMomentsAndSurface(const CellType& cell,
                                                        const ReconstructionType& recon,
                                                        double cell_volume,
                                                        double& out_vfrac,
                                                        Eigen::Vector3d& out_firstMoment) {
            auto vas = IRL::getVolumeMoments<VolumeMomentsAndSurfaceType>(cell, recon);
            out_vfrac = vas.getMoments().volume() / cell_volume;
            out_firstMoment << vas.getMoments().centroid().x(),
                            vas.getMoments().centroid().y(),
                            vas.getMoments().centroid().z();
        }

        // Compute stencil-centered (about global liquid centroid) 2nd moment tensor from the accumulated 2nd moments
        inline Eigen::Matrix3d centeredSecondMomentFromTotal(IRL::GeneralMoments3D<2> totalMoments) {
            const double V = totalMoments.volume();

            // global centroid of stencil liquid in the SAME coordinate system as moments were computed
            IRL::Pt C(0.0, 0.0, 0.0);
            if (V > 0.0) {
                C = IRL::Pt(totalMoments[1] / V, totalMoments[2] / V, totalMoments[3] / V);
            }

            // Shift moments so origin becomes centroid: x_new = x_old - C  -> datum = -C
            totalMoments.moveMoments(IRL::Pt(-C[0], -C[1], -C[2]));

            Eigen::Matrix3d M2;
            M2 << totalMoments[4], totalMoments[5], totalMoments[6],
                totalMoments[5], totalMoments[7], totalMoments[8],
                totalMoments[6], totalMoments[8], totalMoments[9];
            return M2;
        }
        
        
        std::vector<double> generateParaboloid(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double coeff_stddev = 0.1, bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr) 
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
                auto coords = makeCenteredCoords(stencil_size);
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
                IRL::GeneralMoments3D<2> totalMoments = IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For 2nd moment shift later

                // Loop over cells in stencil
                for (int i = 0; i < stencil_size; i++) {
                    for (int j = 0; j < stencil_size; j++) {
                        for (int k = 0; k < stencil_size; k++) {
                            auto cell = RectangularCuboid::fromBoundingPts(
                                Pt(coords[i], coords[j], coords[k]),
                                Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                            // auto volume_and_surface = getVolumeMoments<
                             //   AddSurfaceOutput<VolumeMoments, ParametrizedSurfaceOutput>>(
                              //  cell, paraboloid);

                            if (secondMoment != nullptr) {
                                fillCellFromGeneralMoments2(
                                    cell, paraboloid, cell_volume,
                                    vfrac[i][j][k],
                                    firstMoment[i][j][k],
                                    totalMoments
                                );

                            } else {
                                auto volume_and_surface = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid);
                                vfrac[i][j][k] = volume_and_surface.getMoments().volume() / cell_volume;
                                firstMoment[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                        volume_and_surface.getMoments().centroid().y(),
                                                        volume_and_surface.getMoments().centroid().z();
                            }

                            if (visualize) {
                                auto surface = getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid).getSurface();
                                //surfaces.push_back(volume_and_surface.getSurface());
                                surfaces.push_back(surface);

                                centroid[i][j][k] = computeCentroidFromFirstMoment(
                                    firstMoment[i][j][k], vfrac[i][j][k] * cell_volume);
                            }
                        }
                    }
                }

                // Check central cell
                if (centerCellIsCut(vfrac, stencil_size, machineZero)) {
                    // Accept this sample
                    
                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        *secondMoment = centeredSecondMomentFromTotal(totalMoments);
                    }

                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    // Flatten the 3D vector vfrac into a 1D vector
                    return flattenVfrac(vfrac, stencil_size);
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
            int stencil_size, double coeff_stddev = 0.1, double max_thickness = 0.5, double thickness_stddev = 0.0, bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr) 
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
                auto coords = makeCenteredCoords(stencil_size);
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
                IRL::GeneralMoments3D<2> totalMoments = IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For 2nd moment shift later

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
                            if (secondMoment != nullptr) {
                                auto gm1 = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid1);
                                auto gm2 = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid2);
                                auto gmSheet = gm2 - gm1;

                                // 0th moment
                                double vol = gmSheet.volume();

                                vfrac[i][j][k] = vol / cell_volume;

                                // 1st moments
                                firstMoment[i][j][k] << gmSheet[1], gmSheet[2], gmSheet[3];

                                if (vol < 0.0) {
                                    vol = 0.0;
                                    firstMoment[i][j][k] = Eigen::Vector3d::Zero();
                                    //secondMoment[i][j][k] = Eigen::Matrix3d::Zero(); Causes issues
                                }

                                // accumulate for stencil-centered 2nd moment later
                                totalMoments += gmSheet;

                            } else {
                                auto volume_and_surface1 = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid1);

                                auto volume_and_surface2 = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid2);

                                double V1 = volume_and_surface1.getMoments().volume();
                                double V2 = volume_and_surface2.getMoments().volume();
                                double Vdiff = V2 - V1;



                                Eigen::Vector3d M1(volume_and_surface1.getMoments().centroid().x(),
                                                volume_and_surface1.getMoments().centroid().y(),
                                                volume_and_surface1.getMoments().centroid().z());

                                Eigen::Vector3d M2(volume_and_surface2.getMoments().centroid().x(),
                                                volume_and_surface2.getMoments().centroid().y(),
                                                volume_and_surface2.getMoments().centroid().z());


                                if (Vdiff <= 0.0) {
                                    vfrac[i][j][k] = 0.0;
                                    firstMoment[i][j][k] = Eigen::Vector3d::Zero();
                                } else {
                                    vfrac[i][j][k] = Vdiff / cell_volume;
                                    firstMoment[i][j][k] = M2 - M1;
                                }
                            }

                            if (visualize) {
                                auto surface1 = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid1).getSurface();
                                auto surface2 = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid2).getSurface();
                                //surfaces.push_back(volume_and_surface1.getSurface());
                                //surfaces.push_back(volume_and_surface2.getSurface());
                                surfaces.push_back(surface1);
                                surfaces.push_back(surface2);
                                centroid[i][j][k] = computeCentroidFromFirstMoment(
                                    firstMoment[i][j][k], vfrac[i][j][k] * cell_volume);
                            }
                        }
                    }
                }

                // Check central cell
                if (centerCellIsCut(vfrac, stencil_size, machineZero)) {
                    // Accept this sample
                    
                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        *secondMoment = centeredSecondMomentFromTotal(totalMoments);
                    }

                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    // Flatten the 3D vector vfrac into a 1D vector
                    return flattenVfrac(vfrac, stencil_size);
                }

                // else: reject and regenerate
            }
        }

        std::vector<double> generateCutSheet(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double coeff_stddev = 0.1, double max_thickness = 0.5, double thickness_stddev = 0.0, bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr) 
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
                auto coords = makeCenteredCoords(stencil_size);
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
                IRL::GeneralMoments3D<2> totalMoments = IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For 2nd moment shift later

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
                            if (secondMoment != nullptr) {
                                auto gm1 = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid1);
                                auto gm2 = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid2);
                                auto gmSheet = gm2 - gm1;

                                // 0th moment
                                double vol = gmSheet.volume();

                                vfrac[i][j][k] = vol / cell_volume;

                                // 1st moments
                                firstMoment[i][j][k] << gmSheet[1], gmSheet[2], gmSheet[3];

                                if (vol < 0.0) {
                                    vol = 0.0;
                                    firstMoment[i][j][k] = Eigen::Vector3d::Zero();
                                    //secondMoment[i][j][k] = Eigen::Matrix3d::Zero(); Causes issues
                                }

                                // accumulate for stencil-centered 2nd moment later
                                totalMoments += gmSheet;

                            } else {
                                auto volume_and_surface1 = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid1);

                                auto volume_and_surface2 = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid2);

                                double V1 = volume_and_surface1.getMoments().volume();
                                double V2 = volume_and_surface2.getMoments().volume();
                                double Vdiff = V2 - V1;



                                Eigen::Vector3d M1(volume_and_surface1.getMoments().centroid().x(),
                                                volume_and_surface1.getMoments().centroid().y(),
                                                volume_and_surface1.getMoments().centroid().z());

                                Eigen::Vector3d M2(volume_and_surface2.getMoments().centroid().x(),
                                                volume_and_surface2.getMoments().centroid().y(),
                                                volume_and_surface2.getMoments().centroid().z());


                                if (Vdiff <= 0.0) {
                                    vfrac[i][j][k] = 0.0;
                                    firstMoment[i][j][k] = Eigen::Vector3d::Zero();
                                } else {
                                    vfrac[i][j][k] = Vdiff / cell_volume;
                                    firstMoment[i][j][k] = M2 - M1;
                                }
                            }

                            if (visualize) {
                                auto surface1 = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid1).getSurface();
                                auto surface2 = getVolumeMoments<
                                    VolumeMomentsAndSurface>(cell, paraboloid2).getSurface();
                                //surfaces.push_back(volume_and_surface1.getSurface());
                                //surfaces.push_back(volume_and_surface2.getSurface());
                                surfaces.push_back(surface1);
                                surfaces.push_back(surface2);
                                centroid[i][j][k] = computeCentroidFromFirstMoment(
                                    firstMoment[i][j][k], vfrac[i][j][k] * cell_volume);
                            }
                        }
                    }
                }

                // Check central cell
                if (centerCellIsCut(vfrac, stencil_size, machineZero)) {
                    // Accept this sample
                    
                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        *secondMoment = centeredSecondMomentFromTotal(totalMoments);
                    }

                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    // Flatten the 3D vector vfrac into a 1D vector
                    return flattenVfrac(vfrac, stencil_size);
                }

                // else: reject and regenerate
            }
        }

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

        std::vector<double> generateCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double max_radius = 0.5, double radius_stddev = 0.0, bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr) 
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
                //TESTING 
                //axis_direction = Eigen::Vector3d(1.0, 0.0, 0.0);
                //axis_origin = Eigen::Vector3d(0.0, 0.0, 0.0);

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

                IRL::GeneralMoments3D<2> totalMoments_refined =
                    IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For exact 2nd moment shift later

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

                            if (secondMoment != nullptr) {
                                // Exact raw 0th/1st/2nd moments about global origin, accumulated on refined grid
                                auto gm = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid);
                                totalMoments_refined += gm;
                            }

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

                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        // liquid-centered (about global liquid centroid) second moment matrix from refined accumulated moments
                        *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                    }

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

        inline Eigen::Vector3d centerlinePoint(
            double theta,
            const Eigen::Vector3d& c0,
            double Rc,
            const Eigen::Vector3d& u,
            const Eigen::Vector3d& v)
        {
            return c0 + Rc * (std::cos(theta) * u + std::sin(theta) * v);
        }

        inline void closestPointAndTangentOnCircle(
            const Eigen::Vector3d& p,
            const Eigen::Vector3d& c0,
            double Rc,
            const Eigen::Vector3d& u,
            const Eigen::Vector3d& v,
            Eigen::Vector3d& closest,
            Eigen::Vector3d& tangent_unit)
        {
            // Vector from circle center to point
            Eigen::Vector3d q = p - c0;

            // Project onto circle plane basis
            double qu = q.dot(u);
            double qv = q.dot(v);

            // Closest angle on full circle
            double theta_star = std::atan2(qv, qu);

            // Closest point on circle
            closest = centerlinePoint(theta_star, c0, Rc, u, v);

            // Tangent direction (already unit if u and v are orthonormal)
            tangent_unit = (-std::sin(theta_star) * u + std::cos(theta_star) * v);
        }


        std::vector<double> generateBentCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double max_radius = 0.5, double radius_stddev = 0.0,
            double radius_circle_min = 2.5, double radius_circle_max = 10.0,
            bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr)
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

                // Make a random tube radius
                double tube_radius = max_radius;
                if (radius_stddev > 0.0) {
                    tube_radius = sample_truncated_normal(0, radius_stddev, machineZero, max_radius);
                } else {
                    std::uniform_real_distribution<double> random_thickness(machineZero, max_radius);
                    tube_radius = random_thickness(eng);
                }

                // Random plane (u, v) in which the arc lies
                Eigen::Vector3d u = generateRandomDirection(eng).normalized();
                Eigen::Vector3d tmp = generateRandomDirection(eng).normalized();
                tmp -= tmp.dot(u) * u; // make tmp orthogonal to u
                if (tmp.squaredNorm() < 1e-14) {
                    // rare degeneracy: choose a deterministic perpendicular
                    tmp = (std::abs(u.x()) < 0.9 ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0));
                    tmp -= tmp.dot(u) * u;
                }
                Eigen::Vector3d v = tmp.normalized();

                // plane normal (unit)
                Eigen::Vector3d plane_normal = u.cross(v).normalized();

                // Approximating the unit cell as a sphere of radius equal to half the diagonal of the unit cube, which is ~0.8661
                const double half_diag = 0.8661;

                // Make the plane from two directions and a random offset h along the normal. Sample h in a range that ensures the plane intersects the "spherical" cell
                std::uniform_real_distribution<double> dist_h(-(half_diag + tube_radius), (half_diag + tube_radius));
                double h = dist_h(eng);
                double plane_distance_from_origin = std::abs(h);

                // Point on plane closest to origin (the foot point)
                Eigen::Vector3d p_plane_closest = h * plane_normal;

                // Calculate the radius of a circle that would exist if you cut a sphere centered at 0,0,0 with radius 0.8661 at that height above 0.
                double cut_central_cell_radius = 0.0;
                {
                    // sphere radius = diag, plane at signed distance h from origin
                    double inside = half_diag * half_diag - plane_distance_from_origin * plane_distance_from_origin;
                    if (inside > 0.0) {
                        cut_central_cell_radius = std::sqrt(inside);
                    } else {
                        continue; // plane doesn't intersect the "spherical" cell, try again
                    }
                }

                // Make a random circle_radius 
                std::uniform_real_distribution<double> dist_radius_circle(radius_circle_min, radius_circle_max);
                double radius_circle = dist_radius_circle(eng);

                // A random direction in the plane for the circle to be centered along. This controls the "bending direction" of the tube
                Eigen::Vector3d circle_center_direction = generateRandomDirection(eng).normalized();
                // project onto plane
                Eigen::Vector3d dir_plane = circle_center_direction - circle_center_direction.dot(plane_normal) * plane_normal;
                double dir_plane_norm2 = dir_plane.squaredNorm();
                if (dir_plane_norm2 < 1e-14) {
                    // if circle_center_direction almost parallel to plane normal, just use u as an in-plane direction
                    dir_plane = u;
                } else {
                    dir_plane /= std::sqrt(dir_plane_norm2);
                }

                // Along this direction, choose a distance d from the foot point to place the center of the circle (c0)
                // Sample d in a range that ensures the tube of radius tube_radius around the arc will intersect the central "spherical" cell
                double d_min = radius_circle - (tube_radius + cut_central_cell_radius);
                double d_max = radius_circle + tube_radius + cut_central_cell_radius;
                std::uniform_real_distribution<double> dist_d(d_min, d_max);
                double d = dist_d(eng);

                // Final circle center (c0)
                Eigen::Vector3d c0 = p_plane_closest + d * dir_plane;

                // Quick reject: check distance from stencil center to centerline (arc segment)
                /*
                Eigen::Vector3d center(0.0, 0.0, 0.0);
                Eigen::Vector3d closest_center, tangent_center;
                closestPointAndTangentOnArc(center, closest_center, tangent_center);
                double distance_to_centerline = (closest_center - center).norm();

                if (std::abs(distance_to_centerline - radius) > 0.8661) {
                    continue; // try again
                }
                */

                // Refined mesh
                int refined_stencil_size = refinement_factor * stencil_size;

                std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                    std::vector<std::vector<double>>(refined_stencil_size,
                    std::vector<double>(refined_stencil_size)));

                std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                    std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

                // Refined cell coordinates
                auto coords = std::vector<double>(refined_stencil_size + 1);
                for (int i = 0; i <= refined_stencil_size; i++) {
                    coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                }

                // Fill refined stencil
                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

                IRL::GeneralMoments3D<2> totalMoments_refined =
                    IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For exact 2nd moment

                for (int i = 0; i < refined_stencil_size; i++) {
                    for (int j = 0; j < refined_stencil_size; j++) {
                        for (int k = 0; k < refined_stencil_size; k++) {
                            auto cell = IRL::RectangularCuboid::fromBoundingPts(
                                IRL::Pt(coords[i], coords[j], coords[k]),
                                IRL::Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                            // Find center of cell
                            Eigen::Vector3d cell_center((coords[i + 1] + coords[i]) / 2.0,
                                                        (coords[j + 1] + coords[j]) / 2.0,
                                                        (coords[k + 1] + coords[k]) / 2.0);

                            // Find closest point on circle + tangent at that point
                            Eigen::Vector3d closest_point_on_circle;
                            Eigen::Vector3d circle_tangent;
                            closestPointAndTangentOnCircle(
                                cell_center,
                                c0,
                                radius_circle,
                                u,
                                v,
                                closest_point_on_circle,
                                circle_tangent
                            );


                            // Local radial direction from circle to cell center
                            Eigen::Vector3d radial_dir = (cell_center - closest_point_on_circle);
                            double rn = radial_dir.norm();
                            if (rn < 1e-14) {
                                // pick any direction perpendicular to tangent
                                radial_dir = u - u.dot(circle_tangent) * circle_tangent;
                                if (radial_dir.squaredNorm() < 1e-14)
                                    radial_dir = v - v.dot(circle_tangent) * circle_tangent;
                                radial_dir.normalize();
                            } else {
                                radial_dir.normalize();
                            }

                            // Datum on paraboloid representation of cylinder (local)
                            Eigen::Vector3d datum_paraboloid_eVec = closest_point_on_circle + tube_radius * radial_dir;
                            IRL::Pt datum_paraboloid(datum_paraboloid_eVec.x(),
                                                    datum_paraboloid_eVec.y(),
                                                    datum_paraboloid_eVec.z());

                            // Frame aligned with local tangent
                            Eigen::Vector3d paraboloid_x = circle_tangent.normalized(); // local axis
                            Eigen::Vector3d paraboloid_z = radial_dir; // local radial
                            Eigen::Vector3d paraboloid_y = paraboloid_z.cross(paraboloid_x);
                            double y2 = paraboloid_y.squaredNorm();
                            if (y2 < 1e-14) {
                                // Fallback if radial accidentally parallel to tangent
                                Eigen::Vector3d cand = u - u.dot(paraboloid_x) * paraboloid_x;
                                if (cand.squaredNorm() < 1e-14) cand = v - v.dot(paraboloid_x) * paraboloid_x;
                                paraboloid_z = cand.normalized();
                                paraboloid_y = paraboloid_z.cross(paraboloid_x);
                                paraboloid_y.normalize();
                            } else {
                                paraboloid_y /= std::sqrt(y2);
                            }

                            const auto frame = IRL::ReferenceFrame(
                                IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()),
                                IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()),
                                IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));

                            // Tube = paraboloid with coeffs (0, 1/(2R)) in the local frame
                            const auto paraboloid = IRL::Paraboloid(datum_paraboloid, frame, 0, 1 / (2 * tube_radius));

                            auto volume_and_surface = getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid);

                            volumes_refined[i][j][k] = volume_and_surface.getMoments().volume();
                            firstMoments_refined[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                            volume_and_surface.getMoments().centroid().y(),
                                                            volume_and_surface.getMoments().centroid().z();

                            if (secondMoment != nullptr) {
                                // Exact raw 0th/1st/2nd moments about global origin, accumulated on refined grid
                                auto gm = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid);
                                totalMoments_refined += gm;
                            }

                            if (visualize) {
                                auto surface = volume_and_surface.getSurface();
                                surfaces.push_back(surface);
                            }
                        }
                    }
                }

                // Compress refined to coarse stencil
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
                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {

                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        // liquid-centered (about global liquid centroid) second moment matrix from refined accumulated moments
                        *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                    }

                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    std::vector<double> flattened_vfrac;
                    for (int ii = 0; ii < stencil_size; ++ii) {
                        for (int jj = 0; jj < stencil_size; ++jj) {
                            for (int kk = 0; kk < stencil_size; ++kk) {
                                flattened_vfrac.push_back(vfrac[ii][jj][kk]);
                            }
                        }
                    }
                    return flattened_vfrac; // accept this sample
                }
                // else: reject and regenerate
            }
        }

        inline double wrapToPi(double a) {
            while (a <= -M_PI) a += 2.0 * M_PI;
            while (a >   M_PI) a -= 2.0 * M_PI;
            return a;
        }

        // Shift angle a by ±2π so that it's as close as possible to ref.
        inline double unwrapNear(double a, double ref) {
            return ref + wrapToPi(a - ref);
        }

        inline double angleOnCirclePlane(const Eigen::Vector3d& p,
                                        const Eigen::Vector3d& c0,
                                        const Eigen::Vector3d& u,
                                        const Eigen::Vector3d& v)
        {
            Eigen::Vector3d q = p - c0;
            double qu = q.dot(u);
            double qv = q.dot(v);
            return std::atan2(qv, qu);
        }

        std::vector<double> generateBentTruncatedCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, bool truncateInsideCentralCell, double max_radius = 0.5, double radius_stddev = 0.0, 
            double radius_circle_min = 2.5, double radius_circle_max = 10.0,
            bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr)
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

                // Random tube radius in [0, max_radius], with optional stddev for truncated normal sampling
                double tube_radius = max_radius;
                if (radius_stddev > 0.0) {
                    tube_radius = sample_truncated_normal(0, radius_stddev, machineZero, max_radius);
                } else {
                    std::uniform_real_distribution<double> random_thickness(machineZero, max_radius);
                    tube_radius = random_thickness(eng);
                }

                // random circle plane (u,v)
                Eigen::Vector3d u = generateRandomDirection(eng).normalized();
                Eigen::Vector3d tmp = generateRandomDirection(eng).normalized();
                tmp -= tmp.dot(u) * u;
                if (tmp.squaredNorm() < 1e-14) {
                    tmp = (std::abs(u.x()) < 0.9 ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0));
                    tmp -= tmp.dot(u) * u;
                }
                Eigen::Vector3d v = tmp.normalized();
                Eigen::Vector3d plane_normal = u.cross(v).normalized();

                // plane offset h near origin
                const double half_diag = 0.8661; // diagonal distance from center to corner of unit cube

                std::uniform_real_distribution<double> dist_h(-(half_diag + tube_radius), (half_diag + tube_radius));
                double h = dist_h(eng);
                double plane_dist = std::abs(h);
                Eigen::Vector3d p_plane_closest = h * plane_normal;

                double cut_central_cell_radius = 0.0;
                {
                    double inside = half_diag * half_diag - plane_dist * plane_dist;
                    if (inside > 0.0) cut_central_cell_radius = std::sqrt(inside);
                    else continue; // plane misses the central "sphere"
                }

                // circle/bend radius
                double radius_circle_min = 2.5;
                double radius_circle_max = 10.0;
                std::uniform_real_distribution<double> dist_radius_circle(radius_circle_min, radius_circle_max);
                double radius_circle = dist_radius_circle(eng);

                // place c0 in the plane along a random in-plane direction
                Eigen::Vector3d circle_center_direction = generateRandomDirection(eng).normalized();
                Eigen::Vector3d dir_plane = circle_center_direction - circle_center_direction.dot(plane_normal) * plane_normal;
                if (dir_plane.squaredNorm() < 1e-14) dir_plane = u;
                else dir_plane.normalize();

                double d_min = radius_circle - (tube_radius + cut_central_cell_radius);
                double d_max = radius_circle + tube_radius + cut_central_cell_radius;
                std::uniform_real_distribution<double> dist_d(d_min, d_max);
                double d = dist_d(eng);

                Eigen::Vector3d c0 = p_plane_closest + d * dir_plane;

                // create sphere radius
                double sphere_radius = tube_radius;
                {
                    double sphere_radius_min = tube_radius;
                    double sphere_radius_max = max_radius;
                    std::uniform_real_distribution<double> dist_sr(sphere_radius_min, sphere_radius_max);
                    sphere_radius = dist_sr(eng);
                }

                // pick random theta_cut on full circle and accept only if it satisfies your constraints
                double theta_cut = 0.0;
                Eigen::Vector3d sphere_origin = Eigen::Vector3d::Zero();

                // distance from origin to closest point on the circle (d_min_circle)
                // (use closest-point angle for origin)
                {
                    Eigen::Vector3d O(0.0, 0.0, 0.0);
                    double theta_star_O = angleOnCirclePlane(O, c0, u, v);
                    Eigen::Vector3d closest_O = centerlinePoint(theta_star_O, c0, radius_circle, u, v);
                    double d_min_circle = closest_O.norm();

                    // Try a few theta samples; if none work, restart outer while(true)
                    bool found = false;
                    std::uniform_real_distribution<double> dist_theta(-M_PI, M_PI);
                    for (int tries = 0; tries < 64; ++tries) {
                        double t = dist_theta(eng);
                        Eigen::Vector3d p_cut = centerlinePoint(t, c0, radius_circle, u, v);

                        // within stencil + sphere_radius (use infinity norm / bounding cube check)
                        double half = 0.5 * static_cast<double>(stencil_size);
                        double maxAbs = std::max({std::abs(p_cut.x()), std::abs(p_cut.y()), std::abs(p_cut.z())});
                        if (maxAbs > half + sphere_radius) continue;

                        // If truncateInsideCentralCell == false: require OUTSIDE  => ||p_cut|| > d_min_circle + sphere_radius
                        // If truncateInsideCentralCell == true : require INSIDE   => ||p_cut|| <= d_min_circle + sphere_radius
                        const double thresh = d_min_circle + sphere_radius;
                        const double r_cut  = p_cut.norm();
                        if (!truncateInsideCentralCell) {
                            if (r_cut <= thresh) continue;     // want outside, but it's inside
                        } else {
                            if (r_cut >  thresh) continue;     // want inside, but it's outside
                        }

                        theta_cut = t;
                        sphere_origin = p_cut;
                        found = true;
                        break;
                    }
                    if (!found) continue;
                }

                // range of theta for tube in [theta_start, theta_cut]
                // Randomize span to get diverse truncation lengths.
                double tube_theta_span = 0.0;
                {
                    std::uniform_real_distribution<double> dist_span(0.5 * M_PI, 1.75 * M_PI);
                    tube_theta_span = dist_span(eng);
                }
                double theta_start = theta_cut - tube_theta_span;

                auto thetaInTubeRange = [&](double theta) -> bool {
                    // unwrap theta close to theta_cut so interval comparison works
                    double t = unwrapNear(theta, theta_cut);
                    return (t >= theta_start && t <= theta_cut);
                };

                // --- refined mesh ---
                int refined_stencil_size = refinement_factor * stencil_size;

                std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                    std::vector<std::vector<double>>(refined_stencil_size,
                    std::vector<double>(refined_stencil_size)));

                std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                    std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

                // Refined cell coordinates
                auto coords = std::vector<double>(refined_stencil_size + 1);
                for (int i = 0; i <= refined_stencil_size; i++) {
                    coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                }

                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

                IRL::GeneralMoments3D<2> totalMoments_refined =
                    IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);

                for (int i = 0; i < refined_stencil_size; i++) {
                    for (int j = 0; j < refined_stencil_size; j++) {
                        for (int k = 0; k < refined_stencil_size; k++) {

                            auto cell = IRL::RectangularCuboid::fromBoundingPts(
                                IRL::Pt(coords[i], coords[j], coords[k]),
                                IRL::Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                            Eigen::Vector3d cell_center((coords[i + 1] + coords[i]) / 2.0,
                                                        (coords[j + 1] + coords[j]) / 2.0,
                                                        (coords[k + 1] + coords[k]) / 2.0);

                            // Compute theta of this cell projection in the circle plane
                            double theta_cell = angleOnCirclePlane(cell_center, c0, u, v);

                            // Always compute closest point on circle for tube distance and tube frame if needed
                            Eigen::Vector3d closest_point_on_circle;
                            Eigen::Vector3d circle_tangent;
                            closestPointAndTangentOnCircle(
                                cell_center,
                                c0,
                                radius_circle,
                                u,
                                v,
                                closest_point_on_circle,
                                circle_tangent
                            );

                            double dist_to_circle_centerline = (cell_center - closest_point_on_circle).norm();
                            double d2 = dist_to_circle_centerline - tube_radius; // tube signed-ish

                            double d1 = (cell_center - sphere_origin).norm() - sphere_radius; // sphere signed-ish

                            bool use_sphere = false;

                            // Rule 1: theta outside tube range => sphere
                            if (!thetaInTubeRange(theta_cell)) {
                                use_sphere = true;
                            } else {
                                // Rule 2: compare distance-to-surface values
                                if (d1 < d2) use_sphere = true;
                                else use_sphere = false;
                            }

                            // Build the surface for this cell
                            IRL::Paraboloid paraboloid = IRL::Paraboloid(
                                IRL::Pt(0,0,0),
                                IRL::ReferenceFrame(IRL::Normal(1,0,0), IRL::Normal(0,1,0), IRL::Normal(0,0,1)),
                                0, 0
                            );

                            if (use_sphere) {
                                // --- Sphere interface ---
                                Eigen::Vector3d z = (cell_center - sphere_origin);
                                double zn = z.norm();
                                if (zn < 1e-14) {
                                    z = Eigen::Vector3d(1, 0, 0);
                                } else {
                                    z /= zn;
                                }

                                // pick a helper not parallel to z
                                Eigen::Vector3d helper = (std::abs(z.x()) < 0.9 ? Eigen::Vector3d(1,0,0) : Eigen::Vector3d(0,1,0));
                                Eigen::Vector3d x = z.cross(helper);
                                double xn = x.norm();
                                if (xn < 1e-14) {
                                    helper = Eigen::Vector3d(0,0,1);
                                    x = z.cross(helper);
                                    xn = x.norm();
                                }
                                x /= xn;
                                Eigen::Vector3d y = z.cross(x).normalized();

                                Eigen::Vector3d datum_e = sphere_origin + sphere_radius * z;
                                IRL::Pt datum(datum_e.x(), datum_e.y(), datum_e.z());

                                const auto frame = IRL::ReferenceFrame(
                                    IRL::Normal(x.x(), x.y(), x.z()),
                                    IRL::Normal(y.x(), y.y(), y.z()),
                                    IRL::Normal(z.x(), z.y(), z.z())
                                );

                                paraboloid = IRL::Paraboloid(datum, frame,
                                                            1.0 / (2.0 * sphere_radius),
                                                            1.0 / (2.0 * sphere_radius));
                            } else {
                                // --- Tube interface (local cylinder) ---
                                Eigen::Vector3d radial_dir = (cell_center - closest_point_on_circle);
                                double rn = radial_dir.norm();
                                if (rn < 1e-14) {
                                    radial_dir = u - u.dot(circle_tangent) * circle_tangent;
                                    if (radial_dir.squaredNorm() < 1e-14)
                                        radial_dir = v - v.dot(circle_tangent) * circle_tangent;
                                    radial_dir.normalize();
                                } else {
                                    radial_dir /= rn;
                                }

                                Eigen::Vector3d x = circle_tangent.normalized();
                                Eigen::Vector3d z = radial_dir;
                                Eigen::Vector3d y = z.cross(x);
                                double y2 = y.squaredNorm();
                                if (y2 < 1e-14) {
                                    Eigen::Vector3d cand = u - u.dot(x) * x;
                                    if (cand.squaredNorm() < 1e-14) cand = v - v.dot(x) * x;
                                    z = cand.normalized();
                                    y = z.cross(x).normalized();
                                } else {
                                    y /= std::sqrt(y2);
                                }

                                Eigen::Vector3d datum_e = closest_point_on_circle + tube_radius * z;
                                IRL::Pt datum(datum_e.x(), datum_e.y(), datum_e.z());

                                const auto frame = IRL::ReferenceFrame(
                                    IRL::Normal(x.x(), x.y(), x.z()),
                                    IRL::Normal(y.x(), y.y(), y.z()),
                                    IRL::Normal(z.x(), z.y(), z.z())
                                );

                                paraboloid = IRL::Paraboloid(datum, frame, 0.0, 1.0 / (2.0 * tube_radius));
                            }

                            auto volume_and_surface = getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid);

                            volumes_refined[i][j][k] = volume_and_surface.getMoments().volume();
                            firstMoments_refined[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                            volume_and_surface.getMoments().centroid().y(),
                                                            volume_and_surface.getMoments().centroid().z();

                            if (secondMoment != nullptr) {
                                auto gm = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid);
                                totalMoments_refined += gm;
                            }

                            if (visualize) {
                                auto surface = volume_and_surface.getSurface();
                                surfaces.push_back(surface);
                            }
                        }
                    }
                }

                // Compress refined to coarse stencil
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
                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {

                    if (secondMoment != nullptr) {
                        *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                    }

                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }

                    std::vector<double> flattened_vfrac;
                    flattened_vfrac.reserve(stencil_size * stencil_size * stencil_size);
                    for (int ii = 0; ii < stencil_size; ++ii)
                        for (int jj = 0; jj < stencil_size; ++jj)
                            for (int kk = 0; kk < stencil_size; ++kk)
                                flattened_vfrac.push_back(vfrac[ii][jj][kk]);

                    return flattened_vfrac;
                }

                // else: reject and regenerate
            }
        }


        std::vector<double> generateTruncatedCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double max_radius = 0.5, double radius_stddev = 0.0, bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr) 
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
                IRL::GeneralMoments3D<2> totalMoments_refined =
                    IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For exact 2nd moment shift later

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

                            if (secondMoment != nullptr) {
                                // Exact raw 0th/1st/2nd moments about global origin, accumulated on refined grid
                                auto gm = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid);
                                totalMoments_refined += gm;
                            }

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
                                        // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        // liquid-centered (about global liquid centroid) second moment matrix from refined accumulated moments
                        *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                    }

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
            int stencil_size, double max_radius = 0.5, double radius_stddev = 0.0, bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr) 
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
                    IRL::GeneralMoments3D<2> totalMoments_refined =
                        IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For exact 2nd moment shift later

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

                                if (secondMoment != nullptr) {
                                    // Exact raw 0th/1st/2nd moments about global origin, accumulated on refined grid
                                    auto gm = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid);
                                    totalMoments_refined += gm;
                                }

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
                        
                        // Now calc stencil 2nd moments if requested
                        if (secondMoment != nullptr) {
                            // liquid-centered (about global liquid centroid) second moment matrix from refined accumulated moments
                            *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                        }
                        
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
                                        bool visualize = false, bool exact_2nd_moments = false)
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

            Eigen::Matrix3d secondMoment;
            Eigen::Matrix3d* secondMomentPtr = nullptr;    // default: disabled using exact 2nd moment

            if (exact_2nd_moments) {
                secondMomentPtr = &secondMoment; // This enables calculation of exact 2nd moments for all switch cases below
            }
            std::cout<<"Generating data point of type "<<datapoint_type<<std::endl;
            std::cout<<"Visualize: "<<visualize<<std::endl;

            // Testing
            bool truncateInsideCentrCell = true;
            double radius_circ_min = 2.5;
            double radius_circ_max = 10.0;


            switch (datapoint_type) {
                case 1:
                    if (include_truncated_cylinder) {
                        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
                        double p = prob_dist(eng);  // draw a random number in [0,1)

                        if (p < 0.2) {
                            // 20% chance → truncated cylinder
                            generateTruncatedCylinder(
                                vfrac, firstMoment, stencil_size,
                                max_cylinder_radius, cylinder_radius_stddev, visualize, &secondMoment);
                        } else {
                            // 80% chance → normal cylinder
                            generateCylinder(
                                vfrac, firstMoment, stencil_size,
                                max_cylinder_radius, cylinder_radius_stddev, visualize, &secondMoment);
                        }
                    } else {
                        generateCylinder(vfrac, firstMoment, stencil_size, max_cylinder_radius, cylinder_radius_stddev, visualize, &secondMoment);
                    }
                    break;
                case 2:
                    generateSphere(vfrac, firstMoment, stencil_size, max_sphere_radius, sphere_radius_stddev, visualize, &secondMoment);
                    break;
                case 3:
                    generateSheet(vfrac, firstMoment, stencil_size, sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev, visualize, &secondMoment);
                    break;
                case 4:
                    std::cout<<"Generating Cut Sheet"<<std::endl;
                    //generateBentCylinder(vfrac, firstMoment, stencil_size, max_cylinder_radius, cylinder_radius_stddev, visualize, &secondMoment);
                    generateBentTruncatedCylinder(vfrac, firstMoment, stencil_size, truncateInsideCentrCell, max_cylinder_radius, cylinder_radius_stddev, radius_circ_min, radius_circ_max, visualize, &secondMoment);
                    //generateCutSheet(vfrac, firstMoment, stencil_size, sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev, visualize, &secondMoment);
                    break;
                default:
                    generateParaboloid(vfrac, firstMoment, stencil_size, paraboloid_coeff_stddev, visualize, &secondMoment);
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
                if (exact_2nd_moments) {
                    // Append exact 2nd moments to flattened_state
                    flattened_state.push_back(secondMoment(0, 0)); // Ixx
                    flattened_state.push_back(secondMoment(1, 1)); // Iyy
                    flattened_state.push_back(secondMoment(2, 2)); // Izz
                    flattened_state.push_back(secondMoment(0, 1)); // Ixy
                    flattened_state.push_back(secondMoment(0, 2)); // Ixz
                    flattened_state.push_back(secondMoment(1, 2)); // Iyz
                } else {
                    // Append approximate 2nd moments to flattened_state
                    Eigen::Matrix3d approxSecondMoment = IRL::compute2ndMoment(flattened_state, stencil_size, 1);
                    flattened_state.push_back(approxSecondMoment(0, 0)); // Ixx
                    flattened_state.push_back(approxSecondMoment(1, 1)); // Iyy
                    flattened_state.push_back(approxSecondMoment(2, 2)); // Izz
                    flattened_state.push_back(approxSecondMoment(0, 1)); // Ixy
                    flattened_state.push_back(approxSecondMoment(0, 2)); // Ixz
                    flattened_state.push_back(approxSecondMoment(1, 2)); // Iyz
                }
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

            if (include_Moments >= 4) {
                // This is for testing: compute approximate 2nd moments from flattened_state and compare to exact 2nd moments

                // Compare exact vs approximate 2nd moments
                Eigen::Matrix3d secondMomentFrom1 = IRL::compute2ndMoment(flattened_state, stencil_size, 1);
                
                double frobeniusError = (secondMoment - secondMomentFrom1).norm() / (secondMoment.norm() + 1e-12);
                // Since this in only used for comparison, return a vector with only the fabrenius error as component
                return std::vector<float>(1, static_cast<float>(frobeniusError));
                

                /*
                flattened_state.push_back(inertia(0, 0)); // Ixx
                flattened_state.push_back(inertia(1, 1)); // Iyy
                flattened_state.push_back(inertia(2, 2)); // Izz
                flattened_state.push_back(inertia(0, 1)); // Ixy
                flattened_state.push_back(inertia(0, 2)); // Ixz
                flattened_state.push_back(inertia(1, 2)); // Iyz
                */
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