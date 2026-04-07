#pragma once

#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <array>
#include <numeric>
#include <stdexcept>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/ml_classification/vtk_out.h"
#include "irl/ml_classification/inertia_calc.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

namespace IRL {
    class Data_gen {

        public:
        std::mt19937_64 eng;
        static constexpr double machineZero = 1e-12;
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
        
        
        void generateParaboloid(
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
                //double coeff1 = random_coeff(eng);
                //double coeff2 = random_coeff(eng);
                double coeff1 = sample_truncated_normal(0.0, coeff_stddev, -1.0, 1.0);
                double coeff2 = sample_truncated_normal(0.0, coeff_stddev, -1.0, 1.0);

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
                    return; // done with this function, exit
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


        void generateSheet(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double coeff_stddev = 0.1, double min_thickness = machineZero, double max_thickness = 0.5, double thickness_stddev = 0.0, bool visualize = false,
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
                    thickness = sample_truncated_normal(0, thickness_stddev, min_thickness, max_thickness);
                }else{
                    std::uniform_real_distribution<double> random_thickness(min_thickness, max_thickness);
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
                double coeff1 = sample_truncated_normal(0.0, coeff_stddev, -1.0, 1.0);
                double coeff2 = sample_truncated_normal(0.0, coeff_stddev, -1.0, 1.0);

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
                    return;
                }

                // else: reject and regenerate
            }
        }
        
        void generateCutSheet(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, bool cutInsideCentralCell,
            double coeff_stddev = 0.1, double min_thickness = machineZero, double max_thickness = 0.5, double thickness_stddev = 0.0, bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr) 
        {
            while (true) { // keep trying until center cell has surface crossing

                // Defining cell coordinates
                auto coords = makeCenteredCoords(stencil_size);
                const double cell_volume = 1.0;

                // Random datum anywhere in stencil
                const auto datum = Pt::fromRawDoublePointer(generateRandomPoint(
                    -0.5*static_cast<double>(stencil_size),
                    0.5*static_cast<double>(stencil_size), eng).data());

                // Random unbiased direction
                auto direction = Normal::fromRawDoublePointer(generateRandomDirection(eng).data());
                direction.normalize();

                // Random sheet thickness
                double thickness = max_thickness;
                if (thickness_stddev > 0.0) {
                    thickness = sample_truncated_normal(0, thickness_stddev, min_thickness, max_thickness);
                }else{
                    std::uniform_real_distribution<double> random_thickness(min_thickness, max_thickness);
                    thickness = random_thickness(eng);
                }

                // Two paraboloid datums offset along the direction
                const auto datum_paraboloid1 = Pt(datum - direction * (thickness/2.0));
                const auto datum_paraboloid2 = Pt(datum + direction * (thickness/2.0));

                // Build orthonormal frame aligned with direction
                const auto frame = ReferenceFrame::fromNormal(direction);

                // Random coefficients
                std::normal_distribution<double> random_coeff(0.0, coeff_stddev);
                double coeff1 = random_coeff(eng);
                double coeff2 = random_coeff(eng);
                const auto paraboloid1 = Paraboloid(datum_paraboloid1, frame, coeff1, coeff2);
                const auto paraboloid2 = Paraboloid(datum_paraboloid2, frame, coeff1, coeff2);

                // Random point anywhere in stencil
                const auto sample_pt = Pt::fromRawDoublePointer(generateRandomPoint(
                    -0.5*static_cast<double>(stencil_size),
                    0.5*static_cast<double>(stencil_size), eng).data());

                // Project sample point onto paraboloid 2
                auto tmp_paraboloid = paraboloid2;
                tmp_paraboloid.regenerateAtLocation(sample_pt);
                const auto new_datum = tmp_paraboloid.getDatum();

                // Now choose a random theta
                std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
                double theta = angle_dist(eng);

                IRL::Normal new_normal =
                    std::cos(theta) * tmp_paraboloid.getReferenceFrame()[0] +
                    std::sin(theta) * tmp_paraboloid.getReferenceFrame()[1];

                new_normal.normalize();

                // Quick check if the central cell will be cut with this new plane, if not, reject and try again
                const double half_width = 0.5;                   // inscribed sphere radius of unit cube
                const double half_diag  = 0.86618;                // half diagonal of unit cube
                double dist = std::abs(new_normal * new_datum); 
                bool definitely_cut = (dist < half_width);
                bool definitely_not_cut = (dist > half_diag);
                if (definitely_cut && !cutInsideCentralCell) {
                    continue; // violates condition, reject early
                }
                if (definitely_not_cut && cutInsideCentralCell) {
                    continue; // violates condition, reject early
                }


                // Create localizer plane that will cut both paraboloids
                const auto localizer = PlanarLocalizer::fromOnePlane(Plane(new_normal, new_normal * new_datum));
                const auto paraboloid1_and_plane = LocalizedParaboloid<double>(&localizer, &paraboloid1);
                const auto paraboloid2_and_plane = LocalizedParaboloid<double>(&localizer, &paraboloid2);

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

                            if (secondMoment != nullptr) {
                                auto gm1 = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid1_and_plane);
                                auto gm2 = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid2_and_plane);
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
                                    VolumeMoments>(cell, paraboloid1_and_plane);

                                auto volume_and_surface2 = getVolumeMoments<
                                    VolumeMoments>(cell, paraboloid2_and_plane);

                                double V1 = volume_and_surface1.volume();
                                double V2 = volume_and_surface2.volume();
                                double Vdiff = V2 - V1;



                                Eigen::Vector3d M1(volume_and_surface1.centroid().x(),
                                                volume_and_surface1.centroid().y(),
                                                volume_and_surface1.centroid().z());

                                Eigen::Vector3d M2(volume_and_surface2.centroid().x(),
                                                volume_and_surface2.centroid().y(),
                                                volume_and_surface2.centroid().z());


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
                    }
                    return;
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

        void generateCylinder(
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
                const double cell_volume = 1.0;
                double max_refinement_factor = 6.0;
                double refinement_factor_double = std::ceil(3.0/(2.0*radius)); // want at least ~3 samples across the tube diameter for decent accuracy, can adjust this heuristic as needed
                int refinement_factor = static_cast<int>(refinement_factor_double);
                int refined_stencil_size = refinement_factor * stencil_size;

                // Refined cell coordinates
                auto coords = std::vector<double>(refined_stencil_size+1);
                for (int i = 0; i <= refined_stencil_size; i++) {
                    coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                }

                double totalV = 0.0;
                Eigen::Vector3d totalM1 = Eigen::Vector3d::Zero(); // raw first moment = ∑ V*c
                Eigen::Matrix3d totalM2 = Eigen::Matrix3d::Zero(); // raw second moment = ∑ ∫ x x^T dV

                // Also need to define total moments refined for refinement option
                IRL::GeneralMoments3D<2> totalMoments_refined =
                        IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For exact 2nd moment

                if (refinement_factor_double > max_refinement_factor) {
                    visualize = false;
                    const double max_line_segment_length = 0.25;

                    // Build coarse coordinates
                    std::vector<double> coords_coarse(stencil_size + 1);
                    for (int ii = 0; ii <= stencil_size; ++ii) {
                        coords_coarse[ii] = -0.5 * static_cast<double>(stencil_size) + static_cast<double>(ii);
                    }

                    // Clear coarse outputs
                    for (int i = 0; i < stencil_size; ++i)
                        for (int j = 0; j < stencil_size; ++j)
                            for (int k = 0; k < stencil_size; ++k) {
                                vfrac[i][j][k] = 0.0;
                                firstMoment[i][j][k].setZero();
                            }

                    // Helper: clip segment p(t)=p0+t*(p1-p0), t in [0,1], to cell AABB [bmin,bmax]
                    auto clipSegmentToCell = [&](const Eigen::Vector3d& p0,
                                                const Eigen::Vector3d& p1,
                                                const Eigen::Vector3d& bmin,
                                                const Eigen::Vector3d& bmax,
                                                double& ta, double& tb) -> bool
                    {
                        ta = 0.0;
                        tb = 1.0;
                        const Eigen::Vector3d d = p1 - p0;

                        for (int ax = 0; ax < 3; ++ax) {
                            const double p = d[ax];
                            const double q0 = p0[ax];

                            if (std::abs(p) < 1e-14) {
                                if (q0 < bmin[ax] || q0 > bmax[ax]) return false;
                            } else {
                                const double invp = 1.0 / p;
                                double tNear = (bmin[ax] - q0) * invp;
                                double tFar  = (bmax[ax] - q0) * invp;
                                if (tNear > tFar) std::swap(tNear, tFar);
                                ta = std::max(ta, tNear);
                                tb = std::min(tb, tFar);
                                if (ta > tb) return false;
                            }
                        }
                        return true;
                    };

                    // Central 2nd-moment tensor for a solid cylinder of volume V, length L, axis unit a
                    auto cylinderCentralSecondMoment = [&](double V, double L, const Eigen::Vector3d& a_unit) -> Eigen::Matrix3d {
                        const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
                        const Eigen::Matrix3d aaT = a_unit * a_unit.transpose();
                        const Eigen::Matrix3d P = I - aaT;
                        const double L2 = L * L;
                        const double r2 = radius * radius;
                        return (V * (L2 / 12.0)) * aaT + (V * (r2 / 4.0)) * P;
                    };

                    // Normalize cylinder axis direction
                    Eigen::Vector3d axis_unit = axis_direction.normalized();

                    // Build one long finite line segment along the axis, long enough to cover the whole stencil
                    const double stencil_diag = std::sqrt(3.0) * static_cast<double>(stencil_size);
                    const double total_line_length = 2.0 * stencil_diag + 2.0 * radius + 2.0;
                    const int nSeg = std::max(1, static_cast<int>(std::ceil(total_line_length / max_line_segment_length)));

                    const Eigen::Vector3d line_start = axis_origin - 0.5 * total_line_length * axis_unit;
                    const Eigen::Vector3d line_end   = axis_origin + 0.5 * total_line_length * axis_unit;

                    // Iterate coarse cells and accumulate clipped cylinder pieces
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {

                                const Eigen::Vector3d bmin(coords_coarse[i],   coords_coarse[j],   coords_coarse[k]);
                                const Eigen::Vector3d bmax(coords_coarse[i+1], coords_coarse[j+1], coords_coarse[k+1]);

                                double cellV = 0.0;
                                Eigen::Vector3d cellM1 = Eigen::Vector3d::Zero(); // raw first moment
                                Eigen::Matrix3d cellM2 = Eigen::Matrix3d::Zero(); // raw second moment

                                for (int s = 0; s < nSeg; ++s) {
                                    const double tA = static_cast<double>(s) / static_cast<double>(nSeg);
                                    const double tB = static_cast<double>(s + 1) / static_cast<double>(nSeg);

                                    const Eigen::Vector3d p0 = line_start + tA * (line_end - line_start);
                                    const Eigen::Vector3d p1 = line_start + tB * (line_end - line_start);

                                    double ta, tb;
                                    if (!clipSegmentToCell(p0, p1, bmin, bmax, ta, tb)) continue;

                                    // Part of the line segment inside the cell
                                    const Eigen::Vector3d q0 = p0 + ta * (p1 - p0);
                                    const Eigen::Vector3d q1 = p0 + tb * (p1 - p0);

                                    const Eigen::Vector3d dseg = q1 - q0;
                                    const double L = dseg.norm();
                                    if (L < 1e-14) continue;

                                    const Eigen::Vector3d axis_local = dseg / L;

                                    // Approximate the cell contribution as a solid cylinder piece
                                    const double V = M_PI * radius * radius * L;
                                    const Eigen::Vector3d c = 0.5 * (q0 + q1);

                                    const Eigen::Matrix3d C2   = cylinderCentralSecondMoment(V, L, axis_local);
                                    const Eigen::Matrix3d raw2 = C2 + V * (c * c.transpose());

                                    cellV  += V;
                                    cellM1 += V * c;
                                    cellM2 += raw2;
                                }

                                if (cellV > 0.0) {
                                    vfrac[i][j][k] = std::min(cellV / cell_volume, 1.0);
                                    firstMoment[i][j][k] = cellM1;  // raw first moment

                                    if (secondMoment != nullptr) {
                                        totalV  += cellV;
                                        totalM1 += cellM1;
                                        totalM2 += cellM2;
                                    }
                                } else {
                                    vfrac[i][j][k] = 0.0;
                                    firstMoment[i][j][k].setZero();
                                }
                            }
                        }
                    }
                } else {
                    // refine now
                    std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                        std::vector<std::vector<double>>(refined_stencil_size,
                        std::vector<double>(refined_stencil_size)));

                    std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                        std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                        std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

                    

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
                }    
                
                // --- check central cell ---
                int mid = stencil_size / 2;
                // check central cell
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {

                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        if (refinement_factor_double > max_refinement_factor) {
                            const Eigen::Vector3d xc = totalM1 / totalV;
                            const Eigen::Matrix3d central = totalM2 - totalV * (xc * xc.transpose());
                            *secondMoment = central; // or: central / totalV
                        } else {
                            // liquid-centered (about global liquid centroid) second moment matrix from refined accumulated moments
                            *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                        }
                    }

                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }
                    return; // accept this sample
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


        void generateBentCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double min_radius = machineZero, double max_radius = 0.5, double radius_stddev = 0.0,
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
                // Center cell
                int mid = stencil_size / 2;

                // Make a random tube radius
                double tube_radius = max_radius;
                if (radius_stddev > 0.0) {
                    tube_radius = sample_truncated_normal(0, radius_stddev, min_radius, max_radius);
                } else {
                    std::uniform_real_distribution<double> random_thickness(min_radius, max_radius);
                    tube_radius = random_thickness(eng);
                }

                // Random plane (u, v) in which the arc lies
                Eigen::Vector3d u = generateRandomDirection(eng).normalized();
                Eigen::Vector3d tmp = generateRandomDirection(eng).normalized();
                tmp -= tmp.dot(u) * u; // make tmp orthogonal to u: Gram–Schmidt orthogonalization
                if (tmp.squaredNorm() < 1e-14) {
                    // rare degeneracy: choose a deterministic perpendicular: x or y, then make orthogonal to u just in case
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

                // Refined mesh
                const double cell_volume = 1.0;
                double max_refinement_factor = 6.0;
                double refinement_factor_double = std::ceil(3.0/(2.0*tube_radius)); // want at least ~3 samples across the tube diameter for decent accuracy, can adjust this heuristic as needed
                int refinement_factor = static_cast<int>(refinement_factor_double);
                int refined_stencil_size = refinement_factor * stencil_size;
                // need to declare refined coords outside if below:
                // Refined cell coordinates
                auto coords = std::vector<double>(refined_stencil_size + 1);
                for (int i = 0; i <= refined_stencil_size; i++) {
                    coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                }
                // Also need to define total moments refined for refinement option
                IRL::GeneralMoments3D<2> totalMoments_refined =
                        IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For exact 2nd moment

                // If the tube radius is too small, approximate the circle by lines
                // need the following parameters in this scope first:
                // If we want stencil-level 2nd moments (liquid-centered), accumulate raw moments across the whole stencil
                double totalV = 0.0;
                Eigen::Vector3d totalM1 = Eigen::Vector3d::Zero(); // raw first moment = ∑ V*c
                Eigen::Matrix3d totalM2 = Eigen::Matrix3d::Zero(); // raw second moment = ∑ ∫ x x^T dV
                if (refinement_factor_double > max_refinement_factor) {
                    visualize = false;
                    const double max_line_segment_length = 0.25;

                    // Build coarse coordinates
                    std::vector<double> coords_coarse(stencil_size + 1);
                    for (int ii = 0; ii <= stencil_size; ++ii) {
                        coords_coarse[ii] = -0.5 * stencil_size + static_cast<double>(ii);
                    }

                    // polyline approximation of the circle centerline (full 2*pi)
                    const double circumference = 2.0 * M_PI * radius_circle;
                    const int nSeg = std::max(8, static_cast<int>(std::ceil(circumference / max_line_segment_length)));

                    auto circlePoint = [&](double theta) -> Eigen::Vector3d {
                        // circle lies in plane spanned by u,v with center c0
                        return c0 + radius_circle * (std::cos(theta) * u + std::sin(theta) * v);
                    };

                    // Clear coarse outputs
                    for (int i = 0; i < stencil_size; ++i)
                        for (int j = 0; j < stencil_size; ++j)
                            for (int k = 0; k < stencil_size; ++k) {
                                vfrac[i][j][k] = 0.0;
                                firstMoment[i][j][k].setZero();
                            }

                    
                    // What part of a line segment is within the cell AABB?
                    auto clipSegmentToCell = [&](const Eigen::Vector3d& p0,
                                                const Eigen::Vector3d& p1,
                                                const Eigen::Vector3d& bmin,
                                                const Eigen::Vector3d& bmax,
                                                double& t0, double& t1) -> bool
                    {
                        // Liang–Barsky style parametric clip in 3D for segment p(t)=p0+t*(p1-p0), t in [0,1]
                        t0 = 0.0; t1 = 1.0;
                        Eigen::Vector3d d = p1 - p0;

                        for (int ax = 0; ax < 3; ++ax) {
                            const double p = d[ax];
                            const double q0 = p0[ax];

                            if (std::abs(p) < 1e-14) {
                                // segment parallel to slab: must be within bounds
                                if (q0 < bmin[ax] || q0 > bmax[ax]) return false;
                            } else {
                                double invp = 1.0 / p;
                                double tNear = (bmin[ax] - q0) * invp;
                                double tFar  = (bmax[ax] - q0) * invp;
                                if (tNear > tFar) std::swap(tNear, tFar);
                                t0 = std::max(t0, tNear);
                                t1 = std::min(t1, tFar);
                                if (t0 > t1) return false;
                            }
                        }
                        return true;
                    };

                    auto cylinderCentralSecondMoment = [&](double V, double L, const Eigen::Vector3d& axis_unit) -> Eigen::Matrix3d {
                        // Central second moment (unnormalized): ∫ (x-c)(x-c)^T dV for a solid cylinder.
                        // For axis a, projector P = I - aa^T
                        // C2 = V*(L^2/12) * aa^T + V*(r^2/4) * P
                        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
                        Eigen::Matrix3d aaT = axis_unit * axis_unit.transpose();
                        Eigen::Matrix3d P = I - aaT;
                        const double L2 = L * L;
                        const double r2 = tube_radius * tube_radius;

                        return (V * (L2 / 12.0)) * aaT + (V * (r2 / 4.0)) * P;
                    };

                    // Iterate coarse cells and add contributions from clipped segments
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {

                                const Eigen::Vector3d bmin(coords_coarse[i],   coords_coarse[j],   coords_coarse[k]);
                                const Eigen::Vector3d bmax(coords_coarse[i+1], coords_coarse[j+1], coords_coarse[k+1]);

                                double cellV = 0.0;
                                Eigen::Vector3d cellM1 = Eigen::Vector3d::Zero(); // ∑ V*c within cell
                                Eigen::Matrix3d cellM2 = Eigen::Matrix3d::Zero(); // ∑ ∫ x x^T dV within cell

                                for (int s = 0; s < nSeg; ++s) {
                                    const double tA = (2.0 * M_PI) * (static_cast<double>(s) / nSeg);
                                    const double tB = (2.0 * M_PI) * (static_cast<double>(s + 1) / nSeg);

                                    const Eigen::Vector3d p0 = circlePoint(tA);
                                    const Eigen::Vector3d p1 = circlePoint(tB);

                                    // "Line passes through the cell" => centerline segment intersects the Cell AABB
                                    double ta, tb;
                                    if (!clipSegmentToCell(p0, p1, bmin, bmax, ta, tb)) continue;

                                    // Use only the part of the segment inside the cell
                                    Eigen::Vector3d q0 = p0 + ta * (p1 - p0);
                                    Eigen::Vector3d q1 = p0 + tb * (p1 - p0);

                                    const Eigen::Vector3d d = q1 - q0;
                                    const double L = d.norm();
                                    if (L < 1e-14) continue;

                                    const Eigen::Vector3d axis = d / L;

                                    // Solid cylinder approximation for this clipped piece
                                    const double V = M_PI * tube_radius * tube_radius * L;
                                    const Eigen::Vector3d c = 0.5 * (q0 + q1);

                                    // Central 2nd moment about centroid, then convert to raw about origin
                                    const Eigen::Matrix3d C2 = cylinderCentralSecondMoment(V, L, axis);
                                    const Eigen::Matrix3d raw2 = C2 + V * (c * c.transpose());

                                    cellV  += V;
                                    cellM1 += V * c;
                                    cellM2 += raw2;
                                }

                                // Convert to your stored per-cell format:
                                // vfrac = volume fraction (cell volume = 1)
                                // firstMoment stores centroid (your refinement code does that)
                                if (cellV > 0.0) {
                                    const double vf = std::min(cellV / cell_volume, 1.0);
                                    vfrac[i][j][k] = vf;

                                    // centroid of union-of-cylinders approximation in this cell
                                    firstMoment[i][j][k] = cellM1 / cellV;

                                    if (secondMoment != nullptr) {
                                        // accumulate second moments to get them for the whole stencil
                                        totalV  += cellV;
                                        totalM1 += cellM1;
                                        totalM2 += cellM2;
                                    }
                                } else {
                                    vfrac[i][j][k] = 0.0;
                                    firstMoment[i][j][k].setZero();
                                }
                            }
                        }
                    }
                } else {
                // need to refine now
                    std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                        std::vector<std::vector<double>>(refined_stencil_size,
                        std::vector<double>(refined_stencil_size)));

                    std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                        std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                        std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

                    // Fill refined stencil
                    using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

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

                                // Tube approx with paraboloid with coeffs (0, 1/(2R)) in the local frame- cylinder pieces
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
                }

                // check central cell
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {

                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        if (refinement_factor_double > max_refinement_factor) {
                            const Eigen::Vector3d xc = totalM1 / totalV;
                            const Eigen::Matrix3d central = totalM2 - totalV * (xc * xc.transpose());
                            *secondMoment = central; // or: central / totalV
                        } else {
                            // liquid-centered (about global liquid centroid) second moment matrix from refined accumulated moments
                            *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                        }
                    }

                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }
                    return; // accept this sample
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

        


        void generateTruncatedCylinder(
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
                const double cell_volume = 1.0;
                double max_refinement_factor = 6.0;
                double refinement_factor_double = std::ceil(3.0/(2.0*radius)); // want at least ~3 samples across the tube diameter for decent accuracy, can adjust this heuristic as needed
                int refinement_factor = static_cast<int>(refinement_factor_double);
                int refined_stencil_size = refinement_factor * stencil_size;

                // Refined cell coordinates
                auto coords = std::vector<double>(refined_stencil_size+1);
                for (int i = 0; i <= refined_stencil_size; i++) {
                    coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                }

                double totalV = 0.0;
                Eigen::Vector3d totalM1 = Eigen::Vector3d::Zero(); // raw first moment = ∑ V*c
                Eigen::Matrix3d totalM2 = Eigen::Matrix3d::Zero(); // raw second moment = ∑ ∫ x x^T dV

                // Also need to define total moments refined for refinement option
                IRL::GeneralMoments3D<2> totalMoments_refined =
                        IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For exact 2nd moment

                        
                if (refinement_factor_double > max_refinement_factor) {
                    visualize = false;
                    const double max_line_segment_length = 0.25;

                    // Build coarse coordinates
                    std::vector<double> coords_coarse(stencil_size + 1);
                    for (int ii = 0; ii <= stencil_size; ++ii) {
                        coords_coarse[ii] = -0.5 * static_cast<double>(stencil_size) + static_cast<double>(ii);
                    }

                    // Clear coarse outputs
                    for (int i = 0; i < stencil_size; ++i)
                        for (int j = 0; j < stencil_size; ++j)
                            for (int k = 0; k < stencil_size; ++k) {
                                vfrac[i][j][k] = 0.0;
                                firstMoment[i][j][k].setZero();
                            }

                    // Clip segment p(t)=p0+t*(p1-p0), t in [0,1], to cell AABB [bmin,bmax]
                    auto clipSegmentToCell = [&](const Eigen::Vector3d& p0,
                                                const Eigen::Vector3d& p1,
                                                const Eigen::Vector3d& bmin,
                                                const Eigen::Vector3d& bmax,
                                                double& ta, double& tb) -> bool
                    {
                        ta = 0.0;
                        tb = 1.0;
                        const Eigen::Vector3d d = p1 - p0;

                        for (int ax = 0; ax < 3; ++ax) {
                            const double p = d[ax];
                            const double q0 = p0[ax];

                            if (std::abs(p) < 1e-14) {
                                if (q0 < bmin[ax] || q0 > bmax[ax]) return false;
                            } else {
                                const double invp = 1.0 / p;
                                double tNear = (bmin[ax] - q0) * invp;
                                double tFar  = (bmax[ax] - q0) * invp;
                                if (tNear > tFar) std::swap(tNear, tFar);
                                ta = std::max(ta, tNear);
                                tb = std::min(tb, tFar);
                                if (ta > tb) return false;
                            }
                        }
                        return true;
                    };

                    // Central 2nd-moment tensor for a solid cylinder of volume V, length L, axis unit a
                    auto cylinderCentralSecondMoment = [&](double V, double L, const Eigen::Vector3d& a_unit) -> Eigen::Matrix3d {
                        const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
                        const Eigen::Matrix3d aaT = a_unit * a_unit.transpose();
                        const Eigen::Matrix3d P = I - aaT;
                        const double L2 = L * L;
                        const double r2 = radius * radius;
                        return (V * (L2 / 12.0)) * aaT + (V * (r2 / 4.0)) * P;
                    };

                    // Normalize axis direction
                    const Eigen::Vector3d axis_unit = axis_direction.normalized();

                    // Build one long finite line segment along the axis, long enough to cover the stencil
                    const double stencil_diag = std::sqrt(3.0) * static_cast<double>(stencil_size);
                    const double total_line_length = 2.0 * stencil_diag + 2.0 * radius + 2.0;
                    const int nSeg = std::max(1, static_cast<int>(std::ceil(total_line_length / max_line_segment_length)));

                    const Eigen::Vector3d line_start = axis_origin - 0.5 * total_line_length * axis_unit;
                    const Eigen::Vector3d line_end   = axis_origin + 0.5 * total_line_length * axis_unit;

                    // Iterate coarse cells and accumulate clipped cylinder pieces
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {

                                const Eigen::Vector3d bmin(coords_coarse[i],   coords_coarse[j],   coords_coarse[k]);
                                const Eigen::Vector3d bmax(coords_coarse[i+1], coords_coarse[j+1], coords_coarse[k+1]);

                                double cellV = 0.0;
                                Eigen::Vector3d cellM1 = Eigen::Vector3d::Zero(); // raw first moment
                                Eigen::Matrix3d cellM2 = Eigen::Matrix3d::Zero(); // raw second moment

                                for (int s = 0; s < nSeg; ++s) {
                                    const double tA = static_cast<double>(s) / static_cast<double>(nSeg);
                                    const double tB = static_cast<double>(s + 1) / static_cast<double>(nSeg);

                                    const Eigen::Vector3d p0 = line_start + tA * (line_end - line_start);
                                    const Eigen::Vector3d p1 = line_start + tB * (line_end - line_start);

                                    // First clip to cell
                                    double ta, tb;
                                    if (!clipSegmentToCell(p0, p1, bmin, bmax, ta, tb)) continue;

                                    Eigen::Vector3d q0 = p0 + ta * (p1 - p0);
                                    Eigen::Vector3d q1 = p0 + tb * (p1 - p0);

                                    // Then apply truncation at axis_origin:
                                    // keep only points with (x - axis_origin)·axis_unit >= 0
                                    double s0 = (q0 - axis_origin).dot(axis_unit);
                                    double s1 = (q1 - axis_origin).dot(axis_unit);

                                    if (s0 < 0.0 && s1 < 0.0) {
                                        continue; // fully behind truncation point
                                    }

                                    if (s0 < 0.0 || s1 < 0.0) {
                                        // Segment crosses truncation plane s=0
                                        const double alpha = s0 / (s0 - s1); // in [0,1]
                                        const Eigen::Vector3d qcut = q0 + alpha * (q1 - q0);

                                        if (s0 < 0.0) q0 = qcut;
                                        else          q1 = qcut;
                                    }

                                    const Eigen::Vector3d dseg = q1 - q0;
                                    const double L = dseg.norm();
                                    if (L < 1e-14) continue;

                                    const Eigen::Vector3d axis_local = dseg / L;

                                    // Solid cylinder approximation for this clipped piece
                                    const double V = M_PI * radius * radius * L;
                                    const Eigen::Vector3d c = 0.5 * (q0 + q1);

                                    const Eigen::Matrix3d C2   = cylinderCentralSecondMoment(V, L, axis_local);
                                    const Eigen::Matrix3d raw2 = C2 + V * (c * c.transpose());

                                    cellV  += V;
                                    cellM1 += V * c;
                                    cellM2 += raw2;
                                }

                                if (cellV > 0.0) {
                                    vfrac[i][j][k] = std::min(cellV / cell_volume, 1.0);
                                    firstMoment[i][j][k] = cellM1; // raw first moment

                                    if (secondMoment != nullptr) {
                                        totalV  += cellV;
                                        totalM1 += cellM1;
                                        totalM2 += cellM2;
                                    }
                                } else {
                                    vfrac[i][j][k] = 0.0;
                                    firstMoment[i][j][k].setZero();
                                }
                            }
                        }
                    }

                } else {
                    std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                        std::vector<std::vector<double>>(refined_stencil_size,
                        std::vector<double>(refined_stencil_size)));

                    std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                        std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                        std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));
                

                

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
                }

                int mid = stencil_size / 2;
                // check central cell
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {

                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        if (refinement_factor_double > max_refinement_factor) {
                            const Eigen::Vector3d xc = totalM1 / totalV;
                            const Eigen::Matrix3d central = totalM2 - totalV * (xc * xc.transpose());
                            *secondMoment = central; // or: central / totalV
                        } else {
                            // liquid-centered (about global liquid centroid) second moment matrix from refined accumulated moments
                            *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                        }
                    }

                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }
                    return; // accept this sample
                }
                // else: reject and regenerate
            } 
        }

        void generateSpecificSphere(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& centroid,
            std::vector<IRL::ParaboloidParametrizedSurfaceOutput>& surfaces,
            int stencil_size,
            const Eigen::Vector3d& origin,
            double radius,
            std::vector<double>& coarse_coords,
            bool visualize = false,
            Eigen::Matrix3d* secondMoment = nullptr
            )
        {
            // Analytical description for sphere case: if the sphere is fully contained within the a cell, we can directly set the volume fraction and moments 
            // Analytical description for sphere case: if the sphere is fully contained within ANY coarse cell,
            // we can directly set the volume fraction and moments (no refinement needed).

            // Stencil origin (coarse cell size = 1)
            const double coarse_stencil_min = -0.5 * static_cast<double>(stencil_size);

            // Identify the coarse cell that contains the sphere center (origin).
            // (If origin is outside stencil, clamp to valid range; containment test below will fail anyway.)
            auto clamp_coarse_index = [&](int idx) {
                return std::max(0, std::min(stencil_size - 1, idx));
            };

            int ci = clamp_coarse_index(static_cast<int>(std::floor(origin.x() - coarse_stencil_min)));
            int cj = clamp_coarse_index(static_cast<int>(std::floor(origin.y() - coarse_stencil_min)));
            int ck = clamp_coarse_index(static_cast<int>(std::floor(origin.z() - coarse_stencil_min)));

            // Bounds of that coarse cell
            const double x0 = coarse_stencil_min + static_cast<double>(ci);
            const double x1 = x0 + 1.0;
            const double y0 = coarse_stencil_min + static_cast<double>(cj);
            const double y1 = y0 + 1.0;
            const double z0 = coarse_stencil_min + static_cast<double>(ck);
            const double z1 = z0 + 1.0;

            // Full containment of sphere in this cell:
            const bool sphereFullyInsideSomeCell =
                (origin.x() - radius >= x0) && (origin.x() + radius <= x1) &&
                (origin.y() - radius >= y0) && (origin.y() + radius <= y1) &&
                (origin.z() - radius >= z0) && (origin.z() + radius <= z1);

            if (sphereFullyInsideSomeCell) {
                // zero everything
                for (int i = 0; i < stencil_size; ++i) {
                    for (int j = 0; j < stencil_size; ++j) {
                        for (int k = 0; k < stencil_size; ++k) {
                            vfrac[i][j][k] = 0.0;
                            firstMoment[i][j][k].setZero();
                        }
                    }
                }
                const double V = (4.0 / 3.0) * M_PI * radius * radius * radius; // cell volume is 1
                vfrac[ci][cj][ck] = V;

                // firstMoment stores raw first moment (∫x dV, ∫y dV, ∫z dV)
                firstMoment[ci][cj][ck] = V * origin;

                if (secondMoment != nullptr) {
                    // centered second moment tensor for a solid sphere about its centroid: (V r^2 / 5) * I
                    *secondMoment = (V * radius * radius / 5.0) * Eigen::Matrix3d::Identity();
                }

                if (visualize) {
                    std::cout << "Sphere fully inside coarse cell ("
                            << ci << "," << cj << "," << ck
                            << "), skipping refined intersection and using analytical values.\n";
                }
                return; // accept this sample
            }

            // Analytical is not possible, need to do refined intersection

            const int min_refined_cells_across_diameter = 5;

            // Choose refinement so that (2*radius) / (1/refinement_factor) >= min_refined_cells_across_diameter
            const int refinement_factor = std::max(
                1,
                static_cast<int>(std::ceil(static_cast<double>(min_refined_cells_across_diameter) / (2.0 * radius)))
            );
            const double refined_cell_size = 1.0 / static_cast<double>(refinement_factor);

            // Coarse coordinates
            for (int coarse_index = 0; coarse_index <= stencil_size; ++coarse_index) {
                coarse_coords[coarse_index] = coarse_stencil_min + static_cast<double>(coarse_index);
            }

            // Zero out coarse outputs before accumulation
            for (int coarse_i = 0; coarse_i < stencil_size; ++coarse_i) {
                for (int coarse_j = 0; coarse_j < stencil_size; ++coarse_j) {
                    for (int coarse_k = 0; coarse_k < stencil_size; ++coarse_k) {
                        vfrac[coarse_i][coarse_j][coarse_k] = 0.0;
                        firstMoment[coarse_i][coarse_j][coarse_k].setZero();
                        if (visualize) {
                            centroid[coarse_i][coarse_j][coarse_k].setZero();
                        }
                    }
                }
            }

            // Initialize IRL outputs
            using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;
            IRL::GeneralMoments3D<2> total_general_moments =
                        IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // exact raw 0th/1st/2nd about global origin


            // Sphere Axis aligned bounding box (AABB) in global coordinates
            const double sphere_min_x = origin.x() - radius;
            const double sphere_max_x = origin.x() + radius;
            const double sphere_min_y = origin.y() - radius;
            const double sphere_max_y = origin.y() + radius;
            const double sphere_min_z = origin.z() - radius;
            const double sphere_max_z = origin.z() + radius;

            // Coarse cell index range that intersects sphere AABB
            int coarse_i_start = clamp_coarse_index(static_cast<int>(std::floor(sphere_min_x - coarse_stencil_min)));
            int coarse_i_end   = clamp_coarse_index(static_cast<int>(std::floor(sphere_max_x - coarse_stencil_min)));
            int coarse_j_start = clamp_coarse_index(static_cast<int>(std::floor(sphere_min_y - coarse_stencil_min)));
            int coarse_j_end   = clamp_coarse_index(static_cast<int>(std::floor(sphere_max_y - coarse_stencil_min)));
            int coarse_k_start = clamp_coarse_index(static_cast<int>(std::floor(sphere_min_z - coarse_stencil_min)));
            int coarse_k_end   = clamp_coarse_index(static_cast<int>(std::floor(sphere_max_z - coarse_stencil_min)));

            auto clamp_refined_index = [&](int index_value) {
                return std::max(0, std::min(refinement_factor - 1, index_value));
            };

            for (int coarse_i = coarse_i_start; coarse_i <= coarse_i_end; ++coarse_i) {
                // Coarse cell bounds in x for this coarse_i
                const double coarse_x0 = coarse_stencil_min + static_cast<double>(coarse_i);
                const double coarse_x1 = coarse_x0 + 1.0;

                // Refined index range in x within this coarse cell whose subcells overlap the sphere AABB
                int refined_i_start = clamp_refined_index(
                    static_cast<int>(std::floor((sphere_min_x - coarse_x0) / refined_cell_size))
                );
                int refined_i_end = clamp_refined_index(
                    static_cast<int>(std::floor((sphere_max_x - coarse_x0) / refined_cell_size))
                );
                // If no refined cells in this coarse cell overlap the sphere AABB, skip to next coarse cell
                if (refined_i_end < refined_i_start) continue;

                for (int coarse_j = coarse_j_start; coarse_j <= coarse_j_end; ++coarse_j) {
                    const double coarse_y0 = coarse_stencil_min + static_cast<double>(coarse_j);
                    const double coarse_y1 = coarse_y0 + 1.0;

                    int refined_j_start = clamp_refined_index(
                        static_cast<int>(std::floor((sphere_min_y - coarse_y0) / refined_cell_size))
                    );
                    int refined_j_end = clamp_refined_index(
                        static_cast<int>(std::floor((sphere_max_y - coarse_y0) / refined_cell_size))
                    );
                    if (refined_j_end < refined_j_start) continue;

                    for (int coarse_k = coarse_k_start; coarse_k <= coarse_k_end; ++coarse_k) {
                        const double coarse_z0 = coarse_stencil_min + static_cast<double>(coarse_k);
                        const double coarse_z1 = coarse_z0 + 1.0;

                        int refined_k_start = clamp_refined_index(
                            static_cast<int>(std::floor((sphere_min_z - coarse_z0) / refined_cell_size))
                        );
                        int refined_k_end = clamp_refined_index(
                            static_cast<int>(std::floor((sphere_max_z - coarse_z0) / refined_cell_size))
                        );
                        if (refined_k_end < refined_k_start) continue;

                        // Refine only the relevant subcell index box inside this coarse cell
                        // from knowing refined_i_start/end, refined_j_start/end, refined_k_start/end
                        for (int refined_i_in_coarse = refined_i_start; refined_i_in_coarse <= refined_i_end; ++refined_i_in_coarse) {
                            // Refined cell bounds in x for this refined_i_in_coarse
                            const double refined_x0 = coarse_x0 + static_cast<double>(refined_i_in_coarse) * refined_cell_size;
                            const double refined_x1 = refined_x0 + refined_cell_size;

                            for (int refined_j_in_coarse = refined_j_start; refined_j_in_coarse <= refined_j_end; ++refined_j_in_coarse) {
                                // Refined cell bounds in y for this refined_j_in_coarse
                                const double refined_y0 = coarse_y0 + static_cast<double>(refined_j_in_coarse) * refined_cell_size;
                                const double refined_y1 = refined_y0 + refined_cell_size;

                                for (int refined_k_in_coarse = refined_k_start; refined_k_in_coarse <= refined_k_end; ++refined_k_in_coarse) {
                                    // Refined cell bounds in z for this refined_k_in_coarse
                                    const double refined_z0 = coarse_z0 + static_cast<double>(refined_k_in_coarse) * refined_cell_size;
                                    const double refined_z1 = refined_z0 + refined_cell_size;

                                    // Create refined subcell
                                    auto refined_cell = IRL::RectangularCuboid::fromBoundingPts(
                                        IRL::Pt(refined_x0, refined_y0, refined_z0),
                                        IRL::Pt(refined_x1, refined_y1, refined_z1)
                                    );

                                    // Use refined subcell center for paraboloid construction
                                    const Eigen::Vector3d refined_cell_center(
                                        0.5 * (refined_x0 + refined_x1),
                                        0.5 * (refined_y0 + refined_y1),
                                        0.5 * (refined_z0 + refined_z1)
                                    );

                                    const Eigen::Vector3d sphere_center_to_cell_center = refined_cell_center - origin;

                                    // Datum point on the sphere in the direction of the cell center
                                    const Eigen::Vector3d datum_paraboloid_eigen =
                                        origin + radius * sphere_center_to_cell_center.normalized();

                                    const IRL::Pt datum_paraboloid(
                                        datum_paraboloid_eigen.x(),
                                        datum_paraboloid_eigen.y(),
                                        datum_paraboloid_eigen.z()
                                    );

                                    // Local frame for paraboloid
                                    Eigen::Vector3d paraboloid_z_axis = sphere_center_to_cell_center;
                                    paraboloid_z_axis.normalize();

                                    const Eigen::Vector3d helper_vector(
                                        paraboloid_z_axis.x(),
                                        paraboloid_z_axis.y() + 1.0,
                                        paraboloid_z_axis.z()
                                    );

                                    Eigen::Vector3d paraboloid_x_axis = paraboloid_z_axis.cross(helper_vector);
                                    paraboloid_x_axis.normalize();

                                    Eigen::Vector3d paraboloid_y_axis = paraboloid_z_axis.cross(paraboloid_x_axis);
                                    paraboloid_y_axis.normalize();

                                    const auto reference_frame = IRL::ReferenceFrame(
                                        IRL::Normal(paraboloid_x_axis.x(), paraboloid_x_axis.y(), paraboloid_x_axis.z()),
                                        IRL::Normal(paraboloid_y_axis.x(), paraboloid_y_axis.y(), paraboloid_y_axis.z()),
                                        IRL::Normal(paraboloid_z_axis.x(), paraboloid_z_axis.y(), paraboloid_z_axis.z())
                                    );

                                    const auto paraboloid = IRL::Paraboloid(
                                        datum_paraboloid,
                                        reference_frame,
                                        1.0 / (2.0 * radius), // curvature in x direction
                                        1.0 / (2.0 * radius)  // curvature in y direction, same as x for sphere
                                    );

                                    // Intersect refined cell with paraboloid -> moments (+ surface for visualization)
                                    auto volume_and_surface = getVolumeMoments<VolumeMomentsAndSurface>(refined_cell, paraboloid);

                                    const double refined_volume = volume_and_surface.getMoments().volume();

                                    // IMPORTANT: this assumes getMoments().centroid() returns raw first moments (∫x dV, ∫y dV, ∫z dV)
                                    const Eigen::Vector3d refined_first_moment(
                                        volume_and_surface.getMoments().centroid().x(),
                                        volume_and_surface.getMoments().centroid().y(),
                                        volume_and_surface.getMoments().centroid().z()
                                    );

                                    // Accumulate directly into the parent coarse cell
                                    vfrac[coarse_i][coarse_j][coarse_k] += refined_volume;
                                    firstMoment[coarse_i][coarse_j][coarse_k] += refined_first_moment;

                                    // Accumulate exact raw 0th/1st/2nd moments if requested
                                    if (secondMoment != nullptr) {
                                        auto refined_general_moments =
                                            IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(refined_cell, paraboloid);
                                        total_general_moments += refined_general_moments;
                                    }

                                    if (visualize) {
                                        surfaces.push_back(volume_and_surface.getSurface());
                                    }
                                }
                            }
                        }
                    }
                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        // liquid-centered (about global liquid centroid) second moment matrix from refined accumulated moments
                        *secondMoment = centeredSecondMomentFromTotal(total_general_moments);
                    }
                }
            }
        }

        void generateSphere (
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, double min_radius = machineZero, double max_radius = 0.5, double radius_stddev = 0.0, bool visualize = false,
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

                    // Random radius
                    double radius = max_radius;
                    if (radius_stddev > 0.0) {
                        radius = sample_truncated_normal(0, radius_stddev, min_radius, max_radius);
                    }else{
                        std::uniform_real_distribution<double> random_radius(min_radius, max_radius);
                        radius = random_radius(eng);
                    }

                    Eigen::Vector3d origin = generateRandomPoint(-0.5 - radius , 0.5 + radius, eng);

                    // for visualization option:
                    std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;

                    std::vector<double> coarse_coords(stencil_size + 1);

                    generateSpecificSphere(
                        vfrac,
                        firstMoment,
                        centroid,
                        surfaces,
                        stencil_size,
                        origin,
                        radius,
                        coarse_coords,
                        visualize,
                        secondMoment
                    );

                    // Check central cell
                    int mid = stencil_size / 2;
                    double center_vfrac = vfrac[mid][mid][mid];
                    if (center_vfrac > machineZero && center_vfrac < 1.0-machineZero) {
                        if (visualize) {
                            WriteField(stencil_size, coarse_coords, vfrac, "vfrac");
                            WriteSurface(surfaces, "surface");
                            //printCentroids(centroid); this doesnt work anymore, could be adapted but no need
                        }
                        return; // accept this sample
                    }
                    // else: reject and regenerate
                    
                }
        }

        void generateBentTruncatedCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, bool truncateInsideCentralCell, double min_radius = machineZero, double max_radius = 0.5, double radius_stddev = 0.0, 
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

                // Random tube radius in [0, max_radius], with optional stddev for truncated normal sampling
                double tube_radius = max_radius;
                if (radius_stddev > 0.0) {
                    tube_radius = sample_truncated_normal(0, radius_stddev, min_radius, max_radius);
                } else {
                    std::uniform_real_distribution<double> random_thickness(min_radius, max_radius);
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
                /*
                {
                    double sphere_radius_min = tube_radius;
                    double sphere_radius_max = 2 * tube_radius;
                    sphere_radius = std::abs(sample_truncated_normal(sphere_radius_min, 0.4*sphere_radius_min, sphere_radius_min, sphere_radius_max));
                }
                */

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
                        const double thresh_outside = half_diag;
                        const double thresh_inside = 0.5;
                        const double dist_p_cut  = p_cut.norm(); // distance from origin to point on circle at theta_cut, which is the center of the spherical truncation cap
                        if (!truncateInsideCentralCell) {
                            if (std::abs(dist_p_cut - sphere_radius) < thresh_outside) continue; // want outside, but it's inside
                        } else {
                            if (std::abs(dist_p_cut - sphere_radius) > thresh_inside) continue; // want inside, but it's outside
                        }

                        theta_cut = t;
                        sphere_origin = p_cut;
                        found = true;
                        break;
                    }
                    if (!found) continue;
                }

                // range of theta for tube in [theta_start, theta_cut]
                double tube_theta_span = M_PI; // Half circle, with a high circle radius this should be outside of the stencil

                double theta_start = theta_cut - tube_theta_span;

                auto thetaInTubeRange = [&](double theta) -> bool {
                    // unwrap theta close to theta_cut so interval comparison works
                    double t = unwrapNear(theta, theta_cut);
                    return (t >= theta_start && t <= theta_cut);
                };

                // Refined mesh
                const double cell_volume = 1.0;
                double max_refinement_factor = 8.0;
                double refinement_factor_double = std::ceil(3.0/(2.0*tube_radius)); // want at least ~3 samples across the tube diameter for decent accuracy, can adjust this heuristic as needed
                int refinement_factor = static_cast<int>(refinement_factor_double);
                int refined_stencil_size = refinement_factor * stencil_size;
                // need to declare refined coords outside if below:
                // Refined cell coordinates
                auto coords = std::vector<double>(refined_stencil_size + 1);
                for (int i = 0; i <= refined_stencil_size; i++) {
                    coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
                }
                int mid = stencil_size / 2;
                // Also need to define total moments refined for refinement option
                IRL::GeneralMoments3D<2> totalMoments_refined =
                        IRL::GeneralMoments3D<2>::fromScalarConstant(0.0); // For exact 2nd moment

                // If the tube radius is too small, approximate the circle by lines
                // need the following parameters in this scope first:
                // If we want stencil-level 2nd moments (liquid-centered), accumulate raw moments across the whole stencil
                double totalV = 0.0;
                Eigen::Vector3d totalM1 = Eigen::Vector3d::Zero(); // raw first moment = ∑ V*c
                Eigen::Matrix3d totalM2 = Eigen::Matrix3d::Zero(); // raw second moment = ∑ ∫ x x^T dV

                if (refinement_factor_double > max_refinement_factor) {
                    // Analytical approximation with cylinders and a sphere
                    const double maxSegLen = 0.25;   // max line length in world units (cell size = 1)
                    const double cell_volume = 1.0;
                    const int mid = stencil_size / 2;

                    // Build coarse coordinates for the stencil (cell size = 1)
                    std::vector<double> coords_coarse(stencil_size + 1);
                    for (int ii = 0; ii <= stencil_size; ++ii) {
                        coords_coarse[ii] = -0.5 * static_cast<double>(stencil_size) + static_cast<double>(ii);
                    }

                    // Clear outputs
                    for (int i = 0; i < stencil_size; ++i)
                        for (int j = 0; j < stencil_size; ++j)
                            for (int k = 0; k < stencil_size; ++k) {
                                vfrac[i][j][k] = 0.0;
                                firstMoment[i][j][k].setZero();
                            }
                    
                    // Build sphere

                    // raw 0th/1st/2nd moments for the sphere contribution
                    IRL::GeneralMoments3D<2> sphere_total_general_moments =
                        IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);

                    std::vector<double> sphere_coarse_coords(stencil_size + 1); // Should be the same as coarse coords above

                    generateSpecificSphere(
                        vfrac,
                        firstMoment,
                        centroid,
                        surfaces,
                        stencil_size,
                        sphere_origin,
                        sphere_radius,
                        sphere_coarse_coords,
                        false,     // Lets not use visualize
                        nullptr   // do not compute centered 2nd moment here, this would require more effort to combine
                    );

                    // Clip segment p(t)=p0+t*(p1-p0), t in [0,1], to Cell AABB [bmin,bmax]; returns ta,tb range inside box.
                    auto clipSegmentToCell = [&](const Eigen::Vector3d& p0,
                                                const Eigen::Vector3d& p1,
                                                const Eigen::Vector3d& bmin,
                                                const Eigen::Vector3d& bmax,
                                                double& ta, double& tb) -> bool
                    {
                        ta = 0.0; tb = 1.0;
                        const Eigen::Vector3d d = p1 - p0;

                        for (int ax = 0; ax < 3; ++ax) {
                            const double p = d[ax];
                            const double q0 = p0[ax];

                            if (std::abs(p) < 1e-14) {
                                // Parallel to slab: must already be within bounds on this axis
                                if (q0 < bmin[ax] || q0 > bmax[ax]) return false;
                            } else {
                                const double invp = 1.0 / p;
                                double tNear = (bmin[ax] - q0) * invp;
                                double tFar  = (bmax[ax] - q0) * invp;
                                if (tNear > tFar) std::swap(tNear, tFar);
                                ta = std::max(ta, tNear);
                                tb = std::min(tb, tFar);
                                if (ta > tb) return false;
                            }
                        }
                        return true;
                    };

                    // Central (about centroid) 2nd-moment tensor for a solid cylinder of volume V, length L, axis unit a.
                    // C2 = V*(L^2/12) * (aa^T) + V*(r^2/4) * (I - aa^T)
                    auto cylinderCentralSecondMoment = [&](double V, double L, const Eigen::Vector3d& a_unit) -> Eigen::Matrix3d {
                        const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
                        const Eigen::Matrix3d aaT = a_unit * a_unit.transpose();
                        const Eigen::Matrix3d P = I - aaT;
                        const double L2 = L * L;
                        const double r2 = tube_radius * tube_radius;
                        return (V * (L2 / 12.0)) * aaT + (V * (r2 / 4.0)) * P;
                    };

                    // Centerline point for a given theta on the circle in the (u,v) plane
                    auto circlePoint = [&](double theta) -> Eigen::Vector3d {
                        // You already have centerlinePoint(theta, c0, radius_circle, u, v)
                        return centerlinePoint(theta, c0, radius_circle, u, v);
                    };

                    // --- Discretize ONLY the truncated tube arc [theta_start, theta_cut] (ignore end sphere entirely) ---

                    // Make sure comparisons happen on an unwrapped interval near theta_cut
                    const double theta_start_u = unwrapNear(theta_start, theta_cut);
                    const double theta_cut_u   = theta_cut; // unwrapNear(theta_cut,theta_cut) would be itself

                    const double span = std::max(0.0, theta_cut_u - theta_start_u);            // radians
                    const double arcLen = radius_circle * span;                                 // length along centerline
                    const int nSeg = std::max(1, static_cast<int>(std::ceil(arcLen / maxSegLen)));

                    // --- Iterate coarse cells and accumulate clipped cylinder pieces from polyline segments ---
                    for (int i = 0; i < stencil_size; ++i) {
                        for (int j = 0; j < stencil_size; ++j) {
                            for (int k = 0; k < stencil_size; ++k) {

                                const Eigen::Vector3d bmin(coords_coarse[i],   coords_coarse[j],   coords_coarse[k]);
                                const Eigen::Vector3d bmax(coords_coarse[i+1], coords_coarse[j+1], coords_coarse[k+1]);

                                double cellV = 0.0;
                                Eigen::Vector3d cellM1 = Eigen::Vector3d::Zero(); // ∫ x dV
                                Eigen::Matrix3d cellM2 = Eigen::Matrix3d::Zero(); // ∫ x x^T dV

                                for (int s = 0; s < nSeg; ++s) {
                                    // theta along the truncated arc, linearly in parameter
                                    const double tA = static_cast<double>(s)     / static_cast<double>(nSeg);
                                    const double tB = static_cast<double>(s + 1) / static_cast<double>(nSeg);

                                    const double thA = theta_start_u + tA * span;
                                    const double thB = theta_start_u + tB * span;

                                    const Eigen::Vector3d p0 = circlePoint(thA);
                                    const Eigen::Vector3d p1 = circlePoint(thB);

                                    // Does this centerline segment intersect this cell AABB?
                                    double ta, tb;
                                    if (!clipSegmentToCell(p0, p1, bmin, bmax, ta, tb)) continue;

                                    // Keep only the part inside the cell
                                    const Eigen::Vector3d q0 = p0 + ta * (p1 - p0);
                                    const Eigen::Vector3d q1 = p0 + tb * (p1 - p0);

                                    const Eigen::Vector3d dseg = q1 - q0;
                                    const double L = dseg.norm(); 
                                    if (L < 1e-14) continue;

                                    const Eigen::Vector3d axis = dseg / L;

                                    // Solid cylinder approximation for this clipped piece (radius = tube_radius)
                                    const double V = M_PI * tube_radius * tube_radius * L;
                                    const Eigen::Vector3d c = 0.5 * (q0 + q1);

                                    // Central 2nd moment about centroid -> raw about origin
                                    const Eigen::Matrix3d C2   = cylinderCentralSecondMoment(V, L, axis);
                                    const Eigen::Matrix3d raw2 = C2 + V * (c * c.transpose());

                                    cellV  += V;
                                    cellM1 += V * c;
                                    cellM2 += raw2;
                                }

                                if (cellV > 0.0) {
                                    double addV = cellV / cell_volume;

                                    // limit so final volume fraction does not exceed 1
                                    addV = std::min(addV, 1.0 - vfrac[i][j][k]);

                                    if (addV > 0.0) {
                                        double scale = addV / (cellV / cell_volume);  // fraction of the cylinder volume we actually keep

                                        vfrac[i][j][k] += addV;
                                        firstMoment[i][j][k] += scale * cellM1;
                                    }

                                    if (secondMoment != nullptr) {
                                        // Doenst take sphere into account
                                        totalV  += cellV;
                                        totalM1 += cellM1;
                                        totalM2 += cellM2;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    // Use refined definition
                    std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                        std::vector<std::vector<double>>(refined_stencil_size,
                        std::vector<double>(refined_stencil_size)));

                    std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                        std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                        std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

                    using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

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

                                // Compute closest point and tangent on the full circle (not truncated) for this cell center
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
                                double d2 = dist_to_circle_centerline - tube_radius; // tube
                                double d1 = (cell_center - sphere_origin).norm() - sphere_radius; // sphere

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
                                    // Sphere interface
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
                                    // Tube interface (local cylinder)
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
                }

                // check central cell
                double center_vfrac = vfrac[mid][mid][mid];
                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {

                    // Now calc stencil 2nd moments if requested
                    if (secondMoment != nullptr) {
                        if (refinement_factor_double > max_refinement_factor) {
                            std::cout<<"Warning: second moment doesnt work for underresolved truncated bent cylinder"<<std::endl;
                            // Warning: second moment doesnt work for underresolved truncated bent cylinder
                            const Eigen::Vector3d xc = totalM1 / totalV;
                            const Eigen::Matrix3d central = totalM2 - totalV * (xc * xc.transpose());
                            *secondMoment = central; // or: central / totalV
                        } else {
                            // liquid-centered (about global liquid centroid) second moment matrix from refined accumulated moments
                            *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                        }
                    }

                    if (visualize) {
                        WriteField(stencil_size, coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                        printCentroids(centroid);
                    }
                    return; // accept this sample
                }

                // else: reject and regenerate
            }
        }

        int sampleWeightedIndex(const std::vector<double>& weights) {
            double sum = 0.0;
            for (double w : weights) {
                if (w < 0.0) {
                    throw std::runtime_error("Negative weight in sampleWeightedIndex.");
                }
                sum += w;
            }

            if (sum <= 0.0) {
                throw std::runtime_error("All weights are zero in sampleWeightedIndex.");
            }

            std::uniform_real_distribution<double> dist(0.0, sum);
            double r = dist(eng);

            double cumulative = 0.0;
            for (int i = 0; i < static_cast<int>(weights.size()); ++i) {
                cumulative += weights[i];
                if (r <= cumulative) {
                    return i;
                }
            }

            return static_cast<int>(weights.size()) - 1;
        }


        std::vector<float> generateState(int datapoint_type, int stencil_size, int include_Moments = 0, bool include_Eigenvalues = false,
                                        double paraboloid_coeff_stddev = 0.1,
                                        double sheet_coeff_stddev = 0.1, double max_sheet_thickness = 0.5, double sheet_thickness_stddev = 0.0,
                                        double max_cylinder_radius = 0.5, double cylinder_radius_stddev = 0.0, 
                                        double max_sphere_radius = 0.5, double sphere_radius_stddev = 0.0,
                                        bool visualize = false, bool exact_2nd_moments = false,
                                        double class0_max_characteristic = 2.5)
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

            if (exact_2nd_moments && include_Moments >= 2) {
                std::cout<<"Enable exact 2nd moment calculation"<<std::endl;
                secondMomentPtr = &secondMoment; // This enables calculation of exact 2nd moments for all switch cases below
            }

            const std::array<double, 6>& class0_subclass_weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; // Adjust these weights to control the relative frequency of each subclass within class 0

            std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

            // Testing
            double radius_circ_min = 3.0;
            double radius_circ_max = 10.0;

            double min_cylinder_radius = machineZero;
            double min_sheet_thickness = machineZero;
            double min_sphere_radius = machineZero;
            

            //std::cout<<"Generating datapoint of type "<<datapoint_type<<std::endl;
            switch (datapoint_type) {
                case 1: {
                    double p1 = prob_dist(eng);  // draw a random number in [0,1)

                    if (p1 < 0.2) {
                        // 20% chance → truncated cylinder
                        //generateTruncatedCylinder(vfrac, firstMoment, stencil_size, max_cylinder_radius, cylinder_radius_stddev, visualize, secondMomentPtr);
                        generateBentTruncatedCylinder(vfrac, firstMoment, stencil_size, /*truncateinsidecentralcell*/false, min_cylinder_radius, max_cylinder_radius, cylinder_radius_stddev, radius_circ_min, radius_circ_max, visualize, secondMomentPtr);
                    } else {
                        // 80% chance → normal cylinder
                        //generateCylinder(vfrac, firstMoment, stencil_size, max_cylinder_radius, cylinder_radius_stddev, visualize, secondMomentPtr);
                        generateBentCylinder(vfrac, firstMoment, stencil_size, min_cylinder_radius, max_cylinder_radius, cylinder_radius_stddev, radius_circ_min, radius_circ_max, visualize, secondMomentPtr);
                    }

                    break;
                }
                case 2:
                    generateSphere(vfrac, firstMoment, stencil_size, min_sphere_radius, max_sphere_radius, sphere_radius_stddev, visualize, secondMomentPtr);
                    break;
                case 3: {
                    double p2 = prob_dist(eng);  // draw a random number in [0,1)
                    //generateSheet(vfrac, firstMoment, stencil_size, sheet_coeff_stddev, min_sheet_thickness, max_sheet_thickness, sheet_thickness_stddev, visualize, secondMomentPtr);
                    if (p2 < 0.2) {
                        // 20% chance → cut sheet
                        generateCutSheet(vfrac, firstMoment, stencil_size, /*cutinsidecentralcell*/ false,sheet_coeff_stddev, min_sheet_thickness, max_sheet_thickness, sheet_thickness_stddev, visualize, secondMomentPtr);
                    } else {
                        // 80% chance → normal sheet
                        generateSheet(vfrac, firstMoment, stencil_size, sheet_coeff_stddev, min_sheet_thickness,max_sheet_thickness, sheet_thickness_stddev, visualize, secondMomentPtr);
                    }
                    break;
                }
                case 4:
                    generateBentTruncatedCylinder(vfrac, firstMoment, stencil_size, /*truncateinsidecentralcell*/true, min_cylinder_radius, max_cylinder_radius, cylinder_radius_stddev, radius_circ_min, radius_circ_max, visualize, secondMomentPtr);
                    //generateBentCylinder(vfrac, firstMoment, stencil_size, max_cylinder_radius, cylinder_radius_stddev, visualize, secondMomentPtr);
                    //generateBentTruncatedCylinder(vfrac, firstMoment, stencil_size, truncateInsideCentrCell, max_cylinder_radius, cylinder_radius_stddev, radius_circ_min, radius_circ_max, visualize, secondMomentPtr);
                    //generateCutSheet(vfrac, firstMoment, stencil_size, sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev, visualize, secondMomentPtr);
                    //generateSphere(vfrac, firstMoment, stencil_size, max_sphere_radius, sphere_radius_stddev, visualize, secondMomentPtr);
                    break;
                case 5:
                    generateCutSheet(vfrac, firstMoment, stencil_size, /*cutinsidecentralcell*/ true, sheet_coeff_stddev, min_sheet_thickness, max_sheet_thickness, sheet_thickness_stddev, visualize, secondMomentPtr);
                    break;
                default:
                    std::vector<double> weights(
                        class0_subclass_weights.begin(),
                        class0_subclass_weights.end()
                    );
                    int subtype0 = sampleWeightedIndex(weights);
                    //std::cout<<"Subtype: "<<subtype0<<std::endl;

                    double resolved_min_cylinder_radius = max_cylinder_radius;
                    double resolved_max_cylinder_radius = class0_max_characteristic;

                    double resolved_min_sphere_radius = max_sphere_radius;
                    double resolved_max_sphere_radius = class0_max_characteristic;

                    double resolved_min_sheet_thickness = max_sheet_thickness;
                    double resolved_max_sheet_thickness = class0_max_characteristic;

                    switch (subtype0) {
                        case 0:
                            generateParaboloid(vfrac, firstMoment, stencil_size,
                                            paraboloid_coeff_stddev, visualize, secondMomentPtr);
                            break;

                        case 1: {
                            // Well-resolved ligament-like object
                            double p3 = prob_dist(eng);  // draw a random number in [0,1)

                            if (p3 < 0.2) {
                                generateBentTruncatedCylinder(vfrac, firstMoment, stencil_size,
                                                        /*truncateinsidecentralcell*/ false,
                                                        resolved_min_cylinder_radius,
                                                        resolved_max_cylinder_radius,
                                                        cylinder_radius_stddev,
                                                        radius_circ_min, radius_circ_max,
                                                        visualize, secondMomentPtr);
                            } else {
                                generateBentCylinder(vfrac, firstMoment, stencil_size,
                                                    resolved_min_cylinder_radius,
                                                    resolved_max_cylinder_radius,
                                                    cylinder_radius_stddev,
                                                    radius_circ_min, radius_circ_max,
                                                    visualize, secondMomentPtr);
                            }
                            break;
                        }

                        case 2:
                            // Well-resolved drop-like object
                            
                            generateSphere(vfrac, firstMoment, stencil_size,
                                        resolved_min_sphere_radius,
                                        resolved_max_sphere_radius,
                                        sphere_radius_stddev,
                                        visualize, secondMomentPtr);
                            
                            /*
                            generateParaboloid(vfrac, firstMoment, stencil_size,
                                            paraboloid_coeff_stddev, visualize, secondMomentPtr);
                            */
                            break;

                        case 3: {
                            // Well-resolved sheet-like object
                            double p4 = prob_dist(eng);  // draw a random number in [0,1)
                            
                            if (p4 < 0.2) {
                                generateCutSheet(vfrac, firstMoment, stencil_size,
                                            /*cutinsidecentralcell*/ false,
                                            sheet_coeff_stddev,
                                            resolved_min_sheet_thickness,
                                            resolved_max_sheet_thickness,
                                            sheet_thickness_stddev,
                                            visualize, secondMomentPtr);
                            } else {
                                generateSheet(vfrac, firstMoment, stencil_size,
                                            sheet_coeff_stddev,
                                            resolved_min_sheet_thickness,
                                            resolved_max_sheet_thickness,
                                            sheet_thickness_stddev,
                                            visualize, secondMomentPtr);
                            }
                            break;
                        }
                        case 4:
                            // Well-resolved ligament-tip-like object
                            generateBentTruncatedCylinder(vfrac, firstMoment, stencil_size,
                                                        /*truncateinsidecentralcell*/ true,
                                                        resolved_min_cylinder_radius,
                                                        resolved_max_cylinder_radius,
                                                        cylinder_radius_stddev,
                                                        radius_circ_min, radius_circ_max,
                                                        visualize, secondMomentPtr);
                            break;

                        case 5:
                            // Well-resolved sheet-edge-like object
                            generateCutSheet(vfrac, firstMoment, stencil_size,
                                            /*cutinsidecentralcell*/ true,
                                            sheet_coeff_stddev,
                                            resolved_min_sheet_thickness,
                                            resolved_max_sheet_thickness,
                                            sheet_thickness_stddev,
                                            visualize, secondMomentPtr);
                            break;

                        default:
                            generateParaboloid(vfrac, firstMoment, stencil_size,
                                            paraboloid_coeff_stddev, visualize, secondMomentPtr);
                            break;
                    }

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

            if (include_Eigenvalues) {
                Eigen::Matrix3d I = IRL::computeInertiaTensor(flattened_state, stencil_size, include_Moments, machineZero);

                // Get eigenvalues
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(I);
                Eigen::Vector3d evals = solver.eigenvalues();

                // Sort eigenvalues descending: I1 >= I2 >= I3
                std::sort(evals.data(), evals.data() + 3, std::greater<double>());
                double I1 = evals[0], I2 = evals[1], I3 = evals[2];

                if (include_Moments <= 1) {
                    flattened_state.push_back(I1);
                    flattened_state.push_back(I2);
                    flattened_state.push_back(I3);
                }

                if (include_Moments == 2 || include_Moments == 3) {
                    std::cout<<"WIP: include_Moments == 2 or 3 and include_Eigenvalues == true is not fully implemented yet"<<std::endl;
                }

                if (include_Moments == 4) {
                    // if include_Moments == 4, use only the three eigenvalues
                    flattened_state.clear();
                    flattened_state.push_back(I1);
                    flattened_state.push_back(I2);
                    flattened_state.push_back(I3);
                }

                
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

        void generateData (std::vector<std::vector<float>>* statesV, std::vector<int>* labelsV, int no_datapoints, int stencil_size = 3, int no_datapoint_types_in = 4, 
                                        int include_Moments = 0, bool include_Eigenvalues = false,
                                        double paraboloid_coeff_stddev = 0.1,
                                        double sheet_coeff_stddev = 0.1, double max_sheet_thickness = 0.5, double sheet_thickness_stddev = 0.0,
                                        double max_cylinder_radius = 0.5, double cylinder_radius_stddev = 0.0,
                                        double max_sphere_radius = 0.5, double sphere_radius_stddev = 0.0)
        
        {
            double class0_max_characteristic = 2.5;
            double class0_weight = 2.0;

            std::cout << no_datapoints << std::endl;
            // Initialize random number generator with a seed from current time
            std::srand(std::time(0));

            std::vector<double> class_weights(no_datapoint_types_in, 1.0);
            if (no_datapoint_types_in > 0) {
                class_weights[0] = class0_weight;
            }

            for (int i = 0; i < no_datapoints; i++) {
                int datapoint_type = sampleWeightedIndex(class_weights);

                labelsV->push_back(datapoint_type);
                statesV->push_back(generateState(
                    datapoint_type,
                    stencil_size,
                    include_Moments,
                    include_Eigenvalues,
                    paraboloid_coeff_stddev,
                    sheet_coeff_stddev,
                    max_sheet_thickness,
                    sheet_thickness_stddev,
                    max_cylinder_radius,
                    cylinder_radius_stddev,
                    max_sphere_radius,
                    sphere_radius_stddev,
                    false,   // visualize
                    false,   // exact_2nd_moments
                    class0_max_characteristic
                ));
                if (i % 10000 == 0) {
                    std::cout << "Generated " << i << " datapoints." << std::endl;
                }
            }
        }
    };
}