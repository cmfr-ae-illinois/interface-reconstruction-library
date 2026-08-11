#pragma once

#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <array>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <vector>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/ml_classification/vtk_out.h"
#include "irl/ml_classification/inertia_calc.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

namespace IRL {
    class Data_gen {
        protected:
        int stencil_size = 5;
        int no_datapoints = 4096*64;
        int include_Moments = 1;
        bool include_Surface_Area = false;
        bool include_Eigenvalues = false;
        double paraboloid_coeff_stddev = 0.1;
        double hyperbolic_cylinder_opening_angle_stddev = 20; //degrees
        double sheet_coeff_stddev = 0.1;
        double sheet_thickness_stddev = 0.0;
        double cylinder_radius_stddev = 0.0;
        double radius_circle_min = 2.5;
        double radius_circle_max = 10.0;
        double sphere_radius_stddev = 0.0;
        double ellipsoid_subgrid_stddev = 0.7;
        double min_long_ellipsoid_axis = 3.0;
        double max_long_ellipsoid_axis = 5.0;
        bool exact_2nd_moment = false;  // enable calculation of exact 2nd moments for data generation
        bool visualize = false; // if true, print centroids and / or write surfaces
        double machineZero = 1e-12;
        double lower_limit_subgrid = machineZero;
        double upper_limit_subgrid = std::sqrt(3.0);
        double class0_max_characteristic = 2.5;

        public:
        std::mt19937_64 eng;
        //int stencil_size = 3;
        //std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));


        Data_gen(){
            // Initialize random number generator with a seed from current time
            std::random_device rd;
            eng = std::mt19937_64(rd());  // Use random_device for seeding
            std::cout << "I'm a data generator!" << std::endl;
        }

        void set_stencil_size (int stencil_size_new) {
            stencil_size = stencil_size_new;
        }

        void updateDataParameters(int no_datapoints_in, int include_Moments_in, bool include_Surface_Area_in, bool include_Eigenvalues_in,
                                double paraboloid_coeff_stddev_in, double hyperbolic_cylinder_opening_angle_stddev_in,
                                double sheet_coeff_stddev_in, double sheet_thickness_stddev_in,
                                double cylinder_radius_stddev_in, double radius_circle_min_in, double radius_circle_max_in,
                                double sphere_radius_stddev_in, 
                                double ellipsoid_subgrid_stddev_in, double min_long_ellipsoid_axis_in, double max_long_ellipsoid_axis_in,
                                bool exact_2nd_mom = false, bool visualize_in = false, double machineZero_in = 1e-12, 
                                double lower_limit_subgrid_in = 1e-12, double upper_limit_subgrid_in = 1.732, double class0_max_characteristic_in = 2.5) {
            no_datapoints = no_datapoints_in;
            include_Moments = include_Moments_in;
            include_Surface_Area = include_Surface_Area_in;
            include_Eigenvalues = include_Eigenvalues_in;
            paraboloid_coeff_stddev = paraboloid_coeff_stddev_in;
            hyperbolic_cylinder_opening_angle_stddev = hyperbolic_cylinder_opening_angle_stddev_in;
            sheet_coeff_stddev = sheet_coeff_stddev_in;
            sheet_thickness_stddev = sheet_thickness_stddev_in;
            cylinder_radius_stddev = cylinder_radius_stddev_in;
            radius_circle_min = radius_circle_min_in;
            radius_circle_max = radius_circle_max_in;
            sphere_radius_stddev = sphere_radius_stddev_in;
            ellipsoid_subgrid_stddev = ellipsoid_subgrid_stddev_in;
            min_long_ellipsoid_axis = min_long_ellipsoid_axis_in;
            max_long_ellipsoid_axis = max_long_ellipsoid_axis_in;
            exact_2nd_moment = exact_2nd_mom;
            visualize = visualize_in;
            machineZero = machineZero_in;
            lower_limit_subgrid = lower_limit_subgrid_in;
            upper_limit_subgrid = upper_limit_subgrid_in;
            class0_max_characteristic = class0_max_characteristic_in;
        }

        void setVisualize (bool viz = true) {
            visualize = viz;
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
            std::vector<std::vector<std::vector<double>>>& area,
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
                Eigen::Vector3d datumVec = generateRandomPoint(-2.5, 2.5, eng);

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

                std::normal_distribution<double> random_coeff(0.0, paraboloid_coeff_stddev);
                //double coeff1 = random_coeff(eng);
                //double coeff2 = random_coeff(eng);
                double coeff1 = sample_truncated_normal(0.0, paraboloid_coeff_stddev, -1.0, 1.0);
                double coeff2 = sample_truncated_normal(0.0, paraboloid_coeff_stddev, -1.0, 1.0);

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

                            auto surface = getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid).getSurface();
                            area[i][j][k] = surface.getSurfaceArea();

                            if (visualize) {
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
        
        void generateHyperbolicCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            double min_thickness, double max_thickness, double max_opening_angle = 45,
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
                // TODO make similar to sheet
                // Defining cell coordinates
                auto coords = makeCenteredCoords(stencil_size);
                const double cell_volume = 1.0;

                // Random paraboloid parameters
                Eigen::Vector3d datumVec = generateRandomPoint(-0.5, 0.5, eng);

                Eigen::Vector3d direction = generateRandomDirection(eng);

                // Build orthonormal frame aligned with direction
                Eigen::Vector3d helper = generateRandomDirection(eng); // safe random helper
                Eigen::Vector3d hyperbolic_cylinder_x = direction.cross(helper);
                if (hyperbolic_cylinder_x.norm() < 1e-12) {
                    // fallback: pick fixed axis if helper was parallel
                    hyperbolic_cylinder_x = direction.cross(Eigen::Vector3d(1,0,0));
                }
                hyperbolic_cylinder_x.normalize();
                Eigen::Vector3d hyperbolic_cylinder_y = direction.cross(hyperbolic_cylinder_x);
                hyperbolic_cylinder_y.normalize();
                Eigen::Vector3d hyperbolic_cylinder_z = direction;

                const auto frame = IRL::ReferenceFrame(
                    IRL::Normal(hyperbolic_cylinder_x.x(), hyperbolic_cylinder_x.y(), hyperbolic_cylinder_x.z()), 
                    IRL::Normal(hyperbolic_cylinder_y.x(), hyperbolic_cylinder_y.y(), hyperbolic_cylinder_y.z()), 
                    IRL::Normal(hyperbolic_cylinder_z.x(), hyperbolic_cylinder_z.y(), hyperbolic_cylinder_z.z()));

                IRL::Pt datum(datumVec.x(), datumVec.y(), datumVec.z());

                std::uniform_real_distribution<double> random_thickness(min_thickness, max_thickness);
                double thickness = random_thickness(eng);
                double opening_angle = 0.0;
                if (hyperbolic_cylinder_opening_angle_stddev > 0.0) {
                    opening_angle = sample_truncated_normal(0.0, hyperbolic_cylinder_opening_angle_stddev, 0.0, max_opening_angle);
                } else {
                    //Use uniform distribution
                    std::uniform_real_distribution<double> random_angle(0.0, max_opening_angle);
                    opening_angle = random_angle(eng);
                }

                double coeff_r = 0.25 * thickness * thickness;
                double coeff_b = -std::pow(std::tan(0.5 * opening_angle), 2);

                const auto hyperbolic_cylinder = Cylinder(datum, frame, coeff_b, coeff_r);
                // Initialize field
                std::vector<CylinderParametrizedSurfaceOutput> surfaces;
                using VolumeMomentsAndSurface = AddSurfaceOutput<VolumeMoments, CylinderParametrizedSurfaceOutput>;
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
                                    cell, hyperbolic_cylinder, cell_volume,
                                    vfrac[i][j][k],
                                    firstMoment[i][j][k],
                                    totalMoments
                                );

                            } else {
                                auto volume_and_surface = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, hyperbolic_cylinder);
                                vfrac[i][j][k] = volume_and_surface.getMoments().volume() / cell_volume;
                                firstMoment[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                        volume_and_surface.getMoments().centroid().y(),
                                                        volume_and_surface.getMoments().centroid().z();
                            }

                            auto surface = getVolumeMoments<VolumeMomentsAndSurface>(cell, hyperbolic_cylinder).getSurface();
                            area[i][j][k] = surface.getSurfaceArea();

                            if (visualize) {
                                //surfaces.push_back(volume_and_surface.getSurface());
                                surfaces.push_back(surface);

                                centroid[i][j][k] = computeCentroidFromFirstMoment(
                                    firstMoment[i][j][k], vfrac[i][j][k] * cell_volume);
                            }
                        }
                    }
                }

                // Check central cell
                if (/*centerCellIsCut(vfrac, stencil_size, machineZero)*/true) {
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
                std::cout << "Rejecting hyperbolic cylinder sample, center cell not cut." << std::endl;

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
            std::vector<std::vector<std::vector<double>>>& area,
            //int stencil_size, double coeff_stddev = 0.1, double min_thickness = machineZero, double max_thickness = 0.5, double thickness_stddev = 0.0, 
            double min_thickness,
            double max_thickness,
            bool subgrid = true,
            bool variable_thickness = true, 
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

                Eigen::Vector3d datum = generateRandomPoint(-1.0,1.0,eng);

                // Random unbiased direction
                Eigen::Vector3d direction = generateRandomDirection(eng);

                // Random sheet thickness
                double thickness = max_thickness;
                if (sheet_thickness_stddev > 0.0) {
                    thickness = sample_truncated_normal(0, sheet_thickness_stddev, min_thickness, max_thickness);
                }else{
                    std::uniform_real_distribution<double> random_thickness(min_thickness, max_thickness);
                    thickness = random_thickness(eng);
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

                double coeff1 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);
                double coeff2 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);
                double coeff1p1 = coeff1;
                double coeff1p2 = coeff1;
                double coeff2p1 = coeff2;
                double coeff2p2 = coeff2;
                // Random coefficients
                if (variable_thickness) {
                    // Allow variable thickness sheet to have an even smaller thickness outside of central cell
                    if (sheet_thickness_stddev > 0.0) {
                        thickness = sample_truncated_normal(0, sheet_thickness_stddev, machineZero, max_thickness);
                    }else{
                        std::uniform_real_distribution<double> random_thickness(machineZero, max_thickness);
                        thickness = random_thickness(eng);
                    }

                    // If variable thickness, make sure coeffs are correlated so that sheet doesn't self-intersect
                    // add uniformly distributed noise
                    //std::uniform_real_distribution<double> noise_dist(0.0, thickness/4.0);

                    coeff1p1 += sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                    coeff1p2 -= sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);

                    coeff2p1 += sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                    coeff2p2 -= sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                

                    // Check if central cell is thicker or thinner than the well-resolved threshold
                
                    const int c = stencil_size / 2;

                    const double x0 = coords[c];
                    const double x1 = coords[c + 1];
                    const double y0 = coords[c];
                    const double y1 = coords[c + 1];
                    const double z0 = coords[c];
                    const double z1 = coords[c + 1];

                    bool reject_sample = false;

                    for (double X : {x0, x1}) {
                        for (double Y : {y0, y1}) {
                            for (double Z : {z0, z1}) {
                                const Eigen::Vector3d corner(X, Y, Z);
                                const Eigen::Vector3d rel = corner - datum;

                                const double sheet_x = rel.dot(paraboloid_x);
                                const double sheet_y = rel.dot(paraboloid_y);

                                const double local_thickness =
                                    thickness
                                    + (coeff1p1 - coeff1p2) * sheet_x * sheet_x
                                    + (coeff2p1 - coeff2p2) * sheet_y * sheet_y;

                                if ((subgrid && local_thickness > max_thickness) ||
                                    (!subgrid && local_thickness < min_thickness)) {
                                    reject_sample = true;
                                    break;
                                }
                            }
                            if (reject_sample) break;
                        }
                        if (reject_sample) break;
                    }

                    if (reject_sample) {
                        continue;
                    }
                }

                 // Two paraboloid datums offset along the direction
                Eigen::Vector3d datum_paraboloid1_eVec = datum - direction.normalized() * (thickness/2.0);
                Eigen::Vector3d datum_paraboloid2_eVec = datum + direction.normalized() * (thickness/2.0);


                IRL::Pt datum_paraboloid1(datum_paraboloid1_eVec.x(), datum_paraboloid1_eVec.y(), datum_paraboloid1_eVec.z());
                IRL::Pt datum_paraboloid2(datum_paraboloid2_eVec.x(), datum_paraboloid2_eVec.y(), datum_paraboloid2_eVec.z());

                auto paraboloid1 = Paraboloid(datum_paraboloid1, frame, coeff1p1, coeff2p1);
                auto paraboloid2 = Paraboloid(datum_paraboloid2, frame, coeff1p2, coeff2p2);

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
                            
                            auto surface1 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid1).getSurface();
                            auto surface2 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid2).getSurface();

                            area[i][j][k] = surface1.getSurfaceArea() + surface2.getSurfaceArea();

                            if (visualize) {
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
            std::vector<std::vector<std::vector<double>>>& area,
            bool cutInsideCentralCell,
            double min_thickness = 1e-12, 
            double max_thickness = 0.5, 
            bool subgrid = true,
            bool variable_thickness = false,
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

                Eigen::Vector3d datum = generateRandomPoint(-1.0,1.0,eng);

                // Random unbiased direction
                Eigen::Vector3d direction = generateRandomDirection(eng);

                // Random sheet thickness
                double thickness = max_thickness;
                if (sheet_thickness_stddev > 0.0) {
                    thickness = sample_truncated_normal(0, sheet_thickness_stddev, min_thickness, max_thickness);
                }else{
                    std::uniform_real_distribution<double> random_thickness(min_thickness, max_thickness);
                    thickness = random_thickness(eng);
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

                double coeff1 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);
                double coeff2 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);
                double coeff1p1 = coeff1;
                double coeff1p2 = coeff1;
                double coeff2p1 = coeff2;
                double coeff2p2 = coeff2;
                // Random coefficients
                if (variable_thickness) {
                    // Allow variable thickness sheet to have an even smaller thickness outside of central cell
                    if (sheet_thickness_stddev > 0.0) {
                        thickness = sample_truncated_normal(0, sheet_thickness_stddev, machineZero, max_thickness);
                    }else{
                        std::uniform_real_distribution<double> random_thickness(machineZero, max_thickness);
                        thickness = random_thickness(eng);
                    }

                    // If variable thickness, make sure coeffs are correlated so that sheet doesn't self-intersect
                    // add uniformly distributed noise
                    //std::uniform_real_distribution<double> noise_dist(0.0, thickness/4.0);

                    coeff1p1 += sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                    coeff1p2 -= sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);

                    coeff2p1 += sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                    coeff2p2 -= sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                

                    // Check if central cell is thicker or thinner than the well-resolved threshold
                
                    const int c = stencil_size / 2;

                    const double x0 = coords[c];
                    const double x1 = coords[c + 1];
                    const double y0 = coords[c];
                    const double y1 = coords[c + 1];
                    const double z0 = coords[c];
                    const double z1 = coords[c + 1];

                    bool reject_sample = false;

                    for (double X : {x0, x1}) {
                        for (double Y : {y0, y1}) {
                            for (double Z : {z0, z1}) {
                                const Eigen::Vector3d corner(X, Y, Z);
                                const Eigen::Vector3d rel = corner - datum;

                                const double sheet_x = rel.dot(paraboloid_x);
                                const double sheet_y = rel.dot(paraboloid_y);

                                const double local_thickness =
                                    thickness
                                    + (coeff1p1 - coeff1p2) * sheet_x * sheet_x
                                    + (coeff2p1 - coeff2p2) * sheet_y * sheet_y;

                                if ((subgrid && local_thickness > max_thickness) ||
                                    (!subgrid && local_thickness < min_thickness)) {
                                    reject_sample = true;
                                    break;
                                }
                            }
                            if (reject_sample) break;
                        }
                        if (reject_sample) break;
                    }

                    if (reject_sample) {
                        continue;
                    }
                }

                 // Two paraboloid datums offset along the direction
                Eigen::Vector3d datum_paraboloid1_eVec = datum - direction.normalized() * (thickness/2.0);
                Eigen::Vector3d datum_paraboloid2_eVec = datum + direction.normalized() * (thickness/2.0);


                IRL::Pt datum_paraboloid1(datum_paraboloid1_eVec.x(), datum_paraboloid1_eVec.y(), datum_paraboloid1_eVec.z());
                IRL::Pt datum_paraboloid2(datum_paraboloid2_eVec.x(), datum_paraboloid2_eVec.y(), datum_paraboloid2_eVec.z());

                auto paraboloid1 = Paraboloid(datum_paraboloid1, frame, coeff1p1, coeff2p1);
                auto paraboloid2 = Paraboloid(datum_paraboloid2, frame, coeff1p2, coeff2p2);

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

                            auto surface1 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid1).getSurface();
                            auto surface2 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid2).getSurface();

                            area[i][j][k] = surface1.getSurfaceArea() + surface2.getSurfaceArea();

                            if (visualize) {
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
        /*
        void generateSheet(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            //int stencil_size, double coeff_stddev = 0.1, double min_thickness = machineZero, double max_thickness = 0.5, double thickness_stddev = 0.0, 
            double min_thickness,
            double max_thickness,
            bool subgrid = true,
            bool variable_thickness = true, 
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

                Eigen::Vector3d datum;
                
                if (subgrid && variable_thickness) {
                    // Random datum anywhere in central cell to make sure it is small here
                    datum = generateRandomPoint(
                        -0.5,
                        0.5, eng);
                } else {
                    // Random datum anywhere in stencil
                    datum = generateRandomPoint(
                        -2.5,
                        2.5, eng);
                }

                // Random unbiased direction
                Eigen::Vector3d direction = generateRandomDirection(eng);

                // Random sheet thickness
                double thickness = max_thickness;
                if (sheet_thickness_stddev > 0.0) {
                    thickness = sample_truncated_normal(0, sheet_thickness_stddev, min_thickness, max_thickness);
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

                double coeff1 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);
                double coeff2 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);
                double coeff1p1 = coeff1;
                double coeff1p2 = coeff1;
                double coeff2p1 = coeff2;
                double coeff2p2 = coeff2;
                // Random coefficients
                if (variable_thickness) {
                    // If variable thickness, make sure coeffs are correlated so that sheet doesn't self-intersect
                    // add uniformly distributed noise
                    //std::uniform_real_distribution<double> noise_dist(0.0, thickness/4.0);

                    coeff1p1 += sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                    coeff1p2 -= sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);

                    coeff2p1 += sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                    coeff2p2 -= sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                }

                auto paraboloid1 = Paraboloid(datum_paraboloid1, frame, coeff1p1, coeff2p1);
                auto paraboloid2 = Paraboloid(datum_paraboloid2, frame, coeff1p2, coeff2p2);

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

                                // When in central cell, check if V1 is very small, meaning the central cell does not have 2 interfaces
                                
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
                            
                            auto surface1 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid1).getSurface();
                            auto surface2 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid2).getSurface();

                            area[i][j][k] = surface1.getSurfaceArea() + surface2.getSurfaceArea();

                            if (visualize) {
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
            std::vector<std::vector<std::vector<double>>>& area,
            bool cutInsideCentralCell,
            double min_thickness = 1e-12, 
            double max_thickness = 0.5, 
            bool subgrid = true,
            bool variable_thickness = false,
            Eigen::Matrix3d* secondMoment = nullptr) 
        {
            while (true) { // keep trying until center cell has surface crossing

                // Defining cell coordinates
                auto coords = makeCenteredCoords(stencil_size);
                const double cell_volume = 1.0;

                Eigen::Vector3d datum_eVec;
                
                if (subgrid && variable_thickness) {
                    // Random datum anywhere in central cell to make sure it is small here
                    datum_eVec = generateRandomPoint(
                        -0.5,
                        0.5, eng);
                } else {
                    // Random datum anywhere in stencil
                    datum_eVec = generateRandomPoint(
                        -2.5,
                        2.5, eng);
                }

                // Random datum anywhere in stencil
                const auto datum = Pt::fromRawDoublePointer(datum_eVec.data());

                // Random unbiased direction
                auto direction = Normal::fromRawDoublePointer(generateRandomDirection(eng).data());
                direction.normalize();

                // Random sheet thickness
                double thickness = max_thickness;
                if (sheet_thickness_stddev > 0.0) {
                    thickness = sample_truncated_normal(0, sheet_thickness_stddev, min_thickness, max_thickness);
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
                double coeff1 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);
                double coeff2 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);
                double coeff1p1 = coeff1;
                double coeff1p2 = coeff1;
                double coeff2p1 = coeff2;
                double coeff2p2 = coeff2;
                // Random coefficients
                if (variable_thickness) {
                    // If variable thickness, make sure coeffs are correlated so that sheet doesn't self-intersect
                    // add uniformly distributed noise
                    //std::uniform_real_distribution<double> noise_dist(0.0, thickness/4.0);

                    coeff1p1 += sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                    coeff1p2 -= sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);

                    coeff2p1 += sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                    coeff2p2 -= sample_truncated_normal(0.0, sheet_coeff_stddev/2, 0.0, 0.4);
                }

                auto paraboloid1 = Paraboloid(datum_paraboloid1, frame, coeff1p1, coeff2p1);
                auto paraboloid2 = Paraboloid(datum_paraboloid2, frame, coeff1p2, coeff2p2);

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

                            auto surface1 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid1).getSurface();
                            auto surface2 = getVolumeMoments<
                                VolumeMomentsAndSurface>(cell, paraboloid2).getSurface();

                            area[i][j][k] = surface1.getSurfaceArea() + surface2.getSurfaceArea();

                            if (visualize) {
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
        
        */

        // Function to generate a random vector within given bounds
        Eigen::Vector3d generateRandomPoint(double min_bound, double max_bound, std::mt19937_64& eng) {
            std::uniform_real_distribution<double> dist(min_bound, max_bound);
            return Eigen::Vector3d(dist(eng), dist(eng), dist(eng));  // Random point in 3D space
        }

        void compressStencilRefinedToCoarse(
            const std::vector<std::vector<std::vector<double>>>& volumes_refined,
            const std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoments_refined,
            const std::vector<std::vector<std::vector<double>>>& areas_refined,
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            int refinement_factor,
            double compressed_cell_volume,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>* centroid = nullptr)
        {
            double refinement_factor_double = static_cast<double>(refinement_factor);

            for (int i = 0; i < stencil_size; i++) {
                for (int j = 0; j < stencil_size; j++) {
                    for (int k = 0; k < stencil_size; k++) {
                        double vol_sum = 0.0;
                        Eigen::Vector3d firstMoment_sum = Eigen::Vector3d::Zero();
                        double area_sum = 0.0;

                        for (int m = i * refinement_factor; m < (i + 1) * refinement_factor; m++) {
                            for (int n = j * refinement_factor; n < (j + 1) * refinement_factor; n++) {
                                for (int o = k * refinement_factor; o < (k + 1) * refinement_factor; o++) {
                                    double v = volumes_refined[m][n][o];
                                    vol_sum += v;
                                    firstMoment_sum += firstMoments_refined[m][n][o];
                                    area_sum += areas_refined[m][n][o];
                                }
                            }
                        }


                        if (vol_sum > machineZero) {
                            firstMoment[i][j][k] = firstMoment_sum;
                            vfrac[i][j][k] = vol_sum / compressed_cell_volume;
                            area[i][j][k] = area_sum;
                        } else {
                            firstMoment[i][j][k] = Eigen::Vector3d::Zero();
                            vfrac[i][j][k] = 0.0;
                            area[i][j][k] = 0.0;
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

        void generateSheetTransition(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            bool thick_central_cell = false,
            double min_thickness = 1e-12,
            double max_thickness = 0.5,
            double max_thick_thickness = 2.5,
            Eigen::Matrix3d* secondMoment = nullptr)
        {
            if (max_thickness < min_thickness) {
                throw std::invalid_argument("generateSheetTransition: max_thickness must be >= min_thickness.");
            }
            if (max_thick_thickness < max_thickness) {
                throw std::invalid_argument("generateSheetTransition: max_thick_thickness must be >= max_thickness.");
            }

            while (true) {
                // Coarse mesh
                auto coords_coarse = makeCenteredCoords(stencil_size);
                const double coarse_cell_volume = 1.0;

                // Fixed refinement factor requested for this generator
                const int refinement_factor = 3;
                const int refined_stencil_size = refinement_factor * stencil_size;

                std::vector<double> coords_refined(refined_stencil_size + 1);
                for (int i = 0; i <= refined_stencil_size; ++i) {
                    coords_refined[i] =
                        -0.5 * static_cast<double>(stencil_size)
                        + static_cast<double>(i) / static_cast<double>(refinement_factor);
                }

                // Center cell
                const int mid = stencil_size / 2;

                // Only used for visualization
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                std::vector<ParaboloidParametrizedSurfaceOutput> surfaces;

                // Random datum anywhere in stencil
                const auto datum = Pt::fromRawDoublePointer(generateRandomPoint(-2.5,2.5,eng).data());

                // Random unbiased direction
                auto direction = Normal::fromRawDoublePointer(generateRandomDirection(eng).data());
                direction.normalize();

                // Uniformly distributed thin and thick thicknesses
                std::uniform_real_distribution<double> thin_thickness_dist(min_thickness, max_thickness);
                std::uniform_real_distribution<double> thick_thickness_dist(max_thickness, max_thick_thickness);

                const double thin_thickness = thin_thickness_dist(eng);
                const double thick_thickness = thick_thickness_dist(eng);

                // Build orthonormal frame aligned with direction
                const auto frame = ReferenceFrame::fromNormal(direction);

                // Same coeff1 and coeff2 for both paraboloids.
                // No variable_thickness coeff perturbations here.
                const double coeff1 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);
                const double coeff2 = sample_truncated_normal(0.0, sheet_coeff_stddev, -1.0, 1.0);

                auto makeParaboloidPair = [&](double selected_thickness)
                    -> std::pair<Paraboloid, Paraboloid>
                {
                    const auto datum_paraboloid1 =
                        Pt(datum - direction * (selected_thickness / 2.0));
                    const auto datum_paraboloid2 =
                        Pt(datum + direction * (selected_thickness / 2.0));

                    return std::pair<Paraboloid, Paraboloid>(
                        Paraboloid(datum_paraboloid1, frame, coeff1, coeff2),
                        Paraboloid(datum_paraboloid2, frame, coeff1, coeff2)
                    );
                };

                // Use the central-side thickness to construct the reference outer paraboloid
                // used only for choosing the transition plane, analogous to generateCutSheet.
                const double central_reference_thickness =
                    thick_central_cell ? thick_thickness : thin_thickness;

                const auto reference_paraboloid2_datum =
                    Pt(datum + direction * (central_reference_thickness / 2.0));

                auto reference_paraboloid2 =
                    Paraboloid(reference_paraboloid2_datum, frame, coeff1, coeff2);

                // Random point anywhere in stencil
                const auto sample_pt = Pt::fromRawDoublePointer(generateRandomPoint(
                    -0.5 * static_cast<double>(stencil_size),
                    0.5 * static_cast<double>(stencil_size),
                    eng).data());

                // Project sample point onto the reference outer paraboloid
                auto tmp_paraboloid = reference_paraboloid2;
                tmp_paraboloid.regenerateAtLocation(sample_pt);
                const auto new_datum = tmp_paraboloid.getDatum();

                // Choose a random tangent direction on that paraboloid.
                // This becomes the normal of the transition plane, as in generateCutSheet.
                std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
                const double theta = angle_dist(eng);

                IRL::Normal transition_normal =
                    std::cos(theta) * tmp_paraboloid.getReferenceFrame()[0] +
                    std::sin(theta) * tmp_paraboloid.getReferenceFrame()[1];

                transition_normal.normalize();

                const double transition_plane_constant = transition_normal * new_datum;

                // Never allow the transition plane to cut the central coarse cell.
                // This replaces the old cutInsideCentralCell logic.
                double min_signed = 1.0e300;
                double max_signed = -1.0e300;

                for (int a = 0; a <= 1; ++a) {
                    for (int b = 0; b <= 1; ++b) {
                        for (int c = 0; c <= 1; ++c) {
                            const Pt corner(
                                coords_coarse[mid + a],
                                coords_coarse[mid + b],
                                coords_coarse[mid + c]);

                            const double signed_distance =
                                transition_normal * corner - transition_plane_constant;

                            min_signed = std::min(min_signed, signed_distance);
                            max_signed = std::max(max_signed, signed_distance);
                        }
                    }
                }

                // Reject if the plane intersects or nearly touches the central cell.
                if (min_signed <= machineZero && max_signed >= -machineZero) {
                    continue;
                }

                const double origin_signed = -transition_plane_constant;
                if (std::abs(origin_signed) <= machineZero) {
                    continue;
                }

                auto isOnOriginSideOfTransitionPlane =
                    [&](const Eigen::Vector3d& p) -> bool
                {
                    const Pt pt(p.x(), p.y(), p.z());
                    const double signed_distance =
                        transition_normal * pt - transition_plane_constant;

                    return signed_distance * origin_signed >= 0.0;
                };

                auto usesThickThickness =
                    [&](const Eigen::Vector3d& p) -> bool
                {
                    const bool on_origin_side = isOnOriginSideOfTransitionPlane(p);

                    return  ( thick_central_cell &&  on_origin_side) ||
                            (!thick_central_cell && !on_origin_side);
                };

                auto selectedThickness =
                    [&](const Eigen::Vector3d& p) -> double
                {
                    return usesThickThickness(p) ? thick_thickness : thin_thickness;
                };

                // Quick central-cell check before paying for the full refined stencil.
                auto central_cell = RectangularCuboid::fromBoundingPts(
                    Pt(coords_coarse[mid],     coords_coarse[mid],     coords_coarse[mid]),
                    Pt(coords_coarse[mid + 1], coords_coarse[mid + 1], coords_coarse[mid + 1]));

                const Eigen::Vector3d central_cell_center(
                    0.5 * (coords_coarse[mid] + coords_coarse[mid + 1]),
                    0.5 * (coords_coarse[mid] + coords_coarse[mid + 1]),
                    0.5 * (coords_coarse[mid] + coords_coarse[mid + 1]));

                const bool central_uses_thick = usesThickThickness(central_cell_center);
                const double central_thickness = selectedThickness(central_cell_center);

                const auto central_pair = makeParaboloidPair(central_thickness);

                auto central_moments1 =
                    getVolumeMoments<VolumeMoments>(central_cell, central_pair.first);
                auto central_moments2 =
                    getVolumeMoments<VolumeMoments>(central_cell, central_pair.second);

                const double central_V1 = central_moments1.volume();
                const double central_V2 = central_moments2.volume();
                const double central_Vsheet = central_V2 - central_V1;

                bool central_cell_ok = false;

                //Changed below: No need to require two interfaces in the central cell
                //if (central_uses_thick) {
                    // Thick central cell: only require a partial volume fraction.
                    central_cell_ok =
                        central_Vsheet > machineZero &&
                        central_Vsheet < coarse_cell_volume - machineZero;
                /*} else {
                    // Thin central cell: require two interfaces.
                    //
                    // The requested subtraction test is:
                    // Vsheet = V2 - V1 should differ from V2.
                    // That means V1 is nonzero, so the first paraboloid contributes.
                    //
                    // Also require V2 not to be saturated, so the second paraboloid is
                    // actually cutting the central cell instead of lying completely past it.
                    const bool subtraction_changes_from_V2 =
                        std::abs(central_Vsheet - central_V2) > machineZero;

                    const bool outer_paraboloid_not_saturated =
                        central_V2 < coarse_cell_volume - machineZero;

                    central_cell_ok =
                        central_Vsheet > machineZero &&
                        central_Vsheet < coarse_cell_volume - machineZero &&
                        subtraction_changes_from_V2 &&
                        outer_paraboloid_not_saturated;
                }*/

                if (!central_cell_ok) {
                    continue;
                }

                // Refined fields
                std::vector<std::vector<std::vector<double>>> volumes_refined(
                    refined_stencil_size,
                    std::vector<std::vector<double>>(
                        refined_stencil_size,
                        std::vector<double>(refined_stencil_size, 0.0)));

                std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(
                    refined_stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        refined_stencil_size,
                        std::vector<Eigen::Vector3d>(
                            refined_stencil_size,
                            Eigen::Vector3d::Zero())));

                std::vector<std::vector<std::vector<double>>> surface_areas_refined(
                    refined_stencil_size,
                    std::vector<std::vector<double>>(
                        refined_stencil_size,
                        std::vector<double>(refined_stencil_size, 0.0)));

                using VolumeMomentsAndSurface =
                    AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

                IRL::GeneralMoments3D<2> totalMoments_refined =
                    IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);

                // Fill refined stencil.
                for (int i = 0; i < refined_stencil_size; ++i) {
                    for (int j = 0; j < refined_stencil_size; ++j) {
                        for (int k = 0; k < refined_stencil_size; ++k) {
                            auto cell = RectangularCuboid::fromBoundingPts(
                                Pt(coords_refined[i],     coords_refined[j],     coords_refined[k]),
                                Pt(coords_refined[i + 1], coords_refined[j + 1], coords_refined[k + 1]));

                            const Eigen::Vector3d cell_center(
                                0.5 * (coords_refined[i] + coords_refined[i + 1]),
                                0.5 * (coords_refined[j] + coords_refined[j + 1]),
                                0.5 * (coords_refined[k] + coords_refined[k + 1]));

                            const double local_thickness = selectedThickness(cell_center);
                            const auto paraboloid_pair = makeParaboloidPair(local_thickness);

                            auto volume_and_surface1 =
                                getVolumeMoments<VolumeMomentsAndSurface>(
                                    cell,
                                    paraboloid_pair.first);

                            auto volume_and_surface2 =
                                getVolumeMoments<VolumeMomentsAndSurface>(
                                    cell,
                                    paraboloid_pair.second);

                            const auto moments1 = volume_and_surface1.getMoments();
                            const auto moments2 = volume_and_surface2.getMoments();

                            const double V1 = moments1.volume();
                            const double V2 = moments2.volume();
                            const double Vdiff = V2 - V1;

                            const Eigen::Vector3d M1(
                                moments1.centroid().x(),
                                moments1.centroid().y(),
                                moments1.centroid().z());

                            const Eigen::Vector3d M2(
                                moments2.centroid().x(),
                                moments2.centroid().y(),
                                moments2.centroid().z());

                            auto surface1 = volume_and_surface1.getSurface();
                            auto surface2 = volume_and_surface2.getSurface();

                            surface_areas_refined[i][j][k] =
                                surface1.getSurfaceArea() + surface2.getSurfaceArea();

                            if (visualize) {
                                surfaces.push_back(surface1);
                                surfaces.push_back(surface2);
                            }

                            if (Vdiff <= machineZero) {
                                volumes_refined[i][j][k] = 0.0;
                                firstMoments_refined[i][j][k] = Eigen::Vector3d::Zero();
                            } else {
                                volumes_refined[i][j][k] = Vdiff;
                                firstMoments_refined[i][j][k] = M2 - M1;
                            }

                            if (secondMoment != nullptr) {
                                auto gm1 =
                                    IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(
                                        cell,
                                        paraboloid_pair.first);

                                auto gm2 =
                                    IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(
                                        cell,
                                        paraboloid_pair.second);

                                auto gmSheet = gm2 - gm1;

                                if (gmSheet.volume() > machineZero) {
                                    totalMoments_refined += gmSheet;
                                }
                            }
                        }
                    }
                }

                // Compress refined stencil back to coarse coordinates.
                compressStencilRefinedToCoarse(
                    volumes_refined,
                    firstMoments_refined,
                    surface_areas_refined,
                    vfrac,
                    firstMoment,
                    area,
                    refinement_factor,
                    coarse_cell_volume,
                    &centroid
                );

                if (secondMoment != nullptr) {
                    *secondMoment = centeredSecondMomentFromTotal(totalMoments_refined);
                }

                if (visualize) {
                    std::cout << "Writing coarse stencil to files..." << std::endl;
                    WriteField(stencil_size, coords_coarse, vfrac, "vfrac");
                    WriteSurface(surfaces, "surface");
                    printCentroids(centroid);
                }

                return;
            }
        }
        // This is an old function, might need some adjustments to get to the newer functions
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
                                double cell_area = 0.0;

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
                                    cell_area += M_PI * radius * radius * L;
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

                    std::vector<std::vector<std::vector<double>>> surface_areas_refined(refined_stencil_size,
                        std::vector<std::vector<double>>(refined_stencil_size,
                        std::vector<double>(refined_stencil_size)));

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
                    /*
                    // Compress refined → coarse stencil
                    compressStencilRefinedToCoarse(
                        volumes_refined,
                        firstMoments_refined,
                        surface_areas_refined,
                        vfrac,
                        firstMoment,
                        area,
                        refinement_factor,
                        cell_volume,
                        &centroid
                    );
                    */
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
            std::vector<std::vector<std::vector<double>>>& area,
            double min_radius = 1e-12, double max_radius = 0.5, 
            //double radius_circle_min = 2.5, double radius_circle_max = 10.0,
            //bool visualize = false,
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
                if (cylinder_radius_stddev > 0.0) {
                    tube_radius = sample_truncated_normal(0, cylinder_radius_stddev, min_radius, max_radius);
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
                                double cell_area = 0.0;

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
                                    cell_area += 2 * M_PI * tube_radius * L;
                                }

                                // Convert to your stored per-cell format:
                                // vfrac = volume fraction (cell volume = 1)
                                // firstMoment stores centroid (your refinement code does that)
                                if (cellV > 0.0) {
                                    const double vf = std::min(cellV / cell_volume, 1.0);
                                    vfrac[i][j][k] = vf;

                                    // centroid of union-of-cylinders approximation in this cell
                                    firstMoment[i][j][k] = cellM1;

                                    area[i][j][k] = cell_area;

                                    if (secondMoment != nullptr) {
                                        // accumulate second moments to get them for the whole stencil
                                        totalV  += cellV;
                                        totalM1 += cellM1;
                                        totalM2 += cellM2;
                                    }
                                } else {
                                    vfrac[i][j][k] = 0.0;
                                    firstMoment[i][j][k].setZero();
                                    area[i][j][k] = 0.0;
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

                    std::vector<std::vector<std::vector<double>>> surface_areas_refined(refined_stencil_size,
                                            std::vector<std::vector<double>>(refined_stencil_size,
                                            std::vector<double>(refined_stencil_size)));                    

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

                                auto surface = volume_and_surface.getSurface();
                                surface_areas_refined[i][j][k] = surface.getSurfaceArea();

                                if (visualize) {
                                    surfaces.push_back(surface);
                                }
                            }
                        }
                    }

                    // Compress refined to coarse stencil
                    compressStencilRefinedToCoarse(
                        volumes_refined,
                        firstMoments_refined,
                        surface_areas_refined,
                        vfrac,
                        firstMoment,
                        area,
                        refinement_factor,
                        cell_volume,
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

        // old function, not used anymore. Using it again may need some minor adjustments to make it similar to the newer functions
        void generateTruncatedCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            //double min_radius = machineZero, double max_radius = 0.5,
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
                                double cell_surface_area = 0.0;

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
                                    cell_surface_area += 2.0 * M_PI * radius * L; // lateral surface area
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
                    /*
                    // Compress refined → coarse stencil
                    compressStencilRefinedToCoarse(
                        volumes_refined,
                        firstMoments_refined,

                        vfrac,
                        firstMoment,
                        refinement_factor,
                        cell_volume,
                        &centroid
                    );
                    */
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
            std::vector<std::vector<std::vector<double>>>& area,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& centroid,
            std::vector<IRL::ParaboloidParametrizedSurfaceOutput>& surfaces,
            const Eigen::Vector3d& origin,
            double radius,
            std::vector<double>& coarse_coords,
            //bool visualize = false,
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
                            area[i][j][k] = 0.0;
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

                area[ci][cj][ck] = 4.0 * M_PI * radius * radius;

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
                        area[coarse_i][coarse_j][coarse_k] = 0.0;
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

                                    auto surface = volume_and_surface.getSurface();

                                    area[coarse_i][coarse_j][coarse_k] += surface.getSurfaceArea();

                                    // Accumulate exact raw 0th/1st/2nd moments if requested
                                    if (secondMoment != nullptr) {
                                        auto refined_general_moments =
                                            IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(refined_cell, paraboloid);
                                        total_general_moments += refined_general_moments;
                                    }

                                    if (visualize) {
                                        surfaces.push_back(surface);
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
            std::vector<std::vector<std::vector<double>>>& area,
            double min_radius = 1e-12, double max_radius = 0.5, //double radius_stddev = 0.0, bool visualize = false,
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
                    if (sphere_radius_stddev > 0.0) {
                        radius = sample_truncated_normal(0, sphere_radius_stddev, min_radius, max_radius);
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
                        area,
                        centroid,
                        surfaces,
                        origin,
                        radius,
                        coarse_coords
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

        //Below: ellipsoid helpers
        struct AxisAlignedBoundingBox {
            Eigen::Vector3d min;
            Eigen::Vector3d max;
        };

        static inline IRL::Pt eigenToIRLPt(const Eigen::Vector3d& p)
        {
            return IRL::Pt(p.x(), p.y(), p.z());
        }

        static inline IRL::Normal eigenToIRLNormal(const Eigen::Vector3d& n)
        {
            return IRL::Normal(n.x(), n.y(), n.z());
        }

        static inline IRL::Plane planeFromNormalAndPoint(
            const Eigen::Vector3d& normal,
            const Eigen::Vector3d& point_on_plane)
        {
            Eigen::Vector3d n = normal;

            const double norm = n.norm();
            if (norm < 1.0e-14) {
                throw std::runtime_error(
                    "planeFromNormalAndPoint: plane normal has near-zero magnitude."
                );
            }

            n /= norm;

            // IRL Plane convention:
            //
            //     n dot x - d = 0
            //
            // so:
            //
            //     d = n dot point_on_plane
            const double distance = n.dot(point_on_plane);

            return IRL::Plane(eigenToIRLNormal(n), distance);
        }

        static inline IRL::PlanarLocalizer makeEllipsoidAlignedBoxLocalizer(
            const Eigen::Vector3d& origin,
            const Eigen::Matrix3d& R,
            double u0,
            double u1,
            double v0,
            double v1,
            double w0,
            double w1)
        {
            const Eigen::Vector3d e0 = R.col(0);
            const Eigen::Vector3d e1 = R.col(1);
            const Eigen::Vector3d e2 = R.col(2);

            IRL::PlanarLocalizer localizer;
            localizer.setNumberOfPlanes(6);

            // Local coordinates:
            //
            //     x = origin + u e0 + v e1 + w e2
            //
            // The six planes are outward-facing. The inside of the localizer
            // should lie on the negative side of all six planes.

            const Eigen::Vector3d p_u0 = origin + u0 * e0;
            const Eigen::Vector3d p_u1 = origin + u1 * e0;

            const Eigen::Vector3d p_v0 = origin + v0 * e1;
            const Eigen::Vector3d p_v1 = origin + v1 * e1;

            const Eigen::Vector3d p_w0 = origin + w0 * e2;
            const Eigen::Vector3d p_w1 = origin + w1 * e2;

            // Face order:
            //
            //   0 = -u face
            //   1 = +u face
            //   2 = -v face
            //   3 = +v face
            //   4 = -w face
            //   5 = +w face

            localizer[0] = planeFromNormalAndPoint(-e0, p_u0);
            localizer[1] = planeFromNormalAndPoint( e0, p_u1);

            localizer[2] = planeFromNormalAndPoint(-e1, p_v0);
            localizer[3] = planeFromNormalAndPoint( e1, p_v1);

            localizer[4] = planeFromNormalAndPoint(-e2, p_w0);
            localizer[5] = planeFromNormalAndPoint( e2, p_w1);

            return localizer;
        }

        static inline AxisAlignedBoundingBox makeWorldAABBOfEllipsoidAlignedBox(
            const Eigen::Vector3d& origin,
            const Eigen::Matrix3d& R,
            double u0,
            double u1,
            double v0,
            double v1,
            double w0,
            double w1)
        {
            AxisAlignedBoundingBox box;

            box.min = Eigen::Vector3d::Constant(
                std::numeric_limits<double>::infinity()
            );

            box.max = Eigen::Vector3d::Constant(
                -std::numeric_limits<double>::infinity()
            );

            const std::array<double, 2> us = {u0, u1};
            const std::array<double, 2> vs = {v0, v1};
            const std::array<double, 2> ws = {w0, w1};

            for (double u : us) {
                for (double v : vs) {
                    for (double w : ws) {
                        const Eigen::Vector3d p =
                            origin + R * Eigen::Vector3d(u, v, w);

                        box.min = box.min.cwiseMin(p);
                        box.max = box.max.cwiseMax(p);
                    }
                }
            }

            return box;
        }

        static inline bool aabbOverlap(
            const AxisAlignedBoundingBox& a,
            const AxisAlignedBoundingBox& b,
            double tolerance = 1.0e-14)
        {
            return
                (a.min.x() <= b.max.x() + tolerance) &&
                (a.max.x() + tolerance >= b.min.x()) &&
                (a.min.y() <= b.max.y() + tolerance) &&
                (a.max.y() + tolerance >= b.min.y()) &&
                (a.min.z() <= b.max.z() + tolerance) &&
                (a.max.z() + tolerance >= b.min.z());
        }

        static inline std::vector<double> makePaddedDiameterGridCoordinates(
            double semi_axis,
            int cells_across_diameter,
            double padding)
        {
            if (semi_axis <= 0.0) {
                throw std::runtime_error(
                    "makePaddedDiameterGridCoordinates: semi_axis must be positive."
                );
            }

            if (cells_across_diameter <= 0) {
                throw std::runtime_error(
                    "makePaddedDiameterGridCoordinates: cells_across_diameter must be positive."
                );
            }

            if (padding < 0.0) {
                throw std::runtime_error(
                    "makePaddedDiameterGridCoordinates: padding must be non-negative."
                );
            }

            std::vector<double> coords;

            // Exactly `cells_across_diameter` cells over [-semi_axis, semi_axis],
            // plus one padding cell of thickness `padding` on each side.
            //
            // For cells_across_diameter = 5:
            //
            //   [-a-0.1, -a, -a+2a/5, -a+4a/5,
            //    -a+6a/5, -a+8a/5, a, a+0.1]
            //
            // This gives 5 cells across the ellipsoid diameter and two extra
            // thin padding cells.

            if (padding > 1.0e-14) {
                coords.push_back(-semi_axis - padding);
            }

            coords.push_back(-semi_axis);

            for (int i = 1; i < cells_across_diameter; ++i) {
                const double s =
                    -semi_axis +
                    2.0 * semi_axis *
                    static_cast<double>(i) /
                    static_cast<double>(cells_across_diameter);

                coords.push_back(s);
            }

            coords.push_back(semi_axis);

            if (padding > 1.0e-14) {
                coords.push_back(semi_axis + padding);
            }

            return coords;
        }

        static inline IRL::Paraboloid makeEllipsoidParaboloidAtPoint(
            const Eigen::Vector3d& origin,
            const Eigen::Matrix3d& R,
            const Eigen::Matrix3d& A,
            const Eigen::Vector3d& point_for_direction)
        {
            Eigen::Vector3d direction = point_for_direction - origin;

            if (direction.norm() < 1.0e-14) {
                direction = R.col(0);
            }

            direction.normalize();

            // Ray/ellipsoid intersection:
            //
            //     datum = origin + t direction
            //     t = 1 / sqrt(direction^T A direction)

            const double ray_denominator =
                std::sqrt(direction.dot(A * direction));

            if (ray_denominator < 1.0e-14) {
                throw std::runtime_error(
                    "makeEllipsoidParaboloidAtPoint: invalid ray/ellipsoid intersection."
                );
            }

            const Eigen::Vector3d datum_eigen =
                origin + direction / ray_denominator;

            const IRL::Pt datum_paraboloid =
                eigenToIRLPt(datum_eigen);

            const Eigen::Vector3d centered_datum =
                datum_eigen - origin;

            Eigen::Vector3d normal =
                A * centered_datum;

            const double normal_denominator = normal.norm();

            if (normal_denominator < 1.0e-14) {
                throw std::runtime_error(
                    "makeEllipsoidParaboloidAtPoint: invalid ellipsoid normal."
                );
            }

            normal /= normal_denominator;

            const Eigen::Vector3d tangent0 =
                normal.unitOrthogonal();

            const Eigen::Vector3d tangent1 =
                normal.cross(tangent0).normalized();

            Eigen::Matrix2d curvature_matrix;

            curvature_matrix(0, 0) =
                tangent0.dot(A * tangent0) / normal_denominator;

            curvature_matrix(0, 1) =
                tangent0.dot(A * tangent1) / normal_denominator;

            curvature_matrix(1, 0) =
                curvature_matrix(0, 1);

            curvature_matrix(1, 1) =
                tangent1.dot(A * tangent1) / normal_denominator;

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>
                curvature_solver(curvature_matrix);

            const Eigen::Vector2d eigvec0 =
                curvature_solver.eigenvectors().col(0);

            const Eigen::Vector2d eigvec1 =
                curvature_solver.eigenvectors().col(1);

            Eigen::Vector3d paraboloid_x_axis =
                eigvec0.x() * tangent0 + eigvec0.y() * tangent1;

            Eigen::Vector3d paraboloid_y_axis =
                eigvec1.x() * tangent0 + eigvec1.y() * tangent1;

            paraboloid_x_axis.normalize();
            paraboloid_y_axis.normalize();

            if (paraboloid_x_axis.cross(paraboloid_y_axis).dot(normal) < 0.0) {
                paraboloid_y_axis *= -1.0;
            }

            const double coeff1 =
                0.5 * std::max(0.0, curvature_solver.eigenvalues()(0));

            const double coeff2 =
                0.5 * std::max(0.0, curvature_solver.eigenvalues()(1));

            const auto reference_frame = IRL::ReferenceFrame(
                IRL::Normal(
                    paraboloid_x_axis.x(),
                    paraboloid_x_axis.y(),
                    paraboloid_x_axis.z()
                ),
                IRL::Normal(
                    paraboloid_y_axis.x(),
                    paraboloid_y_axis.y(),
                    paraboloid_y_axis.z()
                ),
                IRL::Normal(
                    normal.x(),
                    normal.y(),
                    normal.z()
                )
            );

            return IRL::Paraboloid(
                datum_paraboloid,
                reference_frame,
                coeff1,
                coeff2
            );
        }

        //Note: this function does not return correct surfaces / areas
        void generateSpecificEllipsoid(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& centroid,
            std::vector<IRL::ParaboloidParametrizedSurfaceOutput>& surfaces,
            const Eigen::Vector3d& origin,
            const Eigen::Vector3d& axis0_direction,
            const Eigen::Vector3d& axis1_direction,
            const Eigen::Vector3d& axis2_direction,
            double axis0_length,
            double axis1_length,
            double axis2_length,
            std::vector<double>& coarse_coords,
            Eigen::Matrix3d* secondMoment = nullptr)
        {
            if (axis0_length <= machineZero ||
                axis1_length <= machineZero ||
                axis2_length <= machineZero) {
                throw std::runtime_error(
                    "generateSpecificEllipsoid: all axis lengths must be positive."
                );
            }

            const double a = axis0_length;
            const double b = axis1_length;
            const double c = axis2_length;

            // Build an orthonormal frame from the three input directions.
            //
            // axis0_direction corresponds to semi-axis a.
            // axis1_direction corresponds to semi-axis b.
            // axis2_direction is mainly used to choose handedness/sign.

            Eigen::Vector3d e0 = axis0_direction.normalized();

            Eigen::Vector3d e1 =
                axis1_direction - axis1_direction.dot(e0) * e0;

            if (e1.norm() < 1.0e-14) {
                e1 = e0.unitOrthogonal();
            } else {
                e1.normalize();
            }

            Eigen::Vector3d e2 = e0.cross(e1).normalized();

            if (axis2_direction.dot(e2) < 0.0) {
                e2 *= -1.0;
                e1 = e2.cross(e0).normalized();
            }

            Eigen::Matrix3d R;
            R.col(0) = e0;
            R.col(1) = e1;
            R.col(2) = e2;

            // Original stencil coordinates.
            const double coarse_stencil_min =
                -0.5 * static_cast<double>(stencil_size);

            for (int coarse_index = 0;
                coarse_index <= stencil_size;
                ++coarse_index) {
                coarse_coords[coarse_index] =
                    coarse_stencil_min + static_cast<double>(coarse_index);
            }

            auto clamp_coarse_index = [&](int idx) {
                return std::max(0, std::min(stencil_size - 1, idx));
            };

            // World-space implicit ellipsoid:
            //
            //     (x - origin)^T A (x - origin) = 1

            const Eigen::Vector3d inv_axes_sq(
                1.0 / (a * a),
                1.0 / (b * b),
                1.0 / (c * c)
            );

            const Eigen::Matrix3d A =
                R * inv_axes_sq.asDiagonal() * R.transpose();

            // Exact global x/y/z AABB of the rotated ellipsoid.

            Eigen::Vector3d half_extent;

            for (int d = 0; d < 3; ++d) {
                half_extent[d] = std::sqrt(
                    std::pow(R(d, 0) * a, 2) +
                    std::pow(R(d, 1) * b, 2) +
                    std::pow(R(d, 2) * c, 2)
                );
            }

            const double ellipsoid_min_x = origin.x() - half_extent.x();
            const double ellipsoid_max_x = origin.x() + half_extent.x();

            const double ellipsoid_min_y = origin.y() - half_extent.y();
            const double ellipsoid_max_y = origin.y() + half_extent.y();

            const double ellipsoid_min_z = origin.z() - half_extent.z();
            const double ellipsoid_max_z = origin.z() + half_extent.z();

            // Zero outputs.

            for (int i = 0; i < stencil_size; ++i) {
                for (int j = 0; j < stencil_size; ++j) {
                    for (int k = 0; k < stencil_size; ++k) {
                        vfrac[i][j][k] = 0.0;
                        firstMoment[i][j][k].setZero();
                        area[i][j][k] = 0.0;

                        if (visualize) {
                            centroid[i][j][k].setZero();
                        }
                    }
                }
            }

            IRL::GeneralMoments3D<2> total_general_moments =
                IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);

            // Analytical shortcut if the whole ellipsoid AABB fits inside one
            // original stencil cell.

            int ci = clamp_coarse_index(
                static_cast<int>(std::floor(origin.x() - coarse_stencil_min))
            );

            int cj = clamp_coarse_index(
                static_cast<int>(std::floor(origin.y() - coarse_stencil_min))
            );

            int ck = clamp_coarse_index(
                static_cast<int>(std::floor(origin.z() - coarse_stencil_min))
            );

            const double x0 = coarse_stencil_min + static_cast<double>(ci);
            const double x1 = x0 + 1.0;

            const double y0 = coarse_stencil_min + static_cast<double>(cj);
            const double y1 = y0 + 1.0;

            const double z0 = coarse_stencil_min + static_cast<double>(ck);
            const double z1 = z0 + 1.0;

            const bool ellipsoid_fully_inside_some_cell =
                (ellipsoid_min_x >= x0) && (ellipsoid_max_x <= x1) &&
                (ellipsoid_min_y >= y0) && (ellipsoid_max_y <= y1) &&
                (ellipsoid_min_z >= z0) && (ellipsoid_max_z <= z1);

            if (ellipsoid_fully_inside_some_cell) {
                const double V = (4.0 / 3.0) * M_PI * a * b * c;

                vfrac[ci][cj][ck] = V;
                firstMoment[ci][cj][ck] = V * origin;

                if (visualize) {
                    centroid[ci][cj][ck] = origin;
                }

                // Knud Thomsen approximation for ellipsoid surface area.
                // This shortcut still fills area, because it does not use the
                // localized-link path.
                const double p = 1.6075;

                area[ci][cj][ck] =
                    4.0 * M_PI *
                    std::pow(
                        (
                            std::pow(a, p) * std::pow(b, p) +
                            std::pow(a, p) * std::pow(c, p) +
                            std::pow(b, p) * std::pow(c, p)
                        ) / 3.0,
                        1.0 / p
                    );

                if (secondMoment != nullptr) {
                    const Eigen::Vector3d axes_sq(a * a, b * b, c * c);

                    // Centered second moment tensor of a solid ellipsoid.
                    *secondMoment =
                        (V / 5.0) *
                        R * axes_sq.asDiagonal() * R.transpose();
                }

                return;
            }

            // ------------------------------------------------------------------
            // Ellipsoid-aligned refined grid
            // ------------------------------------------------------------------
            //
            // New behavior:
            //
            //   1. Choose a target physical refined-cell size from the smallest
            //      ellipsoid diameter.
            //
            //   2. Use more cells along longer ellipsoid axes so that the refined
            //      cell sizes are roughly uniform in physical space.
            //
            //   3. Limit the total number of refined cells to avoid excessive cost.
            //
            //   4. Interpret safety_distance as a margin between the ellipsoid
            //      surface bounding interval and the grid boundary.
            //
            //      The grid spans:
            //
            //          [-a - safety_distance, a + safety_distance]
            //
            //      directly divided into num_u cells. There is no separate padding
            //      cell anymore.

            const int cells_across_smallest_diameter = 5;
            const int minimum_cells_per_direction = 3;
            const int maximum_cells_per_direction = 32;
            const int maximum_total_refined_cells = 2500;
            const double safety_distance = 0.1;

            const double min_axis = std::min({a, b, c});

            double target_refined_cell_size =
                2.0 * min_axis /
                static_cast<double>(cells_across_smallest_diameter);

            auto chooseCellsAcrossPaddedInterval =
                [&](double semi_axis, double target_h) -> int
            {
                const double padded_span =
                    2.0 * (semi_axis + safety_distance);

                int n = static_cast<int>(
                    std::ceil(padded_span / target_h)
                );

                n = std::max(n, minimum_cells_per_direction);
                n = std::min(n, maximum_cells_per_direction);

                return n;
            };

            int refined_cells_u =
                chooseCellsAcrossPaddedInterval(a, target_refined_cell_size);

            int refined_cells_v =
                chooseCellsAcrossPaddedInterval(b, target_refined_cell_size);

            int refined_cells_w =
                chooseCellsAcrossPaddedInterval(c, target_refined_cell_size);

            auto refinedCellProduct = [&]() -> int {
                return refined_cells_u * refined_cells_v * refined_cells_w;
            };

            // Coarsen the target cell size if the first estimate is too expensive.
            //
            // The cube-root factor keeps the coarsening approximately isotropic in
            // the refined-cell-size sense.

            int refinement_adjustment_iterations = 0;

            while (refinedCellProduct() > maximum_total_refined_cells &&
                refinement_adjustment_iterations < 20) {
                const double current_total =
                    static_cast<double>(refinedCellProduct());

                const double scale =
                    1.01 *
                    std::cbrt(
                        current_total /
                        static_cast<double>(maximum_total_refined_cells)
                    );

                target_refined_cell_size *= scale;

                refined_cells_u =
                    chooseCellsAcrossPaddedInterval(a, target_refined_cell_size);

                refined_cells_v =
                    chooseCellsAcrossPaddedInterval(b, target_refined_cell_size);

                refined_cells_w =
                    chooseCellsAcrossPaddedInterval(c, target_refined_cell_size);

                ++refinement_adjustment_iterations;
            }

            auto makeSafetyMarginGridCoordinates =
                [&](double semi_axis, int number_of_cells) -> std::vector<double>
            {
                const double lower = -semi_axis - safety_distance;
                const double upper =  semi_axis + safety_distance;

                std::vector<double> coords(number_of_cells + 1, 0.0);

                for (int index = 0; index <= number_of_cells; ++index) {
                    coords[index] =
                        lower +
                        (upper - lower) *
                        static_cast<double>(index) /
                        static_cast<double>(number_of_cells);
                }

                return coords;
            };

            const std::vector<double> u_coords =
                makeSafetyMarginGridCoordinates(a, refined_cells_u);

            const std::vector<double> v_coords =
                makeSafetyMarginGridCoordinates(b, refined_cells_v);

            const std::vector<double> w_coords =
                makeSafetyMarginGridCoordinates(c, refined_cells_w);

            const int num_u = static_cast<int>(u_coords.size()) - 1;
            const int num_v = static_cast<int>(v_coords.size()) - 1;
            const int num_w = static_cast<int>(w_coords.size()) - 1;

            const int number_of_refined_cells =
                num_u * num_v * num_w;

            if (visualize) {
                std::cout << "Ellipsoid aligned refinement: "
                        << num_u << " x "
                        << num_v << " x "
                        << num_w << " = "
                        << number_of_refined_cells
                        << " cells\n";

                std::cout << "Target refined cell size: "
                        << target_refined_cell_size << "\n";
            }

            std::vector<IRL::PlanarLocalizer> refined_localizers;
            std::vector<IRL::SeparatorVariant> refined_interfaces;
            std::vector<AxisAlignedBoundingBox> refined_aabbs;

            refined_localizers.reserve(number_of_refined_cells);
            refined_interfaces.reserve(number_of_refined_cells);
            refined_aabbs.reserve(number_of_refined_cells);

            for (int iu = 0; iu < num_u; ++iu) {
                const double u0_local = u_coords[iu];
                const double u1_local = u_coords[iu + 1];

                for (int iv = 0; iv < num_v; ++iv) {
                    const double v0_local = v_coords[iv];
                    const double v1_local = v_coords[iv + 1];

                    for (int iw = 0; iw < num_w; ++iw) {
                        const double w0_local = w_coords[iw];
                        const double w1_local = w_coords[iw + 1];

                        refined_localizers.push_back(
                            makeEllipsoidAlignedBoxLocalizer(
                                origin,
                                R,
                                u0_local,
                                u1_local,
                                v0_local,
                                v1_local,
                                w0_local,
                                w1_local
                            )
                        );

                        const Eigen::Vector3d local_cell_center(
                            0.5 * (u0_local + u1_local),
                            0.5 * (v0_local + v1_local),
                            0.5 * (w0_local + w1_local)
                        );

                        const Eigen::Vector3d world_cell_center =
                            origin + R * local_cell_center;

                        refined_interfaces.push_back(
                            makeEllipsoidParaboloidAtPoint(
                                origin,
                                R,
                                A,
                                world_cell_center
                            )
                        );

                        refined_aabbs.push_back(
                            makeWorldAABBOfEllipsoidAlignedBox(
                                origin,
                                R,
                                u0_local,
                                u1_local,
                                v0_local,
                                v1_local,
                                w0_local,
                                w1_local
                            )
                        );
                    }
                }
            }

            const AxisAlignedBoundingBox aligned_grid_aabb =
                makeWorldAABBOfEllipsoidAlignedBox(
                    origin,
                    R,
                    u_coords.front(),
                    u_coords.back(),
                    v_coords.front(),
                    v_coords.back(),
                    w_coords.front(),
                    w_coords.back()
                );

            const int coarse_i_start = clamp_coarse_index(
                static_cast<int>(
                    std::floor(aligned_grid_aabb.min.x() - coarse_stencil_min)
                )
            );

            const int coarse_i_end = clamp_coarse_index(
                static_cast<int>(
                    std::floor(aligned_grid_aabb.max.x() - coarse_stencil_min)
                )
            );

            const int coarse_j_start = clamp_coarse_index(
                static_cast<int>(
                    std::floor(aligned_grid_aabb.min.y() - coarse_stencil_min)
                )
            );

            const int coarse_j_end = clamp_coarse_index(
                static_cast<int>(
                    std::floor(aligned_grid_aabb.max.y() - coarse_stencil_min)
                )
            );

            const int coarse_k_start = clamp_coarse_index(
                static_cast<int>(
                    std::floor(aligned_grid_aabb.min.z() - coarse_stencil_min)
                )
            );

            const int coarse_k_end = clamp_coarse_index(
                static_cast<int>(
                    std::floor(aligned_grid_aabb.max.z() - coarse_stencil_min)
                )
            );

            for (int coarse_i = coarse_i_start;
                coarse_i <= coarse_i_end;
                ++coarse_i) {
                const double coarse_x0 =
                    coarse_stencil_min + static_cast<double>(coarse_i);

                const double coarse_x1 = coarse_x0 + 1.0;

                for (int coarse_j = coarse_j_start;
                    coarse_j <= coarse_j_end;
                    ++coarse_j) {
                    const double coarse_y0 =
                        coarse_stencil_min + static_cast<double>(coarse_j);

                    const double coarse_y1 = coarse_y0 + 1.0;

                    for (int coarse_k = coarse_k_start;
                        coarse_k <= coarse_k_end;
                        ++coarse_k) {
                        const double coarse_z0 =
                            coarse_stencil_min + static_cast<double>(coarse_k);

                        const double coarse_z1 = coarse_z0 + 1.0;

                        const auto coarse_cell =
                            IRL::RectangularCuboid::fromBoundingPts(
                                IRL::Pt(coarse_x0, coarse_y0, coarse_z0),
                                IRL::Pt(coarse_x1, coarse_y1, coarse_z1)
                            );

                        AxisAlignedBoundingBox coarse_aabb;
                        coarse_aabb.min =
                            Eigen::Vector3d(coarse_x0, coarse_y0, coarse_z0);
                        coarse_aabb.max =
                            Eigen::Vector3d(coarse_x1, coarse_y1, coarse_z1);

                        for (int refined_id = 0;
                            refined_id < number_of_refined_cells;
                            ++refined_id) {
                            if (!aabbOverlap(
                                    coarse_aabb,
                                    refined_aabbs[refined_id])) {
                                continue;
                            }

                            IRL::LocalizedSeparatorVariantLink localized_interface(
                                &(refined_localizers[refined_id]),
                                &(refined_interfaces[refined_id])
                            );

                            localized_interface.setId(refined_id);

                            for (int face = 0; face < 6; ++face) {
                                localized_interface.setEdgeConnectivity(
                                    face,
                                    nullptr
                                );
                            }

                            const auto local_moments =
                                IRL::getVolumeMoments<IRL::VolumeMoments>(
                                    coarse_cell,
                                    localized_interface
                                );

                            const double local_volume =
                                local_moments.volume();

                            if (local_volume <= machineZero) {
                                continue;
                            }

                            const Eigen::Vector3d local_first_moment(
                                local_moments.centroid().x(),
                                local_moments.centroid().y(),
                                local_moments.centroid().z()
                            );

                            vfrac[coarse_i][coarse_j][coarse_k] +=
                                local_volume;

                            firstMoment[coarse_i][coarse_j][coarse_k] +=
                                local_first_moment;

                            // Surface area is intentionally not computed in this
                            // localized-link path.
                            area[coarse_i][coarse_j][coarse_k] += 0.0;

                            if (secondMoment != nullptr) {
                                const auto local_general_moments =
                                    IRL::getVolumeMoments<
                                        IRL::GeneralMoments3D<2>
                                    >(
                                        coarse_cell,
                                        localized_interface
                                    );

                                total_general_moments +=
                                    local_general_moments;
                            }
                        }
                    }
                }
            }

            if (visualize) {
                for (int i = 0; i < stencil_size; ++i) {
                    for (int j = 0; j < stencil_size; ++j) {
                        for (int k = 0; k < stencil_size; ++k) {
                            if (vfrac[i][j][k] > machineZero) {
                                centroid[i][j][k] =
                                    firstMoment[i][j][k] /
                                    vfrac[i][j][k];
                            } else {
                                centroid[i][j][k].setZero();
                            }
                        }
                    }
                }
            }

            if (secondMoment != nullptr) {
                *secondMoment =
                    centeredSecondMomentFromTotal(total_general_moments);
            }
        }

        void generateSpecificEllipsoidOLD(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& centroid,
            std::vector<IRL::ParaboloidParametrizedSurfaceOutput>& surfaces,
            const Eigen::Vector3d& origin,
            const Eigen::Vector3d& axis0_direction,
            const Eigen::Vector3d& axis1_direction,
            const Eigen::Vector3d& axis2_direction,
            double axis0_length,
            double axis1_length,
            double axis2_length,
            std::vector<double>& coarse_coords,
            Eigen::Matrix3d* secondMoment = nullptr)
        {
            if (axis0_length <= machineZero ||
                axis1_length <= machineZero ||
                axis2_length <= machineZero) {
                throw std::runtime_error("generateSpecificEllipsoid: all axis lengths must be positive.");
            }

            const double a = axis0_length;
            const double b = axis1_length;
            const double c = axis2_length;

            const Eigen::Vector3d semi_axes(a, b, c);
            const double min_semi_axis = semi_axes.minCoeff();

            // Build an orthonormal frame from the three input directions.
            //
            // axis0_direction corresponds to length a
            // axis1_direction corresponds to length b
            // axis2_direction is used mainly to choose the handedness/sign

            Eigen::Vector3d e0 = axis0_direction.normalized();

            Eigen::Vector3d e1 =
                axis1_direction - axis1_direction.dot(e0) * e0;

            if (e1.norm() < 1.0e-14) {
                e1 = e0.unitOrthogonal();
            } else {
                e1.normalize();
            }

            Eigen::Vector3d e2 = e0.cross(e1).normalized();

            // Try to preserve the sign/orientation suggested by axis2_direction.
            if (axis2_direction.dot(e2) < 0.0) {
                e2 *= -1.0;
                e1 = e2.cross(e0).normalized();
            }

            Eigen::Matrix3d R;
            R.col(0) = e0;
            R.col(1) = e1;
            R.col(2) = e2;

            // Stencil coordinates
            const double coarse_stencil_min = -0.5 * static_cast<double>(stencil_size);

            for (int coarse_index = 0; coarse_index <= stencil_size; ++coarse_index) {
                coarse_coords[coarse_index] =
                    coarse_stencil_min + static_cast<double>(coarse_index);
            }

            auto clamp_coarse_index = [&](int idx) {
                return std::max(0, std::min(stencil_size - 1, idx));
            };

            // World-space implicit ellipsoid:
            //
            //     (x - origin)^T A (x - origin) = 1
            //
            // R stores the ellipsoid frame:
            //
            //     R.col(0) = direction of axis a
            //     R.col(1) = direction of axis b
            //     R.col(2) = direction of axis c

            const Eigen::Vector3d inv_axes_sq(
                1.0 / (a * a),
                1.0 / (b * b),
                1.0 / (c * c)
            );

            const Eigen::Matrix3d A =
                R * inv_axes_sq.asDiagonal() * R.transpose();

            // Exact axis-aligned bounding box half-widths of rotated ellipsoid
            Eigen::Vector3d half_extent;
            for (int d = 0; d < 3; ++d) {
                half_extent[d] = std::sqrt(
                    std::pow(R(d, 0) * a, 2) +
                    std::pow(R(d, 1) * b, 2) +
                    std::pow(R(d, 2) * c, 2)
                );
            }

            const double ellipsoid_min_x = origin.x() - half_extent.x();
            const double ellipsoid_max_x = origin.x() + half_extent.x();
            const double ellipsoid_min_y = origin.y() - half_extent.y();
            const double ellipsoid_max_y = origin.y() + half_extent.y();
            const double ellipsoid_min_z = origin.z() - half_extent.z();
            const double ellipsoid_max_z = origin.z() + half_extent.z();

            // Zero outputs
            for (int i = 0; i < stencil_size; ++i) {
                for (int j = 0; j < stencil_size; ++j) {
                    for (int k = 0; k < stencil_size; ++k) {
                        vfrac[i][j][k] = 0.0;
                        firstMoment[i][j][k].setZero();
                        area[i][j][k] = 0.0;

                        if (visualize) {
                            centroid[i][j][k].setZero();
                        }
                    }
                }
            }

            using VolumeMomentsAndSurface =
                AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;

            IRL::GeneralMoments3D<2> total_general_moments =
                IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);

            // Analytical shortcut if the whole ellipsoid AABB fits inside one coarse cell
            int ci = clamp_coarse_index(
                static_cast<int>(std::floor(origin.x() - coarse_stencil_min))
            );
            int cj = clamp_coarse_index(
                static_cast<int>(std::floor(origin.y() - coarse_stencil_min))
            );
            int ck = clamp_coarse_index(
                static_cast<int>(std::floor(origin.z() - coarse_stencil_min))
            );

            const double x0 = coarse_stencil_min + static_cast<double>(ci);
            const double x1 = x0 + 1.0;
            const double y0 = coarse_stencil_min + static_cast<double>(cj);
            const double y1 = y0 + 1.0;
            const double z0 = coarse_stencil_min + static_cast<double>(ck);
            const double z1 = z0 + 1.0;

            const bool ellipsoid_fully_inside_some_cell =
                (ellipsoid_min_x >= x0) && (ellipsoid_max_x <= x1) &&
                (ellipsoid_min_y >= y0) && (ellipsoid_max_y <= y1) &&
                (ellipsoid_min_z >= z0) && (ellipsoid_max_z <= z1);

            if (ellipsoid_fully_inside_some_cell) {
                const double V = (4.0 / 3.0) * M_PI * a * b * c;

                vfrac[ci][cj][ck] = V;
                firstMoment[ci][cj][ck] = V * origin;

                // Knud Thomsen approximation for ellipsoid surface area
                const double p = 1.6075;
                area[ci][cj][ck] =
                    4.0 * M_PI *
                    std::pow(
                        (
                            std::pow(a, p) * std::pow(b, p) +
                            std::pow(a, p) * std::pow(c, p) +
                            std::pow(b, p) * std::pow(c, p)
                        ) / 3.0,
                        1.0 / p
                    );

                if (secondMoment != nullptr) {
                    const Eigen::Vector3d axes_sq(a * a, b * b, c * c);

                    // Centered second moment tensor of a solid ellipsoid
                    *secondMoment =
                        (V / 5.0) *
                        R * axes_sq.asDiagonal() * R.transpose();
                }

                return;
            }

            // Refined paraboloid intersection path
            const int min_refined_cells_across_smallest_diameter = 5;

            const int refinement_factor = std::max(
                1,
                static_cast<int>(
                    std::ceil(
                        static_cast<double>(min_refined_cells_across_smallest_diameter) /
                        (2.0 * min_semi_axis)
                    )
                )
            );

            const double refined_cell_size =
                1.0 / static_cast<double>(refinement_factor);

            auto clamp_refined_index = [&](int index_value) {
                return std::max(0, std::min(refinement_factor - 1, index_value));
            };

            int coarse_i_start = clamp_coarse_index(
                static_cast<int>(std::floor(ellipsoid_min_x - coarse_stencil_min))
            );
            int coarse_i_end = clamp_coarse_index(
                static_cast<int>(std::floor(ellipsoid_max_x - coarse_stencil_min))
            );

            int coarse_j_start = clamp_coarse_index(
                static_cast<int>(std::floor(ellipsoid_min_y - coarse_stencil_min))
            );
            int coarse_j_end = clamp_coarse_index(
                static_cast<int>(std::floor(ellipsoid_max_y - coarse_stencil_min))
            );

            int coarse_k_start = clamp_coarse_index(
                static_cast<int>(std::floor(ellipsoid_min_z - coarse_stencil_min))
            );
            int coarse_k_end = clamp_coarse_index(
                static_cast<int>(std::floor(ellipsoid_max_z - coarse_stencil_min))
            );

            for (int coarse_i = coarse_i_start; coarse_i <= coarse_i_end; ++coarse_i) {
                const double coarse_x0 =
                    coarse_stencil_min + static_cast<double>(coarse_i);

                int refined_i_start = clamp_refined_index(
                    static_cast<int>(
                        std::floor((ellipsoid_min_x - coarse_x0) / refined_cell_size)
                    )
                );

                int refined_i_end = clamp_refined_index(
                    static_cast<int>(
                        std::floor((ellipsoid_max_x - coarse_x0) / refined_cell_size)
                    )
                );

                if (refined_i_end < refined_i_start) continue;

                for (int coarse_j = coarse_j_start; coarse_j <= coarse_j_end; ++coarse_j) {
                    const double coarse_y0 =
                        coarse_stencil_min + static_cast<double>(coarse_j);

                    int refined_j_start = clamp_refined_index(
                        static_cast<int>(
                            std::floor((ellipsoid_min_y - coarse_y0) / refined_cell_size)
                        )
                    );

                    int refined_j_end = clamp_refined_index(
                        static_cast<int>(
                            std::floor((ellipsoid_max_y - coarse_y0) / refined_cell_size)
                        )
                    );

                    if (refined_j_end < refined_j_start) continue;

                    for (int coarse_k = coarse_k_start; coarse_k <= coarse_k_end; ++coarse_k) {
                        const double coarse_z0 =
                            coarse_stencil_min + static_cast<double>(coarse_k);

                        int refined_k_start = clamp_refined_index(
                            static_cast<int>(
                                std::floor((ellipsoid_min_z - coarse_z0) / refined_cell_size)
                            )
                        );

                        int refined_k_end = clamp_refined_index(
                            static_cast<int>(
                                std::floor((ellipsoid_max_z - coarse_z0) / refined_cell_size)
                            )
                        );

                        if (refined_k_end < refined_k_start) continue;

                        for (int refined_i_in_coarse = refined_i_start;
                            refined_i_in_coarse <= refined_i_end;
                            ++refined_i_in_coarse) {

                            const double refined_x0 =
                                coarse_x0 +
                                static_cast<double>(refined_i_in_coarse) * refined_cell_size;

                            const double refined_x1 = refined_x0 + refined_cell_size;

                            for (int refined_j_in_coarse = refined_j_start;
                                refined_j_in_coarse <= refined_j_end;
                                ++refined_j_in_coarse) {

                                const double refined_y0 =
                                    coarse_y0 +
                                    static_cast<double>(refined_j_in_coarse) * refined_cell_size;

                                const double refined_y1 = refined_y0 + refined_cell_size;

                                for (int refined_k_in_coarse = refined_k_start;
                                    refined_k_in_coarse <= refined_k_end;
                                    ++refined_k_in_coarse) {

                                    const double refined_z0 =
                                        coarse_z0 +
                                        static_cast<double>(refined_k_in_coarse) *
                                        refined_cell_size;

                                    const double refined_z1 =
                                        refined_z0 + refined_cell_size;

                                    auto refined_cell = IRL::RectangularCuboid::fromBoundingPts(
                                        IRL::Pt(refined_x0, refined_y0, refined_z0),
                                        IRL::Pt(refined_x1, refined_y1, refined_z1)
                                    );

                                    const Eigen::Vector3d refined_cell_center(
                                        0.5 * (refined_x0 + refined_x1),
                                        0.5 * (refined_y0 + refined_y1),
                                        0.5 * (refined_z0 + refined_z1)
                                    );

                                    Eigen::Vector3d direction =
                                        refined_cell_center - origin;

                                    if (direction.norm() < 1.0e-14) {
                                        direction = R.col(0);
                                    }

                                    direction.normalize();

                                    // Ray/ellipsoid intersection:
                                    //
                                    //     datum = origin + t direction
                                    //     t = 1 / sqrt(direction^T A direction)
                                    //
                                    const double ray_denominator =
                                        std::sqrt(direction.dot(A * direction));

                                    const Eigen::Vector3d datum_eigen =
                                        origin + direction / ray_denominator;

                                    const IRL::Pt datum_paraboloid(
                                        datum_eigen.x(),
                                        datum_eigen.y(),
                                        datum_eigen.z()
                                    );

                                    const Eigen::Vector3d centered_datum =
                                        datum_eigen - origin;

                                    Eigen::Vector3d normal =
                                        A * centered_datum;

                                    const double normal_denominator = normal.norm();

                                    if (normal_denominator < 1.0e-14) {
                                        continue;
                                    }

                                    normal /= normal_denominator;

                                    const Eigen::Vector3d tangent0 =
                                        normal.unitOrthogonal();

                                    const Eigen::Vector3d tangent1 =
                                        normal.cross(tangent0).normalized();

                                    Eigen::Matrix2d curvature_matrix;
                                    curvature_matrix(0, 0) =
                                        tangent0.dot(A * tangent0) / normal_denominator;
                                    curvature_matrix(0, 1) =
                                        tangent0.dot(A * tangent1) / normal_denominator;
                                    curvature_matrix(1, 0) =
                                        curvature_matrix(0, 1);
                                    curvature_matrix(1, 1) =
                                        tangent1.dot(A * tangent1) / normal_denominator;

                                    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>
                                        curvature_solver(curvature_matrix);

                                    const Eigen::Vector2d eigvec0 =
                                        curvature_solver.eigenvectors().col(0);

                                    const Eigen::Vector2d eigvec1 =
                                        curvature_solver.eigenvectors().col(1);

                                    Eigen::Vector3d paraboloid_x_axis =
                                        eigvec0.x() * tangent0 + eigvec0.y() * tangent1;

                                    Eigen::Vector3d paraboloid_y_axis =
                                        eigvec1.x() * tangent0 + eigvec1.y() * tangent1;

                                    paraboloid_x_axis.normalize();
                                    paraboloid_y_axis.normalize();

                                    if (paraboloid_x_axis.cross(paraboloid_y_axis).dot(normal) < 0.0) {
                                        paraboloid_y_axis *= -1.0;
                                    }

                                    const double coeff1 =
                                        0.5 * std::max(0.0, curvature_solver.eigenvalues()(0));

                                    const double coeff2 =
                                        0.5 * std::max(0.0, curvature_solver.eigenvalues()(1));

                                    const auto reference_frame = IRL::ReferenceFrame(
                                        IRL::Normal(
                                            paraboloid_x_axis.x(),
                                            paraboloid_x_axis.y(),
                                            paraboloid_x_axis.z()
                                        ),
                                        IRL::Normal(
                                            paraboloid_y_axis.x(),
                                            paraboloid_y_axis.y(),
                                            paraboloid_y_axis.z()
                                        ),
                                        IRL::Normal(
                                            normal.x(),
                                            normal.y(),
                                            normal.z()
                                        )
                                    );

                                    const auto paraboloid = IRL::Paraboloid(
                                        datum_paraboloid,
                                        reference_frame,
                                        coeff1,
                                        coeff2
                                    );

                                    auto volume_and_surface =
                                        getVolumeMoments<VolumeMomentsAndSurface>(
                                            refined_cell,
                                            paraboloid
                                        );

                                    const double refined_volume =
                                        volume_and_surface.getMoments().volume();

                                    const Eigen::Vector3d refined_first_moment(
                                        volume_and_surface.getMoments().centroid().x(),
                                        volume_and_surface.getMoments().centroid().y(),
                                        volume_and_surface.getMoments().centroid().z()
                                    );

                                    vfrac[coarse_i][coarse_j][coarse_k] += refined_volume;

                                    firstMoment[coarse_i][coarse_j][coarse_k] +=
                                        refined_first_moment;

                                    auto surface = volume_and_surface.getSurface();

                                    area[coarse_i][coarse_j][coarse_k] +=
                                        surface.getSurfaceArea();

                                    if (secondMoment != nullptr) {
                                        auto refined_general_moments =
                                            IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(
                                                refined_cell,
                                                paraboloid
                                            );

                                        total_general_moments += refined_general_moments;
                                    }

                                    if (visualize) {
                                        surfaces.push_back(surface);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (secondMoment != nullptr) {
                *secondMoment = centeredSecondMomentFromTotal(total_general_moments);
            }
        }

        void generateEllipsoidDroplet(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            double min_axis = 1e-12,
            double max_axis = 0.8,
            Eigen::Matrix3d* secondMoment = nullptr)
        {
            if (min_axis < machineZero || max_axis <= min_axis) {
                throw std::runtime_error("generateEllipsoidDroplet: invalid axis bounds.");
            }

            auto random_rotation_matrix = [&]() {
                std::uniform_real_distribution<double> random_unit(0.0, 1.0);

                const double u1 = random_unit(eng);
                const double u2 = random_unit(eng);
                const double u3 = random_unit(eng);

                const double qx = std::sqrt(1.0 - u1) * std::sin(2.0 * M_PI * u2);
                const double qy = std::sqrt(1.0 - u1) * std::cos(2.0 * M_PI * u2);
                const double qz = std::sqrt(u1)       * std::sin(2.0 * M_PI * u3);
                const double qw = std::sqrt(u1)       * std::cos(2.0 * M_PI * u3);

                Eigen::Quaterniond q(qw, qx, qy, qz);
                q.normalize();

                return q.toRotationMatrix();
            };

            while (true) {
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;
                std::vector<double> coarse_coords(stencil_size + 1);

                double a = max_axis;
                double b = max_axis;
                double c = max_axis;

                if (ellipsoid_subgrid_stddev > 0.0) {
                    // Sample axis lengths from truncated normal distributions
                    a = sample_truncated_normal(max_axis, ellipsoid_subgrid_stddev, min_axis, max_axis);
                    b = sample_truncated_normal(a, ellipsoid_subgrid_stddev*0.5, min_axis, max_axis);
                    c = sample_truncated_normal(a, ellipsoid_subgrid_stddev*0.5, min_axis, max_axis);
                } else {
                    std::uniform_real_distribution<double> random_axis(min_axis, max_axis);
                    a = random_axis(eng);
                    b = random_axis(eng);
                    c = random_axis(eng);
                }

                const Eigen::Matrix3d R = random_rotation_matrix();

                Eigen::Vector3d half_extent;
                for (int d = 0; d < 3; ++d) {
                    half_extent[d] = std::sqrt(
                        std::pow(R(d, 0) * a, 2) +
                        std::pow(R(d, 1) * b, 2) +
                        std::pow(R(d, 2) * c, 2)
                    );
                }

                // Origin within central cell plus 1/2 of the rotated ellipsoid AABB
                // in each coordinate direction.
                std::uniform_real_distribution<double> random_x(
                    -0.5 - half_extent.x(),
                    0.5 + half_extent.x()
                );

                std::uniform_real_distribution<double> random_y(
                    -0.5 - half_extent.y(),
                    0.5 + half_extent.y()
                );

                std::uniform_real_distribution<double> random_z(
                    -0.5 - half_extent.z(),
                    0.5 + half_extent.z()
                );

                const Eigen::Vector3d origin(
                    random_x(eng),
                    random_y(eng),
                    random_z(eng)
                );

                generateSpecificEllipsoid(
                    vfrac,
                    firstMoment,
                    area,
                    centroid,
                    surfaces,
                    origin,
                    R.col(0),
                    R.col(1),
                    R.col(2),
                    a,
                    b,
                    c,
                    coarse_coords,
                    secondMoment
                );

                const int mid = stencil_size / 2;
                const double center_vfrac = vfrac[mid][mid][mid];

                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {
                    if (visualize) {
                        WriteField(stencil_size, coarse_coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                    }

                    return;
                }
            }
        }

        void generateEllipsoidLigament(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            double min_long_axis = 1.5,
            double max_long_axis = 2.5,
            double min_short_axis = 1e-12,
            double max_short_axis = 0.8,
            Eigen::Matrix3d* secondMoment = nullptr)
        {
            if (min_long_axis < machineZero ||
                max_long_axis <= min_long_axis ||
                min_short_axis < machineZero ||
                max_short_axis <= min_short_axis) {
                throw std::runtime_error("generateEllipsoidLigament: invalid axis bounds.");
            }

            auto random_rotation_matrix = [&]() {
                std::uniform_real_distribution<double> random_unit(0.0, 1.0);

                const double u1 = random_unit(eng);
                const double u2 = random_unit(eng);
                const double u3 = random_unit(eng);

                const double qx = std::sqrt(1.0 - u1) * std::sin(2.0 * M_PI * u2);
                const double qy = std::sqrt(1.0 - u1) * std::cos(2.0 * M_PI * u2);
                const double qz = std::sqrt(u1)       * std::sin(2.0 * M_PI * u3);
                const double qw = std::sqrt(u1)       * std::cos(2.0 * M_PI * u3);

                Eigen::Quaterniond q(qw, qx, qy, qz);
                q.normalize();

                return q.toRotationMatrix();
            };

            while (true) {
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;
                std::vector<double> coarse_coords(stencil_size + 1);

                double a = max_long_axis;
                double b = max_short_axis;
                double c = max_short_axis;

                if (ellipsoid_subgrid_stddev > 0.0) {
                    // Sample axis lengths from truncated normal distributions
                    b = sample_truncated_normal(max_short_axis, ellipsoid_subgrid_stddev, min_short_axis, max_short_axis);
                    c = sample_truncated_normal(max_short_axis, ellipsoid_subgrid_stddev, min_short_axis, max_short_axis);
                } else {
                    std::uniform_real_distribution<double> random_axis(min_short_axis, max_short_axis);
                    b = random_axis(eng);
                    c = random_axis(eng);
                }

                std::uniform_real_distribution<double> random_long_axis(
                    min_long_axis,
                    max_long_axis
                );
                a = random_long_axis(eng);

                const double shortest_axis = std::min(b, c);

                const Eigen::Matrix3d R = random_rotation_matrix();

                // Origin within central cell plus shortest axis.
                std::uniform_real_distribution<double> random_origin_coord(
                    -0.5 - shortest_axis,
                    0.5 + shortest_axis
                );

                const Eigen::Vector3d origin(
                    random_origin_coord(eng),
                    random_origin_coord(eng),
                    random_origin_coord(eng)
                );

                generateSpecificEllipsoid(
                    vfrac,
                    firstMoment,
                    area,
                    centroid,
                    surfaces,
                    origin,
                    R.col(0),
                    R.col(1),
                    R.col(2),
                    a,
                    b,
                    c,
                    coarse_coords,
                    secondMoment
                );

                const int mid = stencil_size / 2;
                const double center_vfrac = vfrac[mid][mid][mid];

                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {
                    if (visualize) {
                        WriteField(stencil_size, coarse_coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                    }

                    return;
                }
            }
        }

        void generateEllipsoidLigamentTip(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            double min_long_axis = 1.5,
            double max_long_axis = 2.5,
            double min_short_axis = 1e-12,
            double max_short_axis = 0.8,
            bool truncated_inside_central_cell = true,
            Eigen::Matrix3d* secondMoment = nullptr
            )
        {
            if (min_long_axis < machineZero ||
                max_long_axis <= min_long_axis ||
                min_short_axis < machineZero ||
                max_short_axis <= min_short_axis) {
                throw std::runtime_error(
                    "generateEllipsoidLigamentTip: invalid axis bounds."
                );
            }

            auto random_rotation_matrix = [&]() {
                std::uniform_real_distribution<double> random_unit(0.0, 1.0);

                const double u1 = random_unit(eng);
                const double u2 = random_unit(eng);
                const double u3 = random_unit(eng);

                const double qx =
                    std::sqrt(1.0 - u1) * std::sin(2.0 * M_PI * u2);
                const double qy =
                    std::sqrt(1.0 - u1) * std::cos(2.0 * M_PI * u2);
                const double qz =
                    std::sqrt(u1) * std::sin(2.0 * M_PI * u3);
                const double qw =
                    std::sqrt(u1) * std::cos(2.0 * M_PI * u3);

                Eigen::Quaterniond q(qw, qx, qy, qz);
                q.normalize();

                return q.toRotationMatrix();
            };

            // Returns the positive distance from point in direction until
            // the ray exits the axis-aligned box [box_min, box_max]^3.
            auto distance_to_box_exit = [](
                const Eigen::Vector3d& point,
                const Eigen::Vector3d& direction,
                double box_min,
                double box_max)
            {
                double distance = std::numeric_limits<double>::infinity();

                for (int d = 0; d < 3; ++d) {
                    if (direction[d] > 1.0e-14) {
                        distance = std::min(
                            distance,
                            (box_max - point[d]) / direction[d]
                        );
                    } else if (direction[d] < -1.0e-14) {
                        distance = std::min(
                            distance,
                            (box_min - point[d]) / direction[d]
                        );
                    }
                }

                return distance;
            };

            while (true) {
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(
                            stencil_size,
                            Eigen::Vector3d::Zero()
                        )
                    )
                );

                std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;
                std::vector<double> coarse_coords(stencil_size + 1);

                std::uniform_real_distribution<double> random_long_axis(
                    min_long_axis,
                    max_long_axis
                );

                const double a = random_long_axis(eng);

                double b = max_short_axis;
                double c = max_short_axis;

                if (ellipsoid_subgrid_stddev > 0.0) {
                    b = sample_truncated_normal(
                        max_short_axis,
                        ellipsoid_subgrid_stddev,
                        min_short_axis,
                        max_short_axis
                    );

                    c = sample_truncated_normal(
                        max_short_axis,
                        ellipsoid_subgrid_stddev,
                        min_short_axis,
                        max_short_axis
                    );
                } else {
                    std::uniform_real_distribution<double> random_short_axis(
                        min_short_axis,
                        max_short_axis
                    );

                    b = random_short_axis(eng);
                    c = random_short_axis(eng);
                }

                const Eigen::Matrix3d R = random_rotation_matrix();

                const Eigen::Vector3d long_axis_direction = R.col(0);
                const Eigen::Vector3d b_direction = R.col(1);
                const Eigen::Vector3d c_direction = R.col(2);

                Eigen::Vector3d origin;

                if (truncated_inside_central_cell) {

                    const Eigen::Vector3d ligament_tip =
                        generateRandomPoint(-0.5, 0.5, eng);

                    origin =
                        ligament_tip + a * long_axis_direction;
                } else {
                    // Place the ellipsoid tip outside the central cell while
                    // keeping the ligament axis passing through a randomly
                    // selected point inside the central cell.
                    //
                    // The ellipsoid center is moved out of the central cell in
                    // the negative long-axis direction. The corresponding
                    // negative-long-axis tip is therefore even farther outside:
                    //
                    //     ligament_tip =
                    //         origin - a * long_axis_direction.
                    //
                    // Limiting shift to less than a guarantees that the selected
                    // central-axis point remains inside the ellipsoid.

                    const Eigen::Vector3d central_axis_point =
                        generateRandomPoint(-0.5, 0.5, eng);

                    const Eigen::Vector3d negative_long_axis_direction =
                        -long_axis_direction;

                    const double distance_to_leave_central_cell =
                        distance_to_box_exit(
                            central_axis_point,
                            negative_long_axis_direction,
                            -0.5,
                            0.5
                        );

                    const double stencil_min =
                        -0.5 * static_cast<double>(stencil_size);

                    const double stencil_max =
                        0.5 * static_cast<double>(stencil_size);

                    const double distance_to_stencil_edge =
                        distance_to_box_exit(
                            central_axis_point,
                            negative_long_axis_direction,
                            stencil_min,
                            stencil_max
                        );

                    // Ensure that the ellipsoid center is strictly outside the
                    // central cell rather than exactly on its boundary.
                    const double min_shift =
                        distance_to_leave_central_cell + 1.0e-10;

                    // The center must remain less than one semi-major axis away
                    // from central_axis_point. This guarantees that the selected
                    // point lies strictly inside the ellipsoid along its long
                    // axis.
                    const double max_shift =
                        std::min(
                            distance_to_stencil_edge,
                            0.999 * a
                        );

                    // For some orientations and small values of a, it may not
                    // be possible to move the center outside the central cell
                    // while keeping central_axis_point inside the ellipsoid.
                    if (max_shift <= min_shift) {
                        continue;
                    }

                    std::uniform_real_distribution<double> random_shift(
                        min_shift,
                        max_shift
                    );

                    const double shift = random_shift(eng);

                    origin =
                        central_axis_point
                        + shift * negative_long_axis_direction;
                }

                generateSpecificEllipsoid(
                    vfrac,
                    firstMoment,
                    area,
                    centroid,
                    surfaces,
                    origin,
                    long_axis_direction,
                    b_direction,
                    c_direction,
                    a,
                    b,
                    c,
                    coarse_coords,
                    secondMoment
                );

                const int mid = stencil_size / 2;
                const double center_vfrac = vfrac[mid][mid][mid];

                // Retain only configurations whose interface actually cuts the central cell

                if (center_vfrac > machineZero &&
                    center_vfrac < 1.0 - machineZero) {
                    if (visualize) {
                        WriteField(
                            stencil_size,
                            coarse_coords,
                            vfrac,
                            "vfrac"
                        );

                        WriteSurface(
                            surfaces,
                            "surface"
                        );
                    }

                    return;
                }
            }
        }

        void generateEllipsoidSheet(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            double min_long_axis = 8.66,
            double max_long_axis = 15.0,
            double min_short_axis = 1e-12,
            double max_short_axis = 0.8,
            Eigen::Matrix3d* secondMoment = nullptr)
        {
            if (min_long_axis < machineZero ||
                max_long_axis <= min_long_axis ||
                min_short_axis < machineZero ||
                max_short_axis <= min_short_axis) {
                throw std::runtime_error("generateEllipsoidLigament: invalid axis bounds.");
            }

            auto random_rotation_matrix = [&]() {
                std::uniform_real_distribution<double> random_unit(0.0, 1.0);

                const double u1 = random_unit(eng);
                const double u2 = random_unit(eng);
                const double u3 = random_unit(eng);

                const double qx = std::sqrt(1.0 - u1) * std::sin(2.0 * M_PI * u2);
                const double qy = std::sqrt(1.0 - u1) * std::cos(2.0 * M_PI * u2);
                const double qz = std::sqrt(u1)       * std::sin(2.0 * M_PI * u3);
                const double qw = std::sqrt(u1)       * std::cos(2.0 * M_PI * u3);

                Eigen::Quaterniond q(qw, qx, qy, qz);
                q.normalize();

                return q.toRotationMatrix();
            };

            while (true) {
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;
                std::vector<double> coarse_coords(stencil_size + 1);

                double a = max_long_axis;
                double b = max_long_axis;
                double c = max_short_axis;

                if (ellipsoid_subgrid_stddev > 0.0) {
                    // Sample axis lengths from truncated normal distributions
                    c = sample_truncated_normal(max_short_axis, ellipsoid_subgrid_stddev, min_short_axis, max_short_axis);
                } else {
                    std::uniform_real_distribution<double> random_axis(min_short_axis, max_short_axis);
                    c = random_axis(eng);
                }

                std::uniform_real_distribution<double> random_long_axis(
                    min_long_axis,
                    max_long_axis
                );
                a = random_long_axis(eng);
                b = random_long_axis(eng);

                const double shortest_axis = c;

                const Eigen::Matrix3d R = random_rotation_matrix();

                // Origin within central cell plus shortest axis.
                std::uniform_real_distribution<double> random_origin_coord(
                    -0.5 - shortest_axis,
                    0.5 + shortest_axis
                );

                const Eigen::Vector3d origin(
                    random_origin_coord(eng),
                    random_origin_coord(eng),
                    random_origin_coord(eng)
                );

                generateSpecificEllipsoid(
                    vfrac,
                    firstMoment,
                    area,
                    centroid,
                    surfaces,
                    origin,
                    R.col(0),
                    R.col(1),
                    R.col(2),
                    a,
                    b,
                    c,
                    coarse_coords,
                    secondMoment
                );

                const int mid = stencil_size / 2;
                const double center_vfrac = vfrac[mid][mid][mid];

                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {
                    if (visualize) {
                        WriteField(stencil_size, coarse_coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                    }

                    return;
                }
            }
        }

        void generateEllipsoidSheetEdge(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            double min_large_axis = 3.0,
            double max_large_axis = 5.0,
            double min_short_axis = 1e-12,
            double max_short_axis = 1.6,
            bool cut_inside_central_cell = true,
            Eigen::Matrix3d* secondMoment = nullptr)
        {
            if (min_large_axis < machineZero ||
                max_large_axis <= min_large_axis ||
                min_short_axis < machineZero ||
                max_short_axis <= min_short_axis) {
                throw std::runtime_error("generate_ellipsoid_sheet_edge: invalid axis bounds.");
            }

            auto random_rotation_matrix = [&]() {
                std::uniform_real_distribution<double> random_unit(0.0, 1.0);

                const double u1 = random_unit(eng);
                const double u2 = random_unit(eng);
                const double u3 = random_unit(eng);

                const double qx = std::sqrt(1.0 - u1) * std::sin(2.0 * M_PI * u2);
                const double qy = std::sqrt(1.0 - u1) * std::cos(2.0 * M_PI * u2);
                const double qz = std::sqrt(u1)       * std::sin(2.0 * M_PI * u3);
                const double qw = std::sqrt(u1)       * std::cos(2.0 * M_PI * u3);

                Eigen::Quaterniond q(qw, qx, qy, qz);
                q.normalize();

                return q.toRotationMatrix();
            };

            auto distance_to_box_exit = [](
                const Eigen::Vector3d& point,
                const Eigen::Vector3d& direction,
                double box_min,
                double box_max)
            {
                double distance = std::numeric_limits<double>::infinity();

                for (int d = 0; d < 3; ++d) {
                    if (direction[d] > 1.0e-14) {
                        distance = std::min(
                            distance,
                            (box_max - point[d]) / direction[d]
                        );
                    } else if (direction[d] < -1.0e-14) {
                        distance = std::min(
                            distance,
                            (box_min - point[d]) / direction[d]
                        );
                    }
                }

                return distance;
            };

            while (true) {
                std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                    stencil_size,
                    std::vector<std::vector<Eigen::Vector3d>>(
                        stencil_size,
                        std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                    )
                );

                std::vector<IRL::ParaboloidParametrizedSurfaceOutput> surfaces;
                std::vector<double> coarse_coords(stencil_size + 1);

                std::uniform_real_distribution<double> random_large_axis(
                    min_large_axis,
                    max_large_axis
                );

                double a = random_large_axis(eng);
                double b = random_large_axis(eng);
                double c = max_short_axis;

                if (ellipsoid_subgrid_stddev > 0.0) {
                    c = sample_truncated_normal(
                        max_short_axis,
                        ellipsoid_subgrid_stddev,
                        min_short_axis,
                        max_short_axis
                    );
                } else {
                    std::uniform_real_distribution<double> random_short_axis(
                        min_short_axis,
                        max_short_axis
                    );

                    c = random_short_axis(eng);
                }

                // Ranfom Orientation
                const Eigen::Matrix3d R = random_rotation_matrix();

                const Eigen::Vector3d a_direction = R.col(0);
                const Eigen::Vector3d b_direction = R.col(1);
                const Eigen::Vector3d c_direction = R.col(2);

                Eigen::Vector3d origin;

                // Origin Placement
                if (cut_inside_central_cell) {
                    const Eigen::Vector3d sheet_edge_point =
                        generateRandomPoint(-0.5, 0.5, eng);

                    origin = sheet_edge_point + a * a_direction;
                } else {
                    // Sheet plane cuts the central cell, but the ellipsoid
                    // center is shifted in the negative a-direction.
                    //
                    // We first choose a point in the central cell. This point
                    // lies in the sheet plane spanned by a_direction and
                    // b_direction.
                    //
                    // Then we move the origin negatively along a_direction.
                    //
                    // Minimum distance:
                    //     far enough that the origin leaves the central cell.
                    //
                    // Maximum distance:
                    //     as far as possible before reaching the stencil edge.

                    const Eigen::Vector3d central_plane_point =
                        generateRandomPoint(-0.5, 0.5, eng);

                    const Eigen::Vector3d negative_a_direction =
                        -a_direction;

                    const double distance_to_leave_central_cell =
                        distance_to_box_exit(
                            central_plane_point,
                            negative_a_direction,
                            -0.5,
                            0.5
                        );

                    const double stencil_min =
                        -0.5 * static_cast<double>(stencil_size);

                    const double stencil_max =
                        0.5 * static_cast<double>(stencil_size);

                    const double distance_to_stencil_edge =
                        distance_to_box_exit(
                            central_plane_point,
                            negative_a_direction,
                            stencil_min,
                            stencil_max
                        );

                    // Need a small positive buffer so the origin is actually
                    // outside the central cell, not exactly on its boundary.
                    const double min_shift =
                        distance_to_leave_central_cell + 1.0e-10;

                    // To keep the chosen central point inside the ellipsoid
                    // along the a-direction, the shift should not exceed a.
                    //
                    // The requested geometric maximum is the stencil edge.
                    // This cap keeps the construction useful/rejects less often.
                    const double max_shift =
                        std::min(distance_to_stencil_edge, 0.999 * a);

                    if (max_shift <= min_shift) {
                        continue;
                    }

                    std::uniform_real_distribution<double> random_shift(
                        min_shift,
                        max_shift
                    );

                    const double shift = random_shift(eng);

                    origin = central_plane_point + shift * negative_a_direction;

                    // Equivalent:
                    //
                    //     origin = central_plane_point - shift * a_direction
                    //
                    // so the center is outside the central cell in the
                    // negative-a direction, while central_plane_point remains
                    // in the sheet's a-b plane.
                }

                generateSpecificEllipsoid(
                    vfrac,
                    firstMoment,
                    area,
                    centroid,
                    surfaces,
                    origin,
                    a_direction,
                    b_direction,
                    c_direction,
                    a,
                    b,
                    c,
                    coarse_coords,
                    secondMoment
                );

                const int mid = stencil_size / 2;
                const double center_vfrac = vfrac[mid][mid][mid];

                if (center_vfrac > machineZero && center_vfrac < 1.0 - machineZero) {
                    if (visualize) {
                        WriteField(stencil_size, coarse_coords, vfrac, "vfrac");
                        WriteSurface(surfaces, "surface");
                    }

                    return;
                }
            }
        }

        void generateBentTruncatedCylinder(
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            std::vector<std::vector<std::vector<double>>>& area,
            bool truncateInsideCentralCell, double min_radius = 1e-12, double max_radius = 0.5, //double radius_stddev = 0.0, 
            //double radius_circle_min = 2.5, double radius_circle_max = 10.0,
            //bool visualize = false,
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
                if (cylinder_radius_stddev > 0.0) {
                    tube_radius = sample_truncated_normal(0, cylinder_radius_stddev, min_radius, max_radius);
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

                // Refined mesh
                const double cell_volume = 1.0;
                double max_refinement_factor = 8.0;
                double refinement_factor_double = std::ceil(3.0/(2.0*tube_radius)); // want at least ~3 samples across the tube diameter for decent accuracy, can adjust this heuristic as needed
                int refinement_factor = static_cast<int>(refinement_factor_double);
                int refined_stencil_size = refinement_factor * stencil_size;
                bool approximate_with_cylinders = false;
                if (refinement_factor_double > max_refinement_factor) {
                    approximate_with_cylinders = true;
                    sphere_radius = 0.0; //In this case don't consider a sphere cap
                }

                // pick random theta_cut on full circle and accept only if it satisfies your constraints
                double theta_cut = 0.0;
                Eigen::Vector3d sphere_origin = Eigen::Vector3d::Zero();

                {
                    // Exact central coarse-cell box.
                    // For stencil_size = 5 this is [-0.5, 0.5]^3.
                    auto coords_for_center_cell = makeCenteredCoords(stencil_size);
                    const int mid = stencil_size / 2;

                    const Eigen::Vector3d central_min(
                        coords_for_center_cell[mid],
                        coords_for_center_cell[mid],
                        coords_for_center_cell[mid]);

                    const Eigen::Vector3d central_max(
                        coords_for_center_cell[mid + 1],
                        coords_for_center_cell[mid + 1],
                        coords_for_center_cell[mid + 1]);

                    auto pointAABBDistanceSq = [](
                        const Eigen::Vector3d& p,
                        const Eigen::Vector3d& bmin,
                        const Eigen::Vector3d& bmax)
                    {
                        double d2 = 0.0;

                        for (int a = 0; a < 3; ++a) {
                            if (p[a] < bmin[a]) {
                                const double d = bmin[a] - p[a];
                                d2 += d * d;
                            } else if (p[a] > bmax[a]) {
                                const double d = p[a] - bmax[a];
                                d2 += d * d;
                            }
                        }

                        return d2;
                    };

                    auto pointAABBFarthestDistanceSq = [](
                        const Eigen::Vector3d& p,
                        const Eigen::Vector3d& bmin,
                        const Eigen::Vector3d& bmax)
                    {
                        double d2 = 0.0;

                        for (int a = 0; a < 3; ++a) {
                            const double d_to_min = std::abs(p[a] - bmin[a]);
                            const double d_to_max = std::abs(p[a] - bmax[a]);
                            const double d = std::max(d_to_min, d_to_max);
                            d2 += d * d;
                        }

                        return d2;
                    };

                    auto sphereSurfaceCutsAABB = [&](
                        const Eigen::Vector3d& center,
                        const double radius,
                        const Eigen::Vector3d& bmin,
                        const Eigen::Vector3d& bmax)
                    {
                        const double r2 = radius * radius;
                        const double eps = 1.0e-12;

                        const double dmin2 = pointAABBDistanceSq(center, bmin, bmax);
                        const double dmax2 = pointAABBFarthestDistanceSq(center, bmin, bmax);

                        // Surface cuts the box iff the box contains points both inside and outside
                        // the sphere, including tangential tolerance.
                        return dmin2 <= r2 + eps && dmax2 >= r2 - eps;
                    };

                    auto sphereSolidIntersectsAABB = [&](
                        const Eigen::Vector3d& center,
                        const double radius,
                        const Eigen::Vector3d& bmin,
                        const Eigen::Vector3d& bmax)
                    {
                        const double r2 = radius * radius;
                        const double eps = 1.0e-12;

                        const double dmin2 = pointAABBDistanceSq(center, bmin, bmax);

                        // Solid sphere intersects the box iff the closest point in the box is
                        // within the radius.
                        return dmin2 <= r2 + eps;
                    };

                    // Try a few theta samples; if none work, restart outer while(true).
                    bool found = false;
                    std::uniform_real_distribution<double> dist_theta(-M_PI, M_PI);

                    for (int tries = 0; tries < 64; ++tries) {
                        const double t = dist_theta(eng);
                        const Eigen::Vector3d p_cut = centerlinePoint(t, c0, radius_circle, u, v);

                        // Keep cap center close enough to the stencil that the truncation is relevant.
                        const double half = 0.5 * static_cast<double>(stencil_size);
                        const double maxAbs = p_cut.cwiseAbs().maxCoeff();
                        if (maxAbs > half + sphere_radius) continue;

                        if (truncateInsideCentralCell) {
                            // Tip case:
                            // Require the spherical cap surface to actually cut the central cell.
                            // This excludes cap spheres that are fully outside or fully contain
                            // the central cell.
                            const bool cap_cuts_central_cell =
                                sphereSurfaceCutsAABB(
                                    p_cut,
                                    sphere_radius,
                                    central_min,
                                    central_max);

                            if (!cap_cuts_central_cell) continue;
                        } else {
                            // Non-tip truncated-cylinder case:
                            // Require theta_cut to be outside the central-cell + sphere-radius
                            // region. Equivalently, the cap sphere must not touch the central
                            // cell at all: not fully inside, not partially intersecting, not cutting.
                            const bool cap_touches_central_cell =
                                sphereSolidIntersectsAABB(
                                    p_cut,
                                    sphere_radius,
                                    central_min,
                                    central_max);

                            if (cap_touches_central_cell) continue;
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
                    // Analytical approximation with cylinders
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
                                area[i][j][k] = 0.0;
                            }
                    
                    // Build sphere

                    // raw 0th/1st/2nd moments for the sphere contribution
                    IRL::GeneralMoments3D<2> sphere_total_general_moments =
                        IRL::GeneralMoments3D<2>::fromScalarConstant(0.0);

                    std::vector<double> sphere_coarse_coords(stencil_size + 1); // Should be the same as coarse coords above
                    /* In the past I made a sphere + the cylinders, that led to a confusion with the sphere only class
                    generateSpecificSphere(
                        vfrac,
                        firstMoment,
                        area,
                        centroid,
                        surfaces,
                        sphere_origin,
                        sphere_radius,
                        sphere_coarse_coords,
                        nullptr   // do not compute centered 2nd moment here, this would require more effort to combine
                    );*/

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
                                double cell_area = 0.0;

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
                                    cell_area += 2 * M_PI * tube_radius * L;
                                }

                                if (cellV > 0.0) {
                                    double addV = cellV / cell_volume;
                                    area[i][j][k] += cell_area;

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

                    std::vector<std::vector<std::vector<double>>> areas_refined(refined_stencil_size,
                        std::vector<std::vector<double>>(refined_stencil_size,
                        std::vector<double>(refined_stencil_size)));

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

                                auto surface = volume_and_surface.getSurface();
                                areas_refined[i][j][k] = surface.getSurfaceArea();

                                if (secondMoment != nullptr) {
                                    auto gm = IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid);
                                    totalMoments_refined += gm;
                                }

                                if (visualize) {
                                    surfaces.push_back(surface);
                                }
                            }
                        }
                    }

                    // Compress refined to coarse stencil
                    compressStencilRefinedToCoarse(
                        volumes_refined,
                        firstMoments_refined,
                        areas_refined,
                        vfrac,
                        firstMoment,
                        area,
                        refinement_factor,
                        cell_volume,
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

        std::vector<float> generateState(bool well_resolved, int shape_type)
        {
            // Initialize information containers
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
            std::vector<std::vector<std::vector<double>>> area(
                stencil_size, 
                std::vector<std::vector<double>>(
                    stencil_size, 
                    std::vector<double>(
                        stencil_size
                    )
                )
            );
            Eigen::Matrix3d secondMoment;
            Eigen::Matrix3d* secondMomentPtr = nullptr;    // default: disabled using exact 2nd moment
            if (exact_2nd_moment && include_Moments >= 2) {
                std::cout<<"Enable exact 2nd moment calculation"<<std::endl;
                secondMomentPtr = &secondMoment; // This enables calculation of exact 2nd moments for all switch cases below
            }

            // Initialize characteristic lengths and set min bounds
            double min_cylinder_radius = machineZero;
            double min_sheet_thickness = machineZero;
            double min_sphere_radius = machineZero;
            double min_small_ellipsoid_axis = 0.05;
            
            double max_cylinder_radius = upper_limit_subgrid/2.0;
            double max_sheet_thickness = upper_limit_subgrid;
            double max_sphere_radius = upper_limit_subgrid/2.0;
            double max_small_ellipsoid_axis = upper_limit_subgrid/2.0;

            // Adjust bounds if well-resolved
            if (well_resolved) {
                min_cylinder_radius = max_cylinder_radius + machineZero;
                min_sheet_thickness = max_sheet_thickness + machineZero;
                min_sphere_radius = max_sphere_radius + machineZero;
                min_small_ellipsoid_axis = max_small_ellipsoid_axis + machineZero;

                max_cylinder_radius = class0_max_characteristic/2.0;
                max_sheet_thickness = class0_max_characteristic;
                max_sphere_radius = class0_max_characteristic/2.0;
                max_small_ellipsoid_axis = class0_max_characteristic/2.0;
            }

            double transition_sheet_min_subgrid_thickness = machineZero;
            double transition_sheet_max_subgrid_thickness = upper_limit_subgrid;
            double transition_sheet_well_resolved_thickness = class0_max_characteristic;

            // For sheets
            bool subgrid = !well_resolved;

            // Decide which shapes to generate
            std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
            double subshape_probability = prob_dist(eng);  // draw a random number in (0,1) to decide subtype of shape
            //std::cout<<"Generating datapoint of type "<<datapoint_type<<std::endl;
            switch (shape_type) {
                case 1: {
                    if (subshape_probability < 0.2) {
                        // 10% chance → truncated cylinder
                        if (subshape_probability < 0.1) {
                            generateBentTruncatedCylinder(vfrac, firstMoment, area, /*truncateinsidecentralcell*/false, min_cylinder_radius, max_cylinder_radius, secondMomentPtr);
                        }
                        // 10% chance → ellipsoid
                        else {
                            generateEllipsoidLigamentTip(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, /*truncateinsidecentralcell*/false, secondMomentPtr);
                        }
                    } else {
                        // 80% chance → continuous ligament
                        if (subshape_probability < 0.5) {
                            // 30% chance → ellipsoid ligament
                            generateEllipsoidLigament(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, secondMomentPtr);
                        } else {
                            //50% chance → bent cylinder
                            generateBentCylinder(vfrac, firstMoment, area, min_cylinder_radius, max_cylinder_radius, secondMomentPtr);
                        }
                    }
                    break;
                }
                case 2: {
                    generateEllipsoidDroplet(vfrac, firstMoment, area, min_small_ellipsoid_axis, max_small_ellipsoid_axis, secondMomentPtr);
                    break;
                }
                case 3: {
                    bool cut_inside_central_cell = false;
                    bool variable_sheet_thickness = false;

                    if (subshape_probability < 0.3) {
                        // 30% chance → sheet edge
                        if (subshape_probability < 0.1) {
                            // 10% chance → cut sheet with variable thickness
                            variable_sheet_thickness = true;
                            generateCutSheet(vfrac, firstMoment, area, /*cutinsidecentralcell*/ false, min_sheet_thickness, max_sheet_thickness, subgrid, variable_sheet_thickness, secondMomentPtr);
                        } else {
                            if (subshape_probability < 0.2) {
                                // 10% chance → sheet edge without variable thickness
                            variable_sheet_thickness = false;
                            generateCutSheet(vfrac, firstMoment, area, /*cutinsidecentralcell*/ false, min_sheet_thickness, max_sheet_thickness, subgrid, variable_sheet_thickness, secondMomentPtr);
                            } else {
                                // 10% chance → ellipsoid sheet edge
                                generateEllipsoidSheetEdge(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, cut_inside_central_cell, secondMomentPtr);
                            }
                        }
                    } else {
                        //70% chance → sheet
                        if (subshape_probability < 0.5) {
                            // 20% chance → sheet transition
                            bool thick_central_cell = false;
                            if (well_resolved) {
                                thick_central_cell = true; // for well-resolved, make central cell thick to 
                            }
                            generateSheetTransition(vfrac, firstMoment, area, thick_central_cell, transition_sheet_min_subgrid_thickness, transition_sheet_max_subgrid_thickness, transition_sheet_well_resolved_thickness, secondMomentPtr);
                        } else {
                            if (subshape_probability < 0.7) {
                                // 20% chance → sheet with variable thickness
                                variable_sheet_thickness = true;
                                generateSheet(vfrac, firstMoment, area, min_sheet_thickness,max_sheet_thickness, subgrid, variable_sheet_thickness, secondMomentPtr);
                            } else {
                                if (subshape_probability < 0.9) {
                                // 20% chance → normal sheet
                                variable_sheet_thickness = false;
                                generateSheet(vfrac, firstMoment, area, min_sheet_thickness,max_sheet_thickness, subgrid, variable_sheet_thickness, secondMomentPtr);
                                } else {
                                    // 10% chance → ellipsoid sheet
                                    generateEllipsoidSheet(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, secondMomentPtr);
                                }
                            }
                        }
                    }
                    break;
                }
                case 4: {
                    if (subshape_probability < 0.5) {
                        // 50% chance → ellipsoid ligament tip
                        generateEllipsoidLigamentTip(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, /*truncateinsidecentralcell*/false, secondMomentPtr);
                    } else {
                        // 50% chance → bent truncated cylinder
                        generateBentTruncatedCylinder(vfrac, firstMoment, area, /*truncateinsidecentralcell*/true, min_cylinder_radius, max_cylinder_radius, secondMomentPtr);
                    }
                    break;
                }
                case 5: {
                    if (subshape_probability < 0.5) {
                        // 50% chance → ellipsoid sheet edge
                        bool cut_inside_central_cell = true;
                        generateEllipsoidSheetEdge(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, cut_inside_central_cell, secondMomentPtr);
                    } else {
                        // 50% chance → cut sheet
                        bool variable_sheet_thickness = false;
                        if (subshape_probability < 0.75) {
                            // 25% chance → cut sheet with variable thickness
                            variable_sheet_thickness = true;
                            generateCutSheet(vfrac, firstMoment, area, /*cutinsidecentralcell*/ true, min_sheet_thickness, max_sheet_thickness, subgrid, variable_sheet_thickness, secondMomentPtr);
                        } else {
                            // 25% chance → cut sheet without variable thickness
                            variable_sheet_thickness = false;
                            generateCutSheet(vfrac, firstMoment, area, /*cutinsidecentralcell*/ true, min_sheet_thickness, max_sheet_thickness, subgrid, variable_sheet_thickness, secondMomentPtr);
                        }
                    }
                    break;
                }
                case 6:
                    generateParaboloid(vfrac, firstMoment, area, secondMomentPtr);
                    break;
                case 7: {
                    bool variable_sheet_thickness = false;
                    bool cut_inside_central_cell = true;
                    generateCutSheet(vfrac, firstMoment, area, cut_inside_central_cell, min_sheet_thickness, max_sheet_thickness, subgrid, variable_sheet_thickness, secondMomentPtr);
                    break; }
                case 8: {
                    bool variable_sheet_thickness = true; 
                    generateSheet(vfrac, firstMoment, area, min_sheet_thickness,max_sheet_thickness, subgrid, variable_sheet_thickness, secondMomentPtr);
                    break;}
                case 9:{
                    bool variable_sheet_thickness = false;
                    generateSheet(vfrac, firstMoment, area, min_sheet_thickness,max_sheet_thickness, subgrid, variable_sheet_thickness, secondMomentPtr);
                    break;}
                case 10:{
                    generateBentCylinder(vfrac, firstMoment, area, min_cylinder_radius, max_cylinder_radius, secondMomentPtr);
                    break;}
                case 11:
                    generateBentTruncatedCylinder(vfrac, firstMoment, area, /*truncateinsidecentralcell*/true, min_cylinder_radius, max_cylinder_radius, secondMomentPtr);
                    break;
                case 12: {
                    generateHyperbolicCylinder(vfrac, firstMoment, area, min_sheet_thickness, max_sheet_thickness, 45, secondMomentPtr);
                    break; }
                
                case 13: {
                    bool thick_central_cell = false;
                    generateSheetTransition(vfrac, firstMoment, area, thick_central_cell, transition_sheet_min_subgrid_thickness, transition_sheet_max_subgrid_thickness, transition_sheet_well_resolved_thickness, secondMomentPtr);
                    break;}
                case 14:
                    std::cout << "Generate Ellipsoid Droplet" << std::endl;
                    generateEllipsoidDroplet(vfrac, firstMoment, area, min_small_ellipsoid_axis, max_small_ellipsoid_axis, nullptr);
                    break;
                case 15:
                    std::cout << "Generate Ellipsoid Ligament" << std::endl;
                    generateEllipsoidLigament(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, nullptr);
                    break;
                case 16:
                    std::cout << "Generate Ellipsoid Ligament Tip" << std::endl;
                    generateEllipsoidLigamentTip(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, /*truncateinsidecentralcell*/true, nullptr);
                    break;
                case 17:
                    std::cout << "Generate Ellipsoid Sheet" << std::endl;
                    generateEllipsoidSheet(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, nullptr);
                    break;
                case 18:
                    std::cout << "Generate Ellipsoid Sheet Edge" << std::endl;
                    generateEllipsoidSheetEdge(vfrac, firstMoment, area, min_long_ellipsoid_axis, max_long_ellipsoid_axis, min_small_ellipsoid_axis, max_small_ellipsoid_axis, true, nullptr);
                    break;
                default:
                    // Paraboloid = case 0
                    generateParaboloid(vfrac, firstMoment, area, secondMomentPtr);
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
                        if (include_Surface_Area) {
                            flattened_state.push_back(area[i][j][k]);
                        }
                    }
                }
            }

            // make flattened_state a vector of floats
            std::vector<float> flattened_state_float(flattened_state.begin(), flattened_state.end());

            // Below: Append 2nd moments and eigenvalues. Now moved to preprocessing. This removes the ability to use accurate 2nd moments, but this has not shown significant improvement yet.
            /*
            if (include_Moments >= 2) {
                if (exact_2nd_moment) {
                    // Append exact 2nd moments to flattened_state
                    flattened_state_float.push_back(static_cast<float>(secondMoment(0, 0))); // Ixx
                    flattened_state_float.push_back(static_cast<float>(secondMoment(1, 1))); // Iyy
                    flattened_state_float.push_back(static_cast<float>(secondMoment(2, 2))); // Izz
                    flattened_state_float.push_back(static_cast<float>(secondMoment(0, 1))); // Ixy
                    flattened_state_float.push_back(static_cast<float>(secondMoment(0, 2))); // Ixz
                    flattened_state_float.push_back(static_cast<float>(secondMoment(1, 2))); // Iyz
                } else {
                    // Append approximate 2nd moments to flattened_state
                    Eigen::Matrix3d approxSecondMoment = IRL::compute2ndMoment(flattened_state_float, stencil_size, 1, include_Surface_Area);
                    flattened_state_float.push_back(static_cast<float>(approxSecondMoment(0, 0))); // Ixx
                    flattened_state_float.push_back(static_cast<float>(approxSecondMoment(1, 1))); // Iyy
                    flattened_state_float.push_back(static_cast<float>(approxSecondMoment(2, 2))); // Izz
                    flattened_state_float.push_back(static_cast<float>(approxSecondMoment(0, 1))); // Ixy
                    flattened_state_float.push_back(static_cast<float>(approxSecondMoment(0, 2))); // Ixz
                    flattened_state_float.push_back(static_cast<float>(approxSecondMoment(1, 2))); // Iyz
                }
            }

            if (include_Eigenvalues) {
                IRL::appendInertiaEigenvalues(flattened_state_float, stencil_size, include_Moments, 1, include_Surface_Area, machineZero);
            }
            */

            return flattened_state_float;
        }

        void generateData (std::vector<std::vector<float>>* statesV, std::vector<int>* labelsV)
        {
            int subgrid_subclasses = 5;

            std::cout << "Generating " << no_datapoints << " datapoints." << std::endl;

            // Initialize random number generator with a seed from current time
            std::srand(std::time(0));
            std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

            for (int i = 0; i < no_datapoints; i++) {
                bool well_resolved = true;
                int shape_type = 0;
                int class_type = 0;
                // First decide if this datapoint is well-resolved or subgrid
                double probability_well_resolved = prob_dist(eng);
                if (probability_well_resolved > (static_cast<double>(subgrid_subclasses + 1) / static_cast<double>(2 * subgrid_subclasses + 1))) {
                    well_resolved = true;
                    class_type = 0;
                    std::uniform_int_distribution<int> subclass_dist(0, subgrid_subclasses);
                    shape_type = subclass_dist(eng);
                } else {
                    well_resolved = false;
                    std::uniform_int_distribution<int> subclass_dist(1, subgrid_subclasses);
                    class_type = subclass_dist(eng);
                    shape_type = class_type;
                }

                labelsV->push_back(class_type);
                statesV->push_back(generateState(
                    well_resolved,
                    shape_type
                ));
                if (i % 10000 == 0) {
                    std::cout << "Generated " << i << " datapoints." << std::endl;
                }
            }
        }
    };
}