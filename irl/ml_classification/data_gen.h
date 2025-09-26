#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/ml_classification/vtk_out.h"

namespace IRL {
    class Data_gen {

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

        std::vector<double> generate_Paraboloid (
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, bool visualize = false) 
            {
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
                coords[i] = -0.5*stencil_size + static_cast<double>(i);
            }
            const double cell_volume = 1.0;

            std::uniform_real_distribution<double> random_rotation(-0.5 * M_PI, 0.5 * M_PI);

            std::uniform_real_distribution<double> random_translation(-0.5, 0.5);

            std::uniform_real_distribution<double> random_coeff = std::uniform_real_distribution<double>(0, 5.0); // Modify range
            
            // Define a distribution for picking the sign (either -1 or 1)
            std::uniform_int_distribution<int> signPick(0, 1);

            // Create reference frame
            const auto frame = ReferenceFrame(
                Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
            // Create random rotation angles, datum and coefficiens
            double angles[3] = {random_rotation(eng), random_rotation(eng), 0.0};
            Pt datum(random_translation(eng), random_translation(eng),
                            random_translation(eng));
            double coeff1 = random_coeff(eng) * (signPick(eng) == 0 ? -1 : 1);
            double coeff2 = random_coeff(eng) * (signPick(eng) == 0 ? -1 : 1);
            // Generate rotation matrices
            UnitQuaternion x_rot(angles[0], frame[0]);
            UnitQuaternion y_rot(angles[1], frame[1]);
            UnitQuaternion z_rot(angles[2], frame[2]);
            // Create random paraboloid
            const auto paraboloid = Paraboloid(datum, x_rot * y_rot * z_rot * frame,
                                                    coeff1, coeff2);

            // Initialize field
            std::vector<ParametrizedSurfaceOutput> surfaces;
            for (int i = 0; i < stencil_size; i++) {
                for (int j = 0; j < stencil_size; j++) {
                    for (int k = 0; k < stencil_size; k++) {
                        // Create cell
                        auto cell = RectangularCuboid::fromBoundingPts(
                            Pt(coords[i], coords[j], coords[k]),
                            Pt(coords[i + 1], coords[j + 1], coords[k + 1]));
                        // Intersect cell with paraboloid -- return volume and surface
                        auto volume_fraction = getVolumeFraction(cell, paraboloid);
                        //std::cout << "VFRAC(" << i << ", " << j << ", " << k
                                //<< ") = " << volume_fraction << std::endl;
                        auto volume_and_surface = getVolumeMoments<
                            //AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(
                            AddSurfaceOutput<VolumeMoments, ParametrizedSurfaceOutput>>(
                            cell, paraboloid);
                        // Store surface and volume (fraction)
                        surfaces.push_back(volume_and_surface.getSurface());
                        vfrac[i][j][k] = volume_and_surface.getMoments().volume() / cell_volume;
                        firstMoment[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                             volume_and_surface.getMoments().centroid().y(),
                                             volume_and_surface.getMoments().centroid().z();

                        if (visualize) {
                            centroid[i][j][k] = computeCentroidFromFirstMoment(firstMoment[i][j][k], vfrac[i][j][k] * cell_volume);
                        }
                    }
                }
            }

            // Generate vtk output
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
            //return the result
            return flattened_vfrac;
        }

        std::vector<double> generate_Sheet (
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, bool visualize = false) {

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
                coords[i] = -0.5*stencil_size + static_cast<double>(i);
            }
            const double cell_volume = 1.0;

            Eigen::Vector3d datum = generateRandomPoint(-1.0 , 1.0, eng);

            Eigen::Vector3d direction = generateRandomPoint(-1.0 , 1.0, eng); // direction in chich the paraboloid will open. Will be used to calc datum offset and z-direction of frame as well.

            std::uniform_real_distribution<double> random_thickness(0.1, 0.5);
            
            // OLD: Uniform distribution for coeffs
            //std::uniform_real_distribution<double> random_coeff(0.0, 5.0);

            // New: Normal distribution for coeffs, mean 0, stddev 0.75
            std::normal_distribution<double> random_coeff(0.0, 0.75);

            // Define a distribution for picking the sign (either -1 or 1) for coeffs
            std::uniform_int_distribution<int> signPick(0, 1);

            double thickness = random_thickness(eng);

            //Eigen::Vector3d datum_paraboloid1_eVec = (direction-datum).normalized()*(0.001);
            //Eigen::Vector3d datum_paraboloid2_eVec = (direction-datum).normalized()*(0.001+thickness);
            Eigen::Vector3d datum_paraboloid1_eVec = (direction.normalized()*(-thickness/2.0)-datum);
            Eigen::Vector3d datum_paraboloid2_eVec = (direction.normalized()*(thickness/2.0)-datum);

            // Create frame
            Eigen::Vector3d paraboloid_z = direction;
            paraboloid_z.normalize();
            Eigen::Vector3d helper = generateRandomPoint(-1.0 , 1.0, eng); //Random point, TODO: MAKE SURE IT IS NOT PARALLEL, is sufficient for random rotation???
            Eigen::Vector3d paraboloid_x = paraboloid_z.cross(helper);
            paraboloid_x.normalize();
            Eigen::Vector3d paraboloid_y = paraboloid_z.cross(paraboloid_x);
            paraboloid_y.normalize();
            //frame[1] = IRL::crossProduct(frame[2],frame[0]); frame2 is normal, frame0 is axis
            //frame[1].normalize();

            const auto frame = IRL::ReferenceFrame(IRL::Normal(paraboloid_x.x(), paraboloid_x.y(), paraboloid_x.z()), 
                                    IRL::Normal(paraboloid_y.x(), paraboloid_y.y(), paraboloid_y.z()), 
                                    IRL::Normal(paraboloid_z.x(), paraboloid_z.y(), paraboloid_z.z()));

            //Pt datum1(paraboloid_x.x(), paraboloid_y.y(), paraboloid_z.z() - thickness/2.0);
            //Pt datum2(paraboloid_x.x(), paraboloid_y.y(), paraboloid_z.z() + thickness/2.0);

            //Eigen::Vector3d datum_paraboloid1_eVec = (paraboloid_z-datum).normalized()*(thickness/2.0);
            //Eigen::Vector3d datum_paraboloid2_eVec = (paraboloid_z-datum).normalized()*(-thickness/2.0);

            //Eigen::Vector3d datum_paraboloid1_eVec = (datum-paraboloid_z).normalized()*(thickness/2.0);
            //Eigen::Vector3d datum_paraboloid2_eVec = (datum-paraboloid_z).normalized()*(-thickness/2.0);

            IRL::Pt datum_paraboloid1(datum_paraboloid1_eVec.x(), datum_paraboloid1_eVec.y(), datum_paraboloid1_eVec.z());
            IRL::Pt datum_paraboloid2(datum_paraboloid2_eVec.x(), datum_paraboloid2_eVec.y(), datum_paraboloid2_eVec.z());

            // OLD: Uniform distribution for coeffs
            //double coeff1 = random_coeff(eng) * (signPick(eng) == 0 ? -1 : 1);
            //double coeff2 = random_coeff(eng) * (signPick(eng) == 0 ? -1 : 1);

            // New: Normal distribution for coeffs, mean 0, stddev 1
            //double coeff1 = random_coeff(eng);
            //double coeff2 = random_coeff(eng);

            // now with clamp
            const double max_abs_coeff = 5.0;

            auto sample_bounded_normal = [&](auto& rng, auto& dist, double bound) -> double {
                double value;
                do {
                    value = dist(rng);
                } while (std::abs(value) > bound);
                return value;
            };

            double coeff1 = sample_bounded_normal(eng, random_coeff, max_abs_coeff);
            double coeff2 = sample_bounded_normal(eng, random_coeff, max_abs_coeff);

            // Create random paraboloid

            const auto paraboloid1 = Paraboloid(datum_paraboloid1, frame, coeff1, coeff2);
            const auto paraboloid2 = Paraboloid(datum_paraboloid2, frame, coeff1, coeff2);


            // Initialize field
            std::vector<ParametrizedSurfaceOutput> surfaces;
            for (int i = 0; i < stencil_size; i++) {
                for (int j = 0; j < stencil_size; j++) {
                    for (int k = 0; k < stencil_size; k++) {
                        // Create cell
                        auto cell = RectangularCuboid::fromBoundingPts(
                            Pt(coords[i], coords[j], coords[k]),
                            Pt(coords[i + 1], coords[j + 1], coords[k + 1]));

                        // Intersect cell with paraboloid -- return volume and surface

                        auto volume_fraction1 = getVolumeFraction(cell, paraboloid1);
                        auto volume_fraction2 = getVolumeFraction(cell, paraboloid2);

                        auto volume_and_surface1 = getVolumeMoments<
                            AddSurfaceOutput<VolumeMoments, ParametrizedSurfaceOutput>>(
                            cell, paraboloid1);
                        auto volume_and_surface2 = getVolumeMoments<
                            AddSurfaceOutput<VolumeMoments, ParametrizedSurfaceOutput>>(
                            cell, paraboloid2);

                        // Store surface and volume (fraction)
                        surfaces.push_back(volume_and_surface1.getSurface());
                        surfaces.push_back(volume_and_surface2.getSurface());

                        double V1 = volume_and_surface1.getMoments().volume();
                        double V2 = volume_and_surface2.getMoments().volume();
                        double Vdiff = V2 - V1;

                        //vfrac[i][j][k] = (volume_and_surface2.getMoments().volume() - volume_and_surface1.getMoments().volume()) / cell_volume;

                        // Exclude negative machine zero volumes so that machine learning does not train on them
                        if (Vdiff < 0.0) {
                            Vdiff = 0.0;
                        }

                        vfrac[i][j][k] = Vdiff / cell_volume;
                        
                        Eigen::Vector3d M1(volume_and_surface1.getMoments().centroid().x(),
                                        volume_and_surface1.getMoments().centroid().y(),
                                        volume_and_surface1.getMoments().centroid().z());

                        Eigen::Vector3d M2(volume_and_surface2.getMoments().centroid().x(),
                                        volume_and_surface2.getMoments().centroid().y(),
                                        volume_and_surface2.getMoments().centroid().z());

                        firstMoment[i][j][k] = M2 - M1;

                        if (visualize) {
                            centroid[i][j][k] = computeCentroidFromFirstMoment(firstMoment[i][j][k], vfrac[i][j][k] * cell_volume);
                        }
                    }
                }
            }

            // Generate vtk output
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
            //return the result
            return flattened_vfrac;
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

                        // Store normalized volume fraction
                        vfrac[i][j][k] = vol_sum / compressed_cell_volume;

                        if (vol_sum > 1e-14) {
                            firstMoment[i][j][k] = firstMoment_sum;
                        } else {
                            firstMoment[i][j][k] = Eigen::Vector3d::Zero();
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

        std::vector<double> generate_Cylinder (
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, bool visualize = false) {

            // make centroid, only used for visualization
            std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                stencil_size,
                std::vector<std::vector<Eigen::Vector3d>>(
                    stencil_size,
                    std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                )
            );

            // for visualization option:
            std::vector<IRL::ParametrizedSurfaceOutput> surfaces;

            // define parameters
            //std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));
            const double cell_volume = 1.0;
            int refinement_factor = 3;
            double refinement_factor_double = static_cast<double>(refinement_factor);

            // Define bounds for the random points (you can change these values)
            double min_bound = -refinement_factor_double;  // Minimum bound for each coordinate
            double max_bound = refinement_factor_double;   // Maximum bound for each coordinate

            Eigen::Vector3d axis_point1 = generateRandomPoint(min_bound, max_bound, eng);
            Eigen::Vector3d axis_point2 = generateRandomPoint(min_bound, max_bound, eng);

            Eigen::Vector3d axis_direction = axis_point2-axis_point1;

            std::uniform_real_distribution<double> dist(0.5, stencil_size/2);
            double radius =  dist(eng);

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
                        auto volume_and_surface = IRL::getVolumeMoments<
                            IRL::AddSurfaceOutput<IRL::VolumeMoments, IRL::ParametrizedSurfaceOutput>>(
                            cell, paraboloid);
                        // Store surface and volume (fraction)
                        
                        // use volumes directly instead of fractions
                        volumes_refined[i][j][k] = volume_and_surface.getMoments().volume();
                                         //volume_and_surface.getMoments().centroid() / vol

                        firstMoments_refined[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                        volume_and_surface.getMoments().centroid().y(),
                                                        volume_and_surface.getMoments().centroid().z();

                        if (visualize) {
                            surfaces.push_back(volume_and_surface.getSurface());
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

            // Generate vtk output
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

        std::vector<double> generate_Sphere (
            std::vector<std::vector<std::vector<double>>>& vfrac,
            std::vector<std::vector<std::vector<Eigen::Vector3d>>>& firstMoment,
            int stencil_size, bool visualize = false) {

            // make centroid, only used for visualization
            std::vector<std::vector<std::vector<Eigen::Vector3d>>> centroid(
                stencil_size,
                std::vector<std::vector<Eigen::Vector3d>>(
                    stencil_size,
                    std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                )
            );

            // for visualization option:
            std::vector<IRL::ParametrizedSurfaceOutput> surfaces;

            // define parameters
            //std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));
            const double cell_volume = 1.0;
            int refinement_factor = 3;
            double refinement_factor_double = static_cast<double>(refinement_factor);
            std::uniform_real_distribution<double> dist(0.1, 1); //The minimum size of the sphere is set here

            double radius =  dist(eng);

            Eigen::Vector3d origin = generateRandomPoint(-1.0 + radius , 1.0 - radius, eng);

            // Create refined cells
            int refined_stencil_size = refinement_factor*stencil_size;
            std::vector<std::vector<std::vector<double>>> volumes_refined(refined_stencil_size,
                std::vector<std::vector<double>>(refined_stencil_size,
                std::vector<double>(refined_stencil_size)));

            std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoments_refined(refined_stencil_size,
                std::vector<std::vector<Eigen::Vector3d>>(refined_stencil_size,
                std::vector<Eigen::Vector3d>(refined_stencil_size, Eigen::Vector3d::Zero())));

            //std::vector<std::vector<std::vector<double>>> volumes(refined_stencil_size, std::vector<std::vector<double>>(refined_stencil_size, std::vector<double>(refined_stencil_size)));

            // Use refinement factor on parameters -> the grid is parted, so the parameters have to be scaled to be the same
            radius = radius * refinement_factor_double;
            origin = origin * refinement_factor_double;

            // Defining cell coordinates
            auto coords = std::vector<double>(refined_stencil_size+1);
            for (int i = 0; i <= refined_stencil_size; i++) {
                coords[i] = -0.5 * stencil_size + (static_cast<double>(i) / refinement_factor);
            }

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
                        auto volume_fraction = IRL::getVolumeFraction(cell, paraboloid);
                        auto volume_and_surface = IRL::getVolumeMoments<
                            IRL::AddSurfaceOutput<IRL::VolumeMoments, IRL::ParametrizedSurfaceOutput>>(
                            cell, paraboloid);
                        // Store surface and volume (fraction)
                        
                        volumes_refined[i][j][k] = volume_and_surface.getMoments().volume();
                                         //volume_and_surface.getMoments().centroid() / vol

                        Eigen::Vector3d m1(volume_and_surface.getMoments().centroid().x(),
                                        volume_and_surface.getMoments().centroid().y(),
                                        volume_and_surface.getMoments().centroid().z());

                        firstMoments_refined[i][j][k] << volume_and_surface.getMoments().centroid().x(),
                                                        volume_and_surface.getMoments().centroid().y(),
                                                        volume_and_surface.getMoments().centroid().z();

                        if (visualize) {
                            surfaces.push_back(volume_and_surface.getSurface());
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

            // Generate vtk output
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


        std::vector<double> generate_State(int datapoint_type, int stencil_size, bool include_firstMoment, bool visualize = false) {
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
            /*
            std::vector<std::vector<std::vector<std::vector<double>>>> firstMoment(
                stencil_size,
                std::vector<std::vector<std::vector<double>>>(
                    stencil_size,
                    std::vector<std::vector<double>>(
                        stencil_size,
                        std::vector<double>(3)  // each first moment has 3 components
                    )
                )
            );
            */

            std::vector<std::vector<std::vector<Eigen::Vector3d>>> firstMoment(
                stencil_size,
                std::vector<std::vector<Eigen::Vector3d>>(
                    stencil_size,
                    std::vector<Eigen::Vector3d>(stencil_size, Eigen::Vector3d::Zero())
                )
            );

            switch (datapoint_type) {
                case 1:
                    generate_Cylinder(vfrac, firstMoment, stencil_size, visualize);
                    break;
                case 2:
                    generate_Sphere(vfrac, firstMoment, stencil_size, visualize);
                    break;
                case 3:
                    generate_Sheet(vfrac, firstMoment, stencil_size, visualize);
                    break;
                default:
                    generate_Paraboloid(vfrac, firstMoment, stencil_size, visualize);
                    break;
            }

            // Flatten the 3D vectors into one 1D vector
            std::vector<double> flattened_state;
            for (int i = 0; i < stencil_size; ++i) {
                for (int j = 0; j < stencil_size; ++j) {
                    for (int k = 0; k < stencil_size; ++k) {
                        flattened_state.push_back(vfrac[i][j][k]);
                        if (include_firstMoment) {
                            flattened_state.push_back(firstMoment[i][j][k].x());
                            flattened_state.push_back(firstMoment[i][j][k].y());
                            flattened_state.push_back(firstMoment[i][j][k].z());
                        }
                    }
                }
            }
            return flattened_state;
        }

        void generate_Data (std::vector<std::vector<double>>* statesV, std::vector<int>* labelsV, int no_datapoints, int stencil_size = 3, int no_datapoint_types_in = 4, bool include_centroid = false){
            std::cout << no_datapoints << std::endl;
            // Initialize random number generator with a seed from current time
            std::srand(std::time(0));

            for (int i=0; i<no_datapoints; i++) {
                // Generate the data, init with a random number 0 or 1, 0=paraboloid, 1=cylinder, 2=sphere, 3=sheet

                int datapoint_type = std::rand() % no_datapoint_types_in;

                labelsV->push_back(datapoint_type);
                statesV->push_back(generate_State(datapoint_type, stencil_size, include_centroid));
            }
        }
    };
}