#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>

#include "irl/generic_cutting/generic_cutting.h"

namespace IRL {
    class Data_gen {

        public:
        std::mt19937_64 eng;

        Data_gen(){
            // Initialize random number generator with a seed from current time
            std::random_device rd;
            eng = std::mt19937_64(rd());  // Use random_device for seeding
            std::cout << "I'm a data generator!" << std::endl;
        }

        void generate_Paraboloid (double*** vfrac, int stencil_size, int datapoint_type, int plane_bounds_coefficients) {

            // Defining cell coordinates
            auto coords = std::vector<double>(stencil_size+1);
            for (int i = 0; i < stencil_size+1; i++) {
                coords[i] = -0.5*stencil_size + static_cast<double>(i);
            }
            const double cell_volume = 1.0;

            std::uniform_real_distribution<double> random_rotation(-0.5 * M_PI, 0.5 * M_PI);

            std::uniform_real_distribution<double> random_translation(-0.5, 0.5);
            
            // Bounds of paraboloid parameters, 0=plane, 1=paraboloid from int plane_bounds_coefficients
            // Declare random_coeff outside of the if-else block so it's accessible in both branches
            std::uniform_real_distribution<double> random_coeff(0.0, 0.0); // Declare default (0, 0) range
            if (datapoint_type == 0) {
                random_coeff = std::uniform_real_distribution<double>(0, plane_bounds_coefficients); // Modify range
            } else {
                random_coeff = std::uniform_real_distribution<double>(plane_bounds_coefficients, 5.0); // Modify range
            }
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
                            AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(
                            cell, paraboloid);
                        // Store surface and volume (fraction)
                        surfaces.push_back(volume_and_surface.getSurface());
                        vfrac[i][j][k] = volume_and_surface.getMoments().volume() / cell_volume;
                        //std::cout<<vfrac[i][j][k]<<"\n";
                    }
                }
            }
        }

        // Function to generate a random vector within given bounds
        Eigen::Vector3d generateRandomPoint(double min_bound, double max_bound, std::mt19937_64& eng) {
            std::uniform_real_distribution<double> dist(min_bound, max_bound);
            return Eigen::Vector3d(dist(eng), dist(eng), dist(eng));  // Random point in 3D space
        }

        void generate_Cylinder (double*** vfrac, int stencil_size) {

            const double cell_volume = 1.0;
            int refinement_factor = 3;

            // Define bounds for the random points (you can change these values)
            double min_bound = -refinement_factor;  // Minimum bound for each coordinate
            double max_bound = refinement_factor;   // Maximum bound for each coordinate

            Eigen::Vector3d axis_point1 = generateRandomPoint(min_bound, max_bound, eng);
            Eigen::Vector3d axis_point2 = generateRandomPoint(min_bound, max_bound, eng);

            Eigen::Vector3d axis_direction = axis_point2-axis_point1;

            std::uniform_real_distribution<double> dist(0.5, stencil_size/2);
            double radius =  dist(eng);

            // Create refined cells
            int refined_stencil_size = refinement_factor*stencil_size;
            std::vector<std::vector<std::vector<double>>> refined_vfrac(stencil_size, std::vector<std::vector<double>>(refined_stencil_size, std::vector<double>(refined_stencil_size)));
            // Defining cell coordinates
            auto coords = std::vector<double>(refined_stencil_size+1);
            for (int i = 0; i < refined_stencil_size+1; i++) {
                coords[i] = -0.5*refined_stencil_size + static_cast<double>(i);
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
                            IRL::AddSurfaceOutput<IRL::Volume, IRL::ParametrizedSurfaceOutput>>(
                            cell, paraboloid);
                        // Store surface and volume (fraction)
                        
                        refined_vfrac[i][j][k] = volume_and_surface.getMoments().volume() / (cell_volume/(refinement_factor*refinement_factor*refinement_factor));
                    }
                }
            }

            // Compress refined mesh into original stencil size
            for (int i = 0; i < stencil_size; i++) {
                for (int j = 0; j < stencil_size; j++) {
                    for (int k = 0; k < stencil_size; k++) {
                        double sum = 0.0;

                        for (int m = i*refinement_factor; m < (i+1)*refinement_factor; m++) {
                            for (int n = j*refinement_factor; n < (j+1)*refinement_factor; n++) {
                                for (int o = k*refinement_factor; o < (k+1)*refinement_factor; o++) {
                                    sum += refined_vfrac[m][n][o];
                                }
                            }
                        }

                        vfrac [i][j][k] = sum;
                    }
                }
            }

            // Debug:
            bool valid_size = true;
            // Check if the first level of pointers (vfrac) is valid
            if (vfrac == nullptr) {
                valid_size = false;
                std::cout << "The 3D array has not been allocated properly (1st dimension)." << std::endl;
            }

            // Iterate over the first dimension (vfrac[i]) to check if it's properly allocated
            for (int i = 0; i < stencil_size; ++i) {
                if (vfrac[i] == nullptr) {
                    valid_size = false;
                    std::cout << "The 3D array has not been allocated properly (2nd dimension) at index " << i << "." << std::endl;
                    break; // Exit if a nullptr is found in the second dimension
                }

                // Iterate over the second dimension (vfrac[i][j]) to check if it's properly allocated
                for (int j = 0; j < stencil_size; ++j) {
                    if (vfrac[i][j] == nullptr) {
                        valid_size = false;
                        std::cout << "The 3D array has not been allocated properly (3rd dimension) at indices [" << i << "][" << j << "]." << std::endl;
                        break; // Exit if a nullptr is found in the third dimension
                    }

                    // Now check the third dimension (vfrac[i][j][k]) to make sure it points to valid double values
                    for (int k = 0; k < stencil_size; ++k) {
                        if (std::isnan(vfrac[i][j][k])) {
                            valid_size = false;
                            break; // Exit if an invalid value is found at [i][j][k]
                        }
                    }
                }
            }

            // If size is invalid, print an error message
            if (!valid_size) {
                std::cout << "The 3D array does not have the desired size or was not properly allocated." << std::endl;
            }
        }

        std::vector<double> generate_State(int datapoint_type, int stencil_size, double*** vfrac, double plane_bounds_coefficients){
            // Create a vector of volume fractions and surfaces
            // OLD std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));

            if (datapoint_type == 2) {
                generate_Cylinder(vfrac, stencil_size);
            } else {
                generate_Paraboloid(vfrac, stencil_size, datapoint_type, plane_bounds_coefficients);
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

        void generate_Data (std::vector<std::vector<double>>* statesV, std::vector<int>* labelsV, int no_datapoints, int stencil_size = 3, double plane_bounds_coefficients = 0.5){
            std::cout << no_datapoints << std::endl;
            // Initialize random number generator with a seed from current time
            std::srand(std::time(0));
            
            // Make a pointer, 3D array of volume fractions
            double*** vfrac = new double**[stencil_size];
            for (int i = 0; i < stencil_size; ++i) {
                vfrac[i] = new double*[stencil_size];
                for (int j = 0; j < stencil_size; ++j) {
                    vfrac[i][j] = new double[stencil_size];
                }
            }

            for (int i=0; i<no_datapoints; i++) {
                // Generate the data, init with a random number 0 or 1, 0=plane, 1=paraboloid, 2=cylinder
                int datapoint_type = std::rand() % 3;
                labelsV->push_back(datapoint_type);
                statesV->push_back(generate_State(datapoint_type, stencil_size, vfrac, plane_bounds_coefficients));
            }

            // Free the memory
            for (int i = 0; i < stencil_size; ++i) {
                for (int j = 0; j < stencil_size; ++j) {
                    delete[] vfrac[i][j];  // Free each row of doubles
                }
                delete[] vfrac[i];  // Free each 2D array
            }
            delete[] vfrac;  // Free the top-level array
        }
    };
}