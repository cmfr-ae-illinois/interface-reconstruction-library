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

        std::vector<double> generate(int plane_or_paraboloid_in, int stencil_size,  double plane_bounds_coefficients){
            // Create a vector of volume fractions and surfaces
            std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));

            // Defining cell coordinates
            auto coords = std::vector<double>(stencil_size+1);
            for (int i = 0; i < stencil_size+1; i++) {
                coords[i] = -0.5*stencil_size + static_cast<double>(i);
            }
            const double cell_volume = 1.0;

            std::uniform_real_distribution<double> random_rotation(-0.5 * M_PI,
                                                                    0.5 * M_PI);

            std::uniform_real_distribution<double> random_translation(-0.5, 0.5);
            
            // Bounds of paraboloid parameters, 0=plane, 1=paraboloid from int plane_bounds_coefficients
            // Declare random_coeff outside of the if-else block so it's accessible in both branches
            std::uniform_real_distribution<double> random_coeff(0.0, 0.0); // Declare default (0, 0) range
            if (plane_or_paraboloid_in == 0) {
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
            // Flatten the 3D vector vfrac into a 1D vector
            std::vector<double> flattened_vfrac;
            for (const auto& row : vfrac) {
                for (const auto& col : row) {
                    flattened_vfrac.insert(flattened_vfrac.end(), col.begin(), col.end());
                }
            }
            //return the result
            return flattened_vfrac;
        }

        void generate_Data (std::vector<std::vector<double>>* statesV, std::vector<int>* labelsV, int no_paraboloids, int stencil_size = 3, double plane_bounds_coefficients = 0.5){
            std::cout << no_paraboloids << std::endl;
            // Initialize random number generator with a seed from current time
            std::srand(std::time(0));
            for (int i=0; i<no_paraboloids; i++) {
                // Generate the data, init with a random number 0 or 1, 0=plane, 1=paraboloid
                int plane_or_paraboloid = std::rand() % 2;
                labelsV->push_back(plane_or_paraboloid);
                statesV->push_back(generate(plane_or_paraboloid, stencil_size, plane_bounds_coefficients));
            }
        }
    };
}