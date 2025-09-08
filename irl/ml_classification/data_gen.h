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

        std::vector<double> generate_Paraboloid (int stencil_size, int datapoint_type, int plane_bounds_coefficients, bool visualize = false) {
            
            std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));
  
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

            // Can be used to differantiate planes and paraboloids:
            
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

            // Generate vtk output
            if (visualize) {
                WriteField(stencil_size, coords, vfrac, "vfrac");
                WriteSurface(surfaces, "surface");
                std::cout << "VFRAC:" << std::endl;
                for (int i = 0; i < vfrac.size(); i++) {
                    for (int j = 0; j < vfrac[i].size(); j++) {
                        for (int k = 0; k < vfrac[i][j].size(); k++) {
                            std::cout << "[" << i << "," << j << "," << k << "] " << std::fixed << std::setprecision(3) << vfrac [i][j][k] << "   ";
                        }   
                    std::cout << std::endl;
                    }
                }
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

        std::vector<double> generate_Sheet (int stencil_size, bool visualize = false) {
            
            std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));
  
            // Defining cell coordinates
            auto coords = std::vector<double>(stencil_size+1);
            for (int i = 0; i < stencil_size+1; i++) {
                coords[i] = -0.5*stencil_size + static_cast<double>(i);
            }
            const double cell_volume = 1.0;

            Eigen::Vector3d datum = generateRandomPoint(-1.0 , 1.0, eng);

            Eigen::Vector3d direction = generateRandomPoint(-1.0 , 1.0, eng); // direction in chich the paraboloid will open. Will be used to calc datum offset and z-direction of frame as well.

            std::uniform_real_distribution<double> random_thickness(0.1, 1);
            
            std::uniform_real_distribution<double> random_coeff(0.0, 5.0);

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

            double coeff1 = random_coeff(eng) * (signPick(eng) == 0 ? -1 : 1);
            double coeff2 = random_coeff(eng) * (signPick(eng) == 0 ? -1 : 1);

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
                            AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(
                            cell, paraboloid1);
                        auto volume_and_surface2 = getVolumeMoments<
                            AddSurfaceOutput<Volume, ParametrizedSurfaceOutput>>(
                            cell, paraboloid2);

                        // Store surface and volume (fraction)
                        surfaces.push_back(volume_and_surface1.getSurface());
                        surfaces.push_back(volume_and_surface2.getSurface());

                        vfrac[i][j][k] = (volume_and_surface2.getMoments().volume() - volume_and_surface1.getMoments().volume()) / cell_volume;

                        // Exclude negative machine zero volumes so that machine learning does not train on them
                        if (vfrac[i][j][k] < 0.0) {
                            vfrac[i][j][k] = 0.0;
                        }
                        //vfrac[i][j][k] -= volume_and_surface1.getMoments().volume() / cell_volume;
                    }
                }
            }

            // Generate vtk output
            if (visualize) {
                WriteField(stencil_size, coords, vfrac, "vfrac");
                WriteSurface(surfaces, "surface");
                std::cout << "VFRAC:" << std::endl;
                for (int i = 0; i < vfrac.size(); i++) {
                    for (int j = 0; j < vfrac[i].size(); j++) {
                        for (int k = 0; k < vfrac[i][j].size(); k++) {
                            std::cout << "[" << i << "," << j << "," << k << "] " << std::fixed << std::setprecision(3) << vfrac [i][j][k] << "   ";
                        }   
                    std::cout << std::endl;
                    std::cout << thickness << std::endl;
                    }
                }
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

        std::vector<double> generate_Cylinder (int stencil_size, bool visualize = false) {

            // for visualization option:
            std::vector<IRL::ParametrizedSurfaceOutput> surfaces;

            // define parameters
            std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));
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
            std::vector<std::vector<std::vector<double>>> volumes(refined_stencil_size, std::vector<std::vector<double>>(refined_stencil_size, std::vector<double>(refined_stencil_size)));
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
                            IRL::AddSurfaceOutput<IRL::VolumeMoments, IRL::ParametrizedSurfaceOutput>>(
                            cell, paraboloid);
                        // Store surface and volume (fraction)
                        
                        // use volumes directly instead of fractions
                        volumes[i][j][k] = volume_and_surface.getMoments().volume();
                                         //volume_and_surface.getMoments().centroid() / vol
                        if (visualize) {
                            surfaces.push_back(volume_and_surface.getSurface());
                        }
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
                                    sum += volumes[m][n][o]; // Volumes are between 0 and 1
                                    
                                }
                            }
                        }
                        vfrac [i][j][k] = sum / (cell_volume*(refinement_factor_double*refinement_factor_double*refinement_factor_double));
                    }
                }
            }

            // Generate vtk output
            if (visualize) {
                WriteField(stencil_size, coords, vfrac, "vfrac");
                WriteSurface(surfaces, "surface");
                std::cout << "VFRAC:" << std::endl;
                for (int i = 0; i < vfrac.size(); i++) {
                    for (int j = 0; j < vfrac[i].size(); j++) {
                        for (int k = 0; k < vfrac[i][j].size(); k++) {
                            std::cout << "[" << i << "," << j << "," << k << "] " << std::fixed << std::setprecision(3) << vfrac [i][j][k] << "   ";
                        }   
                    std::cout << std::endl;
                    }
                }
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

        std::vector<double> generate_Sphere (int stencil_size, bool visualize = false) {

            // for visualization option:
            std::vector<IRL::ParametrizedSurfaceOutput> surfaces;

            // define parameters
            std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));
            const double cell_volume = 1.0;
            int refinement_factor = 3;
            double refinement_factor_double = static_cast<double>(refinement_factor);
            std::uniform_real_distribution<double> dist(0.1, 1); //The minimum size of the sphere is set here

            double radius =  dist(eng);

            Eigen::Vector3d origin = generateRandomPoint(-1.0 + radius , 1.0 - radius, eng);

            // Create refined cells
            int refined_stencil_size = refinement_factor*stencil_size;
            std::vector<std::vector<std::vector<double>>> volumes(refined_stencil_size, std::vector<std::vector<double>>(refined_stencil_size, std::vector<double>(refined_stencil_size)));

            // Use refinement factor on parameters -> the grid is parted, so the parameters have to be scaled to be the same
            radius = radius * refinement_factor_double;
            origin = origin * refinement_factor_double;

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
                        
                        volumes[i][j][k] = volume_and_surface.getMoments().volume();
                                         //volume_and_surface.getMoments().centroid() / vol
                        if (visualize) {
                            surfaces.push_back(volume_and_surface.getSurface());
                        }
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
                                    sum += volumes[m][n][o];
                                }
                            }
                        }
                        vfrac [i][j][k] = sum / (cell_volume*(refinement_factor_double*refinement_factor_double*refinement_factor_double));
                    }
                }
            }

            // Generate vtk output
            if (visualize) {
                WriteField(stencil_size, coords, vfrac, "vfrac");
                WriteSurface(surfaces, "surface");
                std::cout << "VFRAC:" << std::endl;
                for (int i = 0; i < vfrac.size(); i++) {
                    for (int j = 0; j < vfrac[i].size(); j++) {
                        for (int k = 0; k < vfrac[i][j].size(); k++) {
                            std::cout << "[" << i << "," << j << "," << k << "] " << std::fixed << std::setprecision(3) << vfrac [i][j][k] << "   ";
                        }   
                    std::cout << std::endl;
                    }
                }
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


        std::vector<double> generate_State(int datapoint_type, int stencil_size, bool include_planes, double plane_bounds_coefficients){
            // Create a vector of volume fractions and surfaces
            // OLD std::vector<std::vector<std::vector<double>>> vfrac(stencil_size, std::vector<std::vector<double>>(stencil_size, std::vector<double>(stencil_size)));

            if (! include_planes) {
                plane_bounds_coefficients = 0.0;
            }

            switch (datapoint_type) {
                case 1:
                    return generate_Cylinder(stencil_size);
                case 2:
                    return generate_Sphere(stencil_size);
                case 3:
                    return generate_Sheet(stencil_size);
                default:
                    return generate_Paraboloid(stencil_size, 1, plane_bounds_coefficients);
            }
        }

        void generate_Data (std::vector<std::vector<double>>* statesV, std::vector<int>* labelsV, int no_datapoints, int stencil_size = 3, int datapoint_types_in = 4, bool include_planes = false, double plane_bounds_coefficients = 0.5){
            std::cout << no_datapoints << std::endl;
            // Initialize random number generator with a seed from current time
            std::srand(std::time(0));

            //int datapoint_types = datapoint_types_in;
            //    if (include_planes) {
            //        datapoint_types += 1;
            //    }

            for (int i=0; i<no_datapoints; i++) {
                // Generate the data, init with a random number 0 or 1, 0=plane, 1=paraboloid, 2=cylinder, 3=sphere, 4=sheet, decide if to include plane

                //int datapoint_type = std::rand() % datapoint_types;
                int datapoint_type = std::rand() % 4;
                //if (! include_planes) {
                //datapoint_type += 1;
                //}

                labelsV->push_back(datapoint_type);
                statesV->push_back(generate_State(datapoint_type, stencil_size, include_planes, plane_bounds_coefficients));
            }
        }
    };
}