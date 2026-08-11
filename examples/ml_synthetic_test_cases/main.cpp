#include "testcases.h"
#include "irl/ml_classification/ml_classifier.h"

int main (int argc, char* argv[]) {

    int stencil_size = 5;

    //Data parameters
    int no_batches = 4096 * 4 * 2;
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
    float epsilon_connectivity = 1e-12f;

    // Net Parameters
    int input_size = stencil_size * stencil_size * stencil_size 
    * (include_Moments >= 1 ? /*(include_Surface_Area ? 5 : 4) ---> remove the 4 again if uncommenting ->*/ 4 : 1)  // 4 if include_Moments >= 1 because we have vfrac + (mx,my,mz) per cell, otherwise just vfrac
    + (include_Moments >= 2 ? 6 : 0)  // +6 if include_Mome´nts >= 2 because we have (xx, yy, zz, xy, xz, yz) components of the 2nd moment tensor; otherwise none
    + (include_Eigenvalues ? 3 : 0);
    int hidden_size1 = 256;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 6; //CHANGED 4 to 6

    IRL::InertiaClassifier inertia_classifier(stencil_size, 1, 0.85, 1.5);
    //IRL::MLClassifier_E3NN ml(stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    IRL::MLClassifier ml(stencil_size, input_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    ml.updateDataParameters(
            no_batches,
            include_Moments,
            include_Surface_Area,
            include_Eigenvalues,
            paraboloid_coeff_stddev,
            hyperbolic_cylinder_opening_angle_stddev,
            sheet_coeff_stddev,
            sheet_thickness_stddev,
            cylinder_radius_stddev,
            radius_circle_min,
            radius_circle_max,
            sphere_radius_stddev,
            ellipsoid_subgrid_stddev,
            min_long_ellipsoid_axis,
            max_long_ellipsoid_axis,
            exact_2nd_moment,
            visualize,
            machineZero,
            lower_limit_subgrid,
            upper_limit_subgrid,
            class0_max_characteristic
        );                    

    shell_testcase(ml);
    return 0;
}