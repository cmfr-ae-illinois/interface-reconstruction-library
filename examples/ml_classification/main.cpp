
#include <iostream>
#include <vector>
#include <string>

#include <vtkCellCenters.h>

#include "irl/ml_classification/ml_classifier.h"
#include "irl/ml_classification/vtk_in.h"
#include "irl/ml_classification/inertia_classifier.h"



int main (int argc, char* argv[]) {
    
    int stencil_size = 5;

    // Net Parameters
    int input_size = stencil_size * stencil_size * stencil_size; // 27 if stencil_size=3 and only vof
    int hidden_size1 = 128;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 4;
    int batch_size = 64;

    //Training parameters
    double learning_rate = 0.001; //was 0.01 for SGD optimizer
    int no_batches = 256;
    int epochs = 20;

    //Data parameters
    int include_Moments = 0;
    double paraboloid_coeff_stddev = 0.1;
    double sheet_coeff_stddev = 0.1;
    double max_sheet_thickness = 0.5;
    double sheet_thickness_stddev = 0.0;
    double max_cylinder_radius = 0.5;
    double cylinder_radius_stddev = 0.0;
    double max_sphere_radius = 0.5;
    double sphere_radius_stddev = 0.0;

    IRL::MLClassifier ml(stencil_size, input_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    ml.generateDataset(no_batches, batch_size, include_Moments,
                        paraboloid_coeff_stddev,
                        sheet_coeff_stddev, max_sheet_thickness, sheet_thickness_stddev,
                        max_cylinder_radius, cylinder_radius_stddev, 
                        max_sphere_radius, sphere_radius_stddev);
    ml.trainModel(epochs, learning_rate, batch_size);

    // vtk reader
    std::string filename = "/home/quirin/mlcfd/Repositories/jet/nga.case";
    IRL::classify_simulation(ml, filename);

    // Creata a text file to store parameters
    std::ofstream file("Parameters.txt");
    if (!file.is_open()) {
        std::cerr << "Error: Could not create Parameters.txt\n";
        return 1;
    }

    // Write parameters
    file << "=== Network Parameters ===\n";
    file << "stencil_size = " << stencil_size << "\n";
    file << "input_size = " << input_size << "\n";
    file << "hidden_size1 = " << hidden_size1 << "\n";
    file << "hidden_size2 = " << hidden_size2 << "\n";
    file << "hidden_size3 = " << hidden_size3 << "\n";
    file << "output_size = " << output_size << "\n";
    file << "batch_size = " << batch_size << "\n\n";

    file << "=== Training Parameters ===\n";
    file << "learning_rate = " << learning_rate << "\n";
    file << "no_batches = " << no_batches << "\n";
    file << "epochs = " << epochs << "\n\n";

    file << "=== Data Parameters ===\n";
    file << "include_Moments = " << include_Moments << "\n";
    file << "paraboloid_coeff_stddev = " << paraboloid_coeff_stddev << "\n";
    file << "sheet_coeff_stddev = " << sheet_coeff_stddev << "\n";
    file << "max_sheet_thickness = " << max_sheet_thickness << "\n";
    file << "sheet_thickness_stddev = " << sheet_thickness_stddev << "\n";
    file << "max_cylinder_radius = " << max_cylinder_radius << "\n";
    file << "cylinder_radius_stddev = " << cylinder_radius_stddev << "\n";
    file << "max_sphere_radius = " << max_sphere_radius << "\n";
    file << "sphere_radius_stddev = " << sphere_radius_stddev << "\n";

    file.close();
    
    /*
    int stencil_size = 3;

    IRL::InertiaClassifier inCl(stencil_size, 0, 0.85, 1.5);

    // vtk reader
    std::string filename = "/home/quirin/mlcfd/Repositories/jet/nga.case";
    IRL::classify_simulation(inCl, filename);
    */

    return 0;
}