
#include <iostream>
#include <vector>
#include <string>

#include <vtkCellCenters.h>

#include "irl/ml_classification/ml_classifier.h"
#include "irl/ml_classification/vtk_in.h"
#include "irl/ml_classification/inertia_classifier.h"
#include "irl/ml_classification/hybid_classifier.h"

#include "irl/ml_classification/data_gen.h"

#include "irl/ml_classification/ml_classifier_e3nn.h"




int main (int argc, char* argv[]) {
    
    int stencil_size = 5;

    //Data parameters
    int no_batches = 4096*4;
    int include_Moments = 1;
    double paraboloid_coeff_stddev = 0.1;
    double sheet_coeff_stddev = 0.1;
    double max_sheet_thickness = 0.5;
    double sheet_thickness_stddev = 0.0;
    double max_cylinder_radius = 0.5;
    double cylinder_radius_stddev = 0.0;
    bool include_truncated_cylinder = true;
    double max_sphere_radius = 0.5;
    double sphere_radius_stddev = 0.0;
    bool exact_2nd_moment = true;  // enable calculation of exact 2nd moments for data generation

    // Net Parameters
    int input_size = stencil_size * stencil_size * stencil_size 
    * (include_Moments >= 1 ? 4 : 1)  // 4 if include_Moments >= 1 because we have vfrac + (mx,my,mz) per cell, otherwise just vfrac
    + (include_Moments >= 2 ? 6 : 0)  // +6 if include_Moments >= 2 because we have (xx, yy, zz, xy, xz, yz) components of the 2nd moment tensor; otherwise none
    + (include_Moments >= 3 ? 3 : 0); // +3 if include_Moments >= 3 because we add the 3 eigenvalues of the inertia matrix; otherwise none
    int hidden_size1 = 256;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 4;

    //Training parameters
    double learning_rate = 0.001; //was 0.01 for SGD optimizer
    int batch_size = 64;
    int epochs = 5;

    //IRL::InertiaClassifier inertia_classifier(stencil_size, 1, 0.85, 1.5);
    //IRL::MLClassifier_E3NN ml(stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    IRL::MLClassifier ml(stencil_size, input_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    
    ml.updateDataParameters(no_batches, include_Moments,
                            paraboloid_coeff_stddev,
                            sheet_coeff_stddev,
                            max_sheet_thickness, sheet_thickness_stddev,
                            max_cylinder_radius, cylinder_radius_stddev, include_truncated_cylinder,
                            max_sphere_radius, sphere_radius_stddev,
                            exact_2nd_moment);                    
    
    //ml.generateDataset();
    ml.loadDataset("/home/quirin/mlcfd/Datasets/float/FirstMoment/s5_1M/data/data.bin");
    //ml.appendDataset("/home/quirin/mlcfd/Datasets/float/From1/s5_2M/data/data.bin", false);
    //ml.saveDataset("data");
    int canonicalize_symmetries = 48;
    bool preProcess = true;
    ml.canonicalize_data(canonicalize_symmetries);

    ml.updateTrainingParameters(learning_rate, batch_size, epochs);
    ml.trainModel();
    //ml.saveModel("model/ml_model.pt");
    //ml.loadModel("/home/quirin/mlcfd/Datasets/float/SecondFrom1/s5_1M/model/ml_model.pt");
    

    // vtk reader
    std::string filename = "/home/quirin/mlcfd/Repositories/jet/nga.case";
    IRL::classify_simulation(ml, filename, canonicalize_symmetries, preProcess, include_Moments);

    return 0;
}