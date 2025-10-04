
#include <iostream>
#include <vector>
#include <string>
//#include <torch/torch.h>

#include <vtkCellCenters.h>

#include "irl/ml_classification/ml_classifier.h"
#include "irl/ml_classification/vtk_in.h"



int main (int argc, char* argv[]) {

    int stencil_size = 3;

    // Net Parameters
    int input_size = stencil_size * stencil_size * stencil_size; // 27 if stencil_size=3 and only vof
    int hidden_size1 = 128;
    int hidden_size2 = 64;
    int hidden_size3 = 32;
    int output_size = 4;
    int batch_size = 64;

    //Training parameters
    double learning_rate = 0.001; //was 0.01 for SGD optimizer
    int no_batches = 512; // Should be divisible by batch size
    int epochs = 20;

    IRL::MLClassifier ml(stencil_size, input_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    ml.generateDataset(no_batches, batch_size);
    ml.trainModel(epochs, learning_rate, batch_size);

    // vtk reader
    std::string filename = "/home/quirin/mlcfd/Repositories/jet/nga.case";
    IRL::classify_simulation(ml, filename);

    return 0;
}