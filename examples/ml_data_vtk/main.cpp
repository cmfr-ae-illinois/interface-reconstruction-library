#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <random>
#include <Eigen/Dense>
#include <vector>
#include <vtkCellCenters.h>

#include "irl/ml_classification/data_gen.h"
#include "irl/ml_classification/inertia_classifier.h"
#include "irl/ml_classification/vtk_in.h"

int main(int argc, char* argv[]) {
    /*
    int stencil_size = 5;

    IRL::InertiaClassifier inCl(stencil_size, 0, 0.85, 1.5);


    IRL::Data_gen data_gen;

    // Map integers to class names
    std::map<int, std::string> classNames = {
        {0, "paraboloid"},
        {1, "ligament"},
        {2, "sphere"},
        {3, "sheet"}
    };

    // 100 random examples per class without visualize
    int numSamples = 100;
    for (int trueClass = 0; trueClass < 4; ++trueClass) {
        // Counters
        std::map<int, int> counts = { {0,0}, {1,0}, {2,0}, {3,0} };

        for (int s = 0; s < numSamples; ++s) {
            std::vector<double> flattened_state = data_gen.generate_State(trueClass, stencil_size, true, false);
            int detectedClass = inCl.classify(flattened_state);
            counts[detectedClass]++;
        }

        std::cout << "True Class: " << classNames[trueClass] << ", Detected classes: "
                  << counts[0] << " paraboloids, "
                  << counts[1] << " ligaments, "
                  << counts[2] << " spheres, "
                  << counts[3] << " sheets.\n";
    }
    */
    int stencil_size = 3;

    IRL::InertiaClassifier inCl(stencil_size, 0, 0.85, 1.5);

    // vtk reader
    std::string filename = "/home/quirin/mlcfd/Repositories/jet/nga.case";
    IRL::classify_simulation(inCl, filename);

    return 0;
}