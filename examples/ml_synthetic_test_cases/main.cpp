#include "testcases.h"
#include "irl/ml_classification/ml_classifier.h"

int main (int argc, char* argv[]) {
    IRL::InertiaClassifier inertia_classifier(stencil_size, 1, 0.85, 1.5);
    //IRL::MLClassifier_E3NN ml(stencil_size, hidden_size1, hidden_size2, hidden_size3, output_size);
    IRL::MLClassifier ml;

    ml.loadmodel("ml_classifier_model.pt");

    shell_testcase(ml);
    return 0;
}