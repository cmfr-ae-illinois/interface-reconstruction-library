# ml_synthetic_test_cases example

This example can be used to evaluate the performance of classifiers based not on simulation data, which can be difficult to judge, but on synthetically generated test cases. These use paraboloids to represent a certain shape, similarly to how the data generation works. But instead of generating a shape in only a stencil, a larger domain is used, of which stencils are extracted and classified. Since we know what to expect, a visualization is unneccessary. Instead, this provides the possibility of numerical and interpretable outputs. Note that many of the common functions from irl/ml_classification/common_functions.h are used to simplify this.

## Spherical shell testcase

This test case creates two nested spheres with a varying thickness between them. The classifier should correctly identify the whole domain as subgrid when below the set decision characterstic length, and classify the whole domain as well-resolved when above. Tweaking this example can yield a good resolution at this critical point. Furthermore, the data does not contain a pure spherical sheet, so this example further provides an out of distribution evaluation of the classifier.