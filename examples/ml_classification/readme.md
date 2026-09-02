# ml_classification example

This example demonstrates the basic functionality of the classifier.
There are two classifiers in working order right now:
1. The inerta classifier: From an information of either purely VOF data (bad since that means the first moments are approximated) or VOF + barycenter data, this classifier constructs the eigenvalues of the inertia tensor and uses heuristic rules to classify them.
2. The ml classifier: From the same two possible sets of information, this classifier uses ML for classification.

The focus of this example is the latter. 
First, a whole bunch of parameters are defined. The definitions in this example are consistent with the standard values, but are assigned nevertheless to give a user an easy way to change them and experiment with them. The most important parameters for any classifier, regardless if inertia or ML-based, are:

```
    // Classifier parameters
    int stencil_size = 5; // sets the size of the neighbourhood of the classified cell
    int canonicalize_symmetries = 48; // Defines what rules of cannonicalization are used. Irellevant for inertia classifier
    float noise_stddev = 0.0f; // If larger than 0.0f, noise is added during the preprocessing of a stencil. Most useful for the ML classifier before training
    int include_Moments = 1; // If 0, only uses VOF information of each cell, if 1, uses the x,y,z value of the first moment of each cell as well, this code has partial functionality with second moments per cell as well
    bool include_Surface_Area = false; // Option to include the surface area of each cell
    bool include_Eigenvalues = false; // Option to include the three inertia tensor eigenvalues of the whole stencil, calculated from the existing information, at the end of the vector
    float epsilon_connectivity = 1e-12f; // If a VOF fraction is below this, it is counted as 0, leading to potential disconnections from the central cell
```

While the inertia classifier doesn't use some of this information, these parameters are pulled from the classifier when assembling the information that ultimately gets passed to it, in one big vector, based on which it performs the classification (see irl/ml_classification/vtk_in.h for example, where this vector is constructed from simulation data). The ml classifier further requires the net architecture as input. Data parameters are used purely for data generation, and training parameters for training. These two have to be updated before using the classifier when any value is different from the standard parameters.

In this example, a dataset is generated with includeMoments = 1, and it is then converted to a dataset containing only VOF information (zeroth moments for unit cells) using retain_only_0th_moments();. When comparing datasets, it is instead best practice to generate the data with includeMoments = 1, save it somewhere it is not overwritten, then load it, use retain_only_0th_moments();, and save it to another place or under another name. This results in two datasets containing the same generated shapes, but using different amounts of information. If the eigenvalues are desired as well, they will be calculated from the existing information during the preprocessing step. This step also applies the local CCL, canonicalization, and has the possibility to add noise.

The net can now be trained and the resulting model saved for analysis, and optinally transferred into a torchless classifier by exporting the weights and biases header file.

For the classification of a simulation to make sense, it needs to have some way of representing subgridscale features. This can be done using a downsampling facor of 2 or greater, or subgrid reconstruction models, for which the downsampling factor = 0.

Now that basic functionality of the classification is understood, the minimal required dataset size for a custom dataset should be determined. Go to the example ml_dataset_size.