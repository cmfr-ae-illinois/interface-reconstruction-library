# ml_stable_classification example

This example should be used once a dataset is decided. Since training is probablistic, different training seeds result in slightly different classifiers. The best way to decide between those is to look at the results on real simulation data. 

Given a dataset, a classifier is trained, and a number of simulations set in the code are evaluated every timestep. A second classifier is trained based on the same dataset and evaluated using the same simulations. The classifications of each classified cell of each timestep are recorded. It is therfore important that the simulations stay the same between these two runs, which is checked by comparing the total classifications made, if they are the same, it is likely that the correct simulations are compared. In practice, this should always be the case. This cylcle is repreated 10 times, but can be adjusted. Note that this runs rather long, up to multiple hours.

After all classifiers are trained and all simulations classified, the agreement per cell is calculated between each classifier and all other classifiers. The classifier with the highest average agreement to all others should be used in further analysis, like the test cases in the example ml_synthetic_testcases.

The input to this function is set in stable_classification_input.txt. Replace the dataset and simulation paths with the correct ones in your filesystem. The format of the textfile is:

`
dataset_path
number_of_simulations
simulation_1_name
simulation_1_data_directory
simulation_1_plic_directory
simulation_1_downsample_factor
simulation_2_name
simulation_2_data_directory
simulation_2_plic_directory
simulation_2_downsample_factor
`

More simulations can be added at the bottom of the textfile.