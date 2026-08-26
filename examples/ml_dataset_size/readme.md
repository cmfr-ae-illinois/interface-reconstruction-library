# ml_dataset_size example

This example is used to show the minimum dataset size for good training results based on specific dataset parameters.

The dataset parameters for this example are the standard values and therefore don't show up in the code, but can be adjust in the same way as shown in the ml_classification example. Data is generated, the classifier trained, the loss evaluated, new data appended to the previous data and so on. First, a small step size is used, which is later coarsened as the projected improvement of the loss decreases. Once satisfied with a dataset size, multiple classifiers can be trained and the most agreeing one found by using the ml_stable_classification example.