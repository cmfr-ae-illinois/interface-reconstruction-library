# ml_data_vtk example

This example is sued to test the data generation and visualize some of the shapes. A couple of important notes:

- Other than the sphere class, all other classes are made up of different shapes with different rules. To test them seperately, the swtich statement ´switch (shape_type)´ in ´generate_shape()´ in data_gen.h provides more options to generate a spcific shape. Before adding a shape to the dataset, it should be added to the end of this switch and visualized and tested.
- Not all classes have the ability to produce a vtk output to visualize the surfaces and volume fractions. For example all ellipsoid classes require the intersection of the grid aligned with the ellipsoid with the stencil grid. While coding a visual output for the stencil grid could be possible, it is not programmed at this time.
