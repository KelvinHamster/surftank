# Surftank

*This project was created under Falk Feddersen's research group at the Scripps Institution of Oceanography and supported by the Mark Walk Wolfinger Surfzone Processes Research Fund.*

This is a Matlab implementation of the boundary element method (BEM), used in order to advect the free surface of a fluid. 
Representing the fluid using potential flow, we use the BEM to solve Laplace's equation.
$$\nabla^2 \phi = \nabla \cdot (\nabla \phi) = 0$$
where the fluid velocity is $\nabla \phi$. This forces the fluid to be irrotational. Additionally, restrictions on the types of boundaries allowed in the BEM means that it must be describable by a simply connected domain -- that is, completely contiguous with no self intersections in the boundary. These restrictions do not inhibit the ability to model various shallow water phenomenon.
Additionally, since the Laplace's equation solver is core to the model, one can use it in isolation to solve Laplace boundary value problems.

While a Lagrangian mesh is used, the BEM is strictly Eulerian. Switching between the formulations requires the advection term $\big((\nabla \phi)\cdot \nabla\big)\phi = |\nabla\phi|^2$, which is calculated using the BEM. This model was built upon the work of shallow-water researchers like Stephan T. Grilli.

## Using the model

The following directories need to be included in `path`:
```
'tests'
'src'
'src/boundary'
'src/curve_functions'
'examples'
```
This list is also in `setup_testenv.m`, and running it from the repository root directory will add these to the path.
A few examples are provided in the `examples` directory. The main parts of the simulation are as follows

## Setting up the Simulation

The entirety of the simulation can be self-contained within the [`bem_sim` object](doc/bem_sim.md).
- For the constructor, a list of [boundaries](doc/boundary.md) is specified.
- For each boundary, [regridding](doc/boundary.md#property-regridding) parameters should be set.
- While the default parameters are most likely fine, corner (boundary intersection) resolution can be set in the [`boundary.property-intersection_ccw`](doc/boundary.md#property-intersection_ccw) struct for each boundary.
- [stepping parameters](doc/bem_sim.md#property-stepping) should be set.
- [Global events](doc/bem_sim.md#global-events), if necessary should be set.

### Running the simulation

A full step can be performed using [`s.full_step()`](doc/bem_sim.md#method-full_step). More information on how each step is handled can be found there.