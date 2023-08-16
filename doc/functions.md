# Functions

## bem_integrate()
`result = bem_integrate(sim,value,mode)`<br>
`result = bem_integrate(sim,value,r_sing)`

integrates `value` over the boundary
`value` is one of the following:
- a function handle `value(x,z)` returning the integrand at a given
point `(x,z)`
- a constant number for a constant integrand (good for
measuring length by using constant 1)
- a cell array describing the integrand along each boundary
(`value{i}` is one of the following to specify it along boundary `i`)
    - array, where `value{i}(j)` is the integrand of node j, where the integral computes the interpolation of `value` (same shape as $\phi$).
    - function handle `value{i}(x,z)` giving the integrand at point
`(x,z)`
    - a constant number for a constant integrand. Use `value{i} = 0` to quickly exclude boundary `i` from the integral.

Function handles are only sampled at node points, so quickly
varying functions do not work well.
`mode` takes one of 4 values:

0. default - integrate value d Gamma
1. integrate value dx
2. integrate value dz
3. integrate value $G((x,z),(x_l,z_l)) ~d \Gamma$
4. integrate value $G_n((x,z),(x_l,z_l)) ~d \Gamma$

If `mode` 3 or 4 is specified, `r_sing` (1x2 array: (x,z) ) is used to
locate the singularity $(x_l,z_l)$.

If only the free surface should be integrated, value can be set to 0 on
non-free-surface nodes.



## bdry_handle_regrid()
`did_change = bdry_handle_regrid(bdry,force)`

Checks if regridding should be done on the given boundary `bdry`, doing it if it should. `force` skips the check and immediately regrids. The strategy used is based on the boundary [`regridding`](boundary.md#property-regridding) struct.

There are currently 3 modes for regridding. The mode for a given boundary is set in `bdry.regridding.mode`, and is expected to be one of:
- `"none"` - No regridding is done. When called, bdry_handle_regrid does nothing and returns `did_change = 0` for no change.
- `"linearize_uniform"` - If [`bdry.formulation`](boundary.md#property-formulation)`.type` is `"parameterized"`, then the boundary nodes are set according to a linspace over the parameterization. That is, node `i+1` is set to `bdry.formulation.surface_defn(i * (param_CCW - param_CW)/(n-1))`, where `param_CCW` and `param_CW` are specified in `bdry.formulation`, and `n` is the number of nodes. <br>
Otherwise, for other formulation types, the boundary nodes are set to linspace in x and z. That is,
node `i+1` is set to `lastnode * i/(n-1) + firstnode * (n-1-i)/n`.
- `"nodeshift"` - Shifts nodes along the interpolation of [`bdry.boundary_nodes`](boundary.md#property-boundary_nodes) such that they are evenly spaced according to [weighed length](regridding.md#regridding-heuristics). A detailed explanation is provided below.

### handle_regrid `"nodeshift"`

This mode is useful for the free surface. Nodes are shifted along the interpolation of [`bdry.boundary_nodes`](boundary.md#property-boundary_nodes) such that they are evenly spaced according to [weighed length](regridding.md#regridding-heuristics). For [`"shift_by_curve2d_normu"`](regridding.md#length-weight-curve2d_normu), a spline is used.

For [`"shift_by_curve3d_normu"` and `"shift_by_curve3d_maxphis"`](regridding.md#length-weight-curve3d), the type is specified by [`bdry.regridding`](boundary.md#property-regridding)`.interp_type`. This can be
- `'linear'` - does linear interpolation (order 1; 2 nodes)
- `'sliding3'` - 3rd order (4 node) lagrange interpolation over a [sliding](boundary.md#interpolation-sliding-segments) element.
- `'sliding5'` - 5rd order (6 node) lagrange interpolation over a [sliding](boundary.md#interpolation-sliding-segments) element.
- `'spline'` - does spline interpolation

The value of $a$ is taken from [`bdry.regridding`](boundary.md#property-regridding)`.a`. Issues have arisen based on regridding near piston-type boundaries, so there are parameters for ignoring nodes near the corner:
- `bdry.regridding.startnode` - a number specifying the index of the first regridded node. If not set, `handle_regrid()` defaults to a value of 1 (the corner node).
- `bdry.regridding.endnode` - a number specifying the index of the last regridded node. If not set, `handle_regrid()` defaults to a value of [`bdry.node_count`](boundary.md#property-node_count) (the corner node), unless `endpadding` is specified, as below:
- `bdry.regridding.endpadding` - a number specifying how far from the end does regridding start. This is only used if `endnode` (above) is not set. Effectively, this option sets `endnode = node_count - endpadding + 1`, as if this parameter was `startnode` on the other side (Setting `startnode` and `endpadding` to the same value ignores the same number of nodes on either side of the boundary).

Since lengths over ignored nodes aren't taken into account the weighted length is modified to decrease the difference in lengths towards either side of the boundary. The optional parameter `bdry.regridding.buffer_length` (defaults to 0 if not specified) gives the number of nodes used to transition from the unregridded spacing to the regridded spacing. The details of this can be found [here](regridding.md#transition-region). An additional optional parameter `bdry.regridding.reweight_min` is used to bound the weights below (default: `0.0`) in the transition region.


## BIE_integ_calc()
`BIE_integ_calc(sim)`

Calculates $K_d$ and $K_n$, the [boundary integral over shape functions](bem_method.md#discretization). These are placed into [`boundary.characteristics`](boundary.md#property-characteristics) as fields `Kd` and `Kn`. The specifics of these arrays are detailed in the [`characteristics` documentation](boundary.md#property-characteristics).

## get_interpolation_struct()


Returns an interpolation struct given the order and bounds. Nodes are spaced evenly, so, given $x_0$ and $x_M$ as the first and last nodes, the interior nodes are $x_i = x_0 + \frac{i}{M}(x_M - x_0)$, so that
$$L_j(x_i) = \delta_{i,j}$$
where $\delta_{i,j}$ is the Kronecker delta (1 if $i=j$ and 0 otherwise). Information on how the interpolation is used is provided [here](boundary.md#interpolation).

### Interpolation struct format

The struct stores coefficients of the Lagrange polynomials, as below.

- `M` - the order of the interpolation. This is the polynomial degree, which is used to interpolate `M+1` points.
- `Lagrange` - the coefficients of the Lagrange polynomials. The first array index references the polynomial, while the second index references the term of the polynomial. Thus, this array is `(M+1) x (M+1)`.
$$
L_j(x)  = \sum_{k=0}^M \verb+Lagrange+(j+1,k+1)~x^k
$$
- `Lagrange_prime` - the coefficients of the first derivative of the Lagrange polynomials. This `(M+1) x M` array is structured similar to `Lagrange`:
$$
    L_j'(x) = \sum_{k=0}^{M-1} \verb+Lagrange_prime+(j+1,k+1)~x^k
$$
- `Lagrange_primeprime` - the coefficients of the second derivative of the Lagrange polynomials. This `(M+1) x (M-1)` array is structured similar to `Lagrange`:
$$
     L_j''(x)= \sum_{k=0}^{M-2} \verb+Lagrange_primeprime+(j+1,k+1)~x^k 
$$

### get_interpolation_struct() arguments

| Argument | Type | Description |
|----------|------|-------------|
|`M`|positive integer| The order of the interpolation. This is the same as `M` in the struct.|
|`x0`| float, less than `xm`| The position of the first node.|
|`xm`| float, greater than `x0`| The position of the last node.|

## intersection_build_default()
`intersection_build_default(CW,CCW)`

Initializes the intersection struct between the two [`boundary`](boundary.md) objects. `CW` is the boundary on the clockwise side, and `CCW` is the boundary on the counterclockwise side.

In the [`boundary`](boundary.md) objects,
`CW.`[`intersection_CCW`](boundary.md#property-intersection_ccw) is built as below,
and `CCW.`[`intersect_CW_link`](boundary.md#property-intersect_cw_link)
is set to be.

The `intersection_CCW` struct has two fields:
- `CCW_link` - A reference to the boundary on the counterclockwise side (this is set to `CCW`).
- `resolution` - The strategy used to resolve the intersection in [`intersection_enforce()`](#intersection_enforce). The details of each are provided in the [`boundary.intersection_CCW`](boundary.md#property-intersection_ccw) struct.

For a corner involving the intersection between a [free surface](boundary.md#free-surface-boundary) and [`"nodes"` formulation type](boundary.md#property-formulation), `resolution` is set to match to the free surface (`"match_CW"` if `CW` is the free surface, or `"match_CCW"` if `CCW` is the free surface). If the corner is the intersection between a solid boundary and a piston-type boundary ([vertical wavemaker](boundary.md#vertical-wavemaker-boundary) or [piston absorber](boundary.md#absorbing-piston-boundary)), then `resolution` is set to match the piston (`"match_CW"` if `CW` is the piston, or `"match_CCW"` if `CCW` is the piston). Otherwise, `resolution` is set to `"ignore"`.

This function is called in the [`bem_sim constructor`](bem_sim.md#method-constructor).

## intersection_enforce()
`intersection_enforce(CW)`

Since corner nodes are written in the [`boundary`](boundary.md) object on both sides, their positions may not match. For example, the free surface is advected, so after [`step()`](bem_sim.md#method-step), an adjacent solid boundary will have a mismatched corner, giving a discontinuity in the boundary. `intersection_enforce()` forces agreement using a strategy specified in the [`intersection_CCW`](boundary.md#property-intersection_ccw) struct. The argument `CW` specifies the boundary on the clockwise side of the corner in question. The strategy used is given from `intersection_CCW.resolution`, detailed in the [`boundary.intersection_CCW`](boundary.md#property-intersection_ccw) documentation.


## matrix_build()
`[A,b] = matrix_build(sim)`<br>
`[A,b] = matrix_build(sim,'b_fieldsrc_u','phi_t','b_fieldsrc_q','phi_tn')`<br>
`b = matrix_build(sim,'build_A',0)`<br>
`A = matrix_build(sim,'build_b',0)`<br>

Generates the matrix / vector for the [BEM linear system](bem_method.md#discretization). The relevant [`bem_sim`](bem_sim.md) object is given in `sim`, while additional optional arguments can be provided. The system that should be solved is `A*x = b`, for `x`.

|Argument|Type|Default|Description|
|--------|----|-------|-----------|
|`build_A`|logical|1|Whether or not the matrix `A` should be generated. Setting this to zero is good for when the boundary values change, but the overall boundary conditions do not.|
|`build_b`|logical|1|Whether or not the vector `b` should be generated.|
|`b_fieldsrc_u`|string/char array| `'phi'`| The field in [`boundary.characteristics`](boundary.md#property-characteristics) that should be sampled for the harmonic function on the boundary. This is usually $\phi$ or $\phi_t$.|
|`b_fieldsrc_q`|string/char array| `'phi_n'`| The field in [`boundary.characteristics`](boundary.md#property-characteristics) that should be sampled for the normal derivative of the harmonic function on the boundary. This is usually $\phi_n = \frac{\partial \phi}{\partial n}$ or $\phi_{tn}$.|

The first return value is the matrix `A` if `build_A == 1`. Otherwise, it is `b` (if `build_b == 1`). The second return value is given as `b` only when both `build_A` and `build_B` are true.