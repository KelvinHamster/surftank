# `bem_sim.m`

A handle class that stores all relevant information to a simulation.

### Properties

- [`boundaries`](#property-boundaries) - A cell array of [boundary objects](boundary.md).

- [`meta`](#property-meta) - A struct handling simulation-relevant information.

- [`stepping`](#property-stepping) - A struct handling time-stepping-relevant information.

- [`global_step_events`](#global-events) - A list of function handles representing the `step()` events not intrinsic to a certain boundary.

- [`global_step_event_orders`](#global-events) - A list of numbers representing the priorities/ordering of the `step()` events in `global_step_events`.

- [`global_update_events`](#global-events) - A list of function handles representing the `update_characteristics()` events not intrinsic to a certain boundary.

- [`global_update_event_orders`](#global-events) - A list of numbers representing the priorities/ordering of the `update_characteristics()` events in `global_update_events`.

### Methods

- [`constructor()`](#method-constructor)

- [`get_nodelist()`](#method-get_nodelist) - Builds the list of nodes across the whole boundary. In the [discretization of the method](bem_method.md#discretization), this is the set $X$.

- [`step()`](#method-step) - Calculates the step size `dt` that should be taken, then calls the [step events](#event-system).

- [`enforce_intersections()`](#method-enforce_intersections) - Ensures that corner nodes agree between boundaries.

- [`update_characteristics()`](#method-update_characteristics) - Calculates all necessary values around the boundary and calls the [update events](#event-system).

- [`regrid()`](#method-regrid)

- [`full_step()`](#method-full_step) - Performs all necessary actions to take one step, including `update_characteristics()`, `step()`, `enforce_intersections()`, and `regrid()`.

- [`append_update_event()`](#method-append_update_event) - Adds the given [global update event](#global-events) to this simulation.

- [`append_step_event()`](#method-append_step_event) - Adds the given [global step event](#global-events) to this simulation.

- [`interior_eval()`](#method-interior_eval) - Calculates $\phi$ or $\phi_t$ on the interior of the domain.

- [`is_point_inside()`](#method-is_point_inside) - Checks if a point is inside the domain.

## Property: boundaries

A cell array of [boundary objects](boundary.md), separating the domain into distinct sections. The sections must traverse the boundary in a counterclockwise motion, where the last node of `boundaries{i}` must match with the first node of `boundaries{i+1}`. Discontinuities in the first derivative (sharp corners) must occur exclusively on the intersection between two boundary sections.


## Property: meta

A struct to store information regarding the simulation. While the values below are universal, other values can be set. For example, sample code for the absorbing layer sets the variable `meta.L`, the domain length,
which is used within the event for the absorbing layer in order to find the start of the absorber.
<!---Reference example code for absorbing layer--->

- `parallel_workers` - specifies the number of worker threads to be used when calculating the boundary integrals. This must be a nonnegative integer -- any value passable into `parfor()`.
<!---Reference BEM Kn Kd; parfor() documentation--->

- `char_t` - the last simulation time [`update_characteristics()`](#method-update_characteristics) was called. This is used to ensure that the integrals are not calculated when they do not need to be. This takes the value of a float or the string/char 'null'.

- `g` - acceleration due to gravity. This is used when forcing the [free surface](boundary.md#free-surface-boundary).
<!---reference free surface forcing term--->


## Property: stepping

A struct to store information regarding how the simulation is time-stepped.

- `t` - the current time of the simulation.

- `dt` - the time step size. This is calculated and set in [`step()`](#method-step).

- `courant_number` - the Courant number, a nondimensional value $\frac{\Delta t}{\Delta x}u$, where $\Delta t$ is the time step size `dt`, $\Delta x$ is the characteristic mesh size (for us, this is node spacing), and $u$ is a velocity.
<!---Probably link to something like wikipedia? Potentially set function handle to calculate this number--->

- `courant_lock` - 0 if a fixed time step size `dt` should be taken. 1 if `dt` should be specified by `courant_target`.

- `courant_target` - the desired Courant number, only needs to be defined if `courant_lock == 1`.

- `dtmax` - the largest time step size possible. If `dt` goes above this value, it is clamped back down to `dtmax`.

- `dtmin` - If `dt` goes below this value, an error is thrown. If `dtmin` is not defined, then there is no minimum threshold.

## Global Events

Relevant to the properties `global_step_events`, `global_step_event_orders`,
`global_update_events`, `global_update_event_orders`. These store events (and their priority/order) that are not tied to a specific boundary. Step events are called within the [`step()`](#method-step) method, and update events are called within the [`update_characteristics()`](#method-update_characteristics) method. Elements of the `global_..._events` arrays are function handles that take only the `bem_sim` object. Elements of the `global_..._event_orders` arrays are floats specifying the ordering of the events. At the corresponding order, the function handle is called.

## Event System

There are two main event passes:

- Step events are called during [`step()`](#method-step) after the timestep size has been calculated. Free surface boundaries, for example, use a step event in order to update the free surface.
<!---Link to FS--->

- Update events are called during [`update_characteristics()`](#method-update_characteristics). These events are used to handle the calculation of all information relevant to the current time of the simulation, before, during, and after the BEM solving step. The first solve, for $\phi$ and $\frac{\partial\phi}{\partial n}$, occurs at ordering `0`, while the second solve, for $\frac{\partial\phi}{\partial t}$
and $\frac{\partial^2\phi}{\partial t\partial n}$, occurs at ordering `1`.
<!---Link to BEM solve--->

Each event has a number, which specifies where in the order that event should be called. Events are called in increasing order, for example, an order `0.5` event occurs before an event of order `1`.

## Method: Constructor
`sim = bem_sim(boundaries)`<br>
`sim = bem_sim(boundaries,varargin)`

Initializes the `bem_sim` object with `boundaries` as a cell array of boundary objects. The [`boundaries`](#property-boundaries) property will be set to this value.

Intersections are set up between adjacent boundaries, through the use of [`intersection_build_default()`](functions.md#intersection_build_default), which set up the [`intersection_CCW`](boundary.md#property-intersection_ccw) and [`intersection_CW_link`](boundary.md#property-intersection_cw_link) properties in [`boundary`](boundary.md).

Additional optional keyword arguments are permitted using parser. The only supported argument is

|Argument|Type|Default|Description|
|--------|----|-------|-----------|
|`parallel_workers`|nonnegative integer|0|The argument that gets passed into `parfor()` for the integral calculations in [`BIE_integ_calc()`](functions.md#bie_integ_calc). This is the number of parallel processes that get delegated.|

## Method: get_nodelist()

Generates a $N \times 2$ array of nodes, across all of the boundaries.
Ignores the last node of each boundary in order to prevent
repetitions on double-nodes.

This method returns the values `[nodes,starts]`, where `nodes` is the afforementioned array, and `starts` is an array that gives the index corresponding to the first node in each boundary. That is,
`nodes(starts(i),:) == boundaries{i}.boundary_nodes(1,:)`.


## Method: step()

Calculates the step size `dt` that should be taken, then calls boundary [step events](#event-system). These events are related to advecting the boundary, so corner nodes may not agree between boundary sections.
Hence, [`enforce_intersections()`](#method-enforce_intersections) is used to make sure they do agree.
<!---TODO stepsize decision--->


## Method: enforce_intersections()

Intersections (corners) are handled utilizing the [intersection structs](boundary.md#property-intersection_ccw). For each boundary, [`intersection_enforce()`](functions.md#intersection_enforce) is called.

## Method: update_characteristics()

Calculates all necessary values around the boundary. This is done by handling all [update events](#event-system). `update_characteristics()` does nothing if this method was already called this time step. To check this, the variable [`meta`](#property-meta)`.char_t` is used. If regridding occurs, `update_characteristics()` may need to be called more than once. Because of this, events should be designed so that being called multiple times does not invoke unexpected behavior.

At ordering 0, a BEM is used to populate `phi` and `phi_n` values within [`boundary.characteristics`](boundary.md#property-characteristics), where each boundary has one already filled, and the other is set after solving the matrix problem. To do this, [`BIE_integ_calc()`](functions.md#bie_integ_calc) is called to solve for [`Kd` and `Kn`](bem_method.md#discretization), the values being placed in [`boundary.characteristics`](#method-update_characteristics). Then, [`matrix_build()`](functions.md#matrix_build) is called to build the system of equations, which is then solved for the missing values among `phi` and `phi_n`.

At ordering 1, a second BEM is used to populate `phi_t` and `phi_tn`. This time, we do not need to call [`BIE_integ_calc()`](functions.md#bie_integ_calc).

Boundary, and global, events should be called in order to facilitate setting `phi`, `phi_n`, `phi_t`, and/or `phi_tn`, as well as other information potentially necessary for analysis.

## Method: regrid()

This method calls for a regrid on every boundary of the simulation, delegating to [`bdry_handle_regrid()`](functions.md#bdry_handle_regrid) for each one. Some regridding techniques require [`update_characteristics()`](#method-update_characteristics) to have been called. Setting [`boundary.regridding.precede_update`](boundary.md#property-regridding) to 1 tells `regrid()` that it is not necessary. Such boundaries are regridded first, in order specified by `boundary.regridding.order`. This ordering is in the same scheme as the [event system](#event-system). Afterwards, `update_characteristics()` is called, and boundaries with `precede_update` set to false are regridded. If any boundary is modified ([`bdry_handle_regrid()`](functions.md#bdry_handle_regrid) returns 1), then [`meta`](#property-meta)`.char_t` is set to `'null'`, flagging the model to require an update.

## Method: full_step()

This method performs everything necessary to step the simulation forward. If one wants to run a simulation, it suffices to properly set up the `bem_sim` object, and then calling `full_step()` in a loop.

Algorithmically, the typical loop is as follows:

1. Build and solve the boundary integral to calculate unknown boundary values (This is done in [`update_characteristics()`](#method-update_characteristics).)
2. Use these values to step the boundary to the next time step ([`step()`](#method-step)). Enforce compatibilities between boundary sections. (Ensuring the boundary is continuous, also in `step()`.)

`full_step()` does (1) if necessary, does (2), and then makes a [`regrid()`](#method-regrid) call, before calling [`update_characteristics()`](#method-update_characteristics). This way, `full_step()` can be called, then values can be obtained, such as energy, which require calculations from `update_characteristics()`.

## Method: append_update_event()
`bem_sim.append_update_event(handle,order)`

Adds the given [global update event](#global-events) to this simulation.
The event is at the given `order`, and calls `handle` during [`update_characteristics()`](#method-update_characteristics).

## Method: append_step_event()
`bem_sim.append_step_event(handle,order)`

Adds the given [global step event](#global-events) to this simulation.
The event is at the given `order`, and calls `handle` during [`step()`](#method-step).



## Method: interior_eval()
`phi = bem_sim.interior_eval(x,z)`

Calculates $\phi$ or $\phi_t$ on the interior of the domain, specified by the given point `(x,z)`. This calculates the [boundary integral](bem_method.md#the-boundary-integral) assuming $\alpha(x,z) = 1$, using [`bem_integrate()`](functions.md#bem_integrate). Instead of $\phi$, one can calculate $\phi_t$, by setting the optional argument `solve_phi_t` to true:

`phi_t = bem_sim.interior_eval(x,z,1)`



## Method: is_point_inside()
`inside = bem_sim.is_point_inside(x,z)`

Checks if a point is inside the domain. This uses the [even-odd rule](https://en.wikipedia.org/wiki/Even-odd_rule) over the segments around the boundary, using line segments (linearly interpolating the boundary). Because of this, the truthiness of the interior condition is not well defined close to the boundary, since a point may be on the interior using linear interpolation, but outside when using a spline.