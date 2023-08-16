# bem_sim_file_manager.m

This class handles saving simulation data to a netCDF file. A buffer is used to store data before writing, so make sure that the buffer is flushed when complete.

The file manager object operates on a directory, which contains different stages. Each stage has an attached `.mat` file storing the [`bem_sim`](bem_sim.md) object, and a netCDF `.nc` data file storing boundary information.

When initializing a stage, the `.mat` file is saved to store a template for what the stage's object should look like. Then, for every snapshot to save, boundary information, namely $x,z,\phi,\phi_n,\phi_t,\phi_{tn}$ and the simulation time, is stored into the `.nc` file.

Then when loading, the [`bem_sim`](bem_sim.md) object is loaded from the `.mat` file, and every boundary has its corresponding values set from the `.nc` file.

### Methods

- [Constructor](#method-constructor)
- [`flush_buffer()`](#method-flush_buffer) - Pushes all snapshots in a a given stage from the buffer into the data file.
- [`get_sim()`](#method-get_sim) - Loads a [`bem_sim`](bem_sim.md) object from a file at the given time.
- [`get_stages()`](#method-get_stages)
- [`save_time()`](#method-save_time)
## Stages

Loading data references the stage's `.mat` file to get the [`bem_sim`](bem_sim.md) object, so different stages are used in order to deal with cases where the object changes. For example, if the number of nodes on a boundary changes, or if a boundary is replaced with something else, then it may be useful to use a different stage.





## Method: Constructor
`man = bem_sim_file_manager(dir)` <br>
`man = bem_sim_file_manager(dir,buffer_size)`

Builds a file manager, operating in the provided directory. This is the directory that will contain the `.mat` and `.nc` files. If the directory does not exist, it is made. The argument `buffer_size` can be provided, which sets the number of snapshots any given stage/boundary can store before needing to [flush](#method-flush_buffer) to a file. By default, this is `50`.

## Method: flush_buffer()
`manager.flush_buffer(stage,bd_mask)`

Pushes all snapshots in a a given stage from the buffer into the data file.
Within the given `stage`, all of the boundaries within `bd_mask` have their own buffer with maximum size `buffer_size`. When calling `flush_buffer()`,
this data is pushed into the corresponding stage `.nc` file. If the `.nc` or `.mat` files do not exist, then they are created. Moreover, if the boundary groups within the `.nc` do not exist, they will be created.

`bd_mask` must be given as a $1\times n$ array of indices, in order for the for-loop syntax `for bd_index = bd_mask` to work.

## Method: get_sim()
`s = manager.get_sim(stage,t)`

Loads a [`bem_sim`](bem_sim.md) object from a file at the given time. Within the given `stage`, the `.mat` file loads the [`bem_sim`](bem_sim.md) object that serves as a template for the sim. For each boundary, the closest snapshot in the `.nc` file to the given time `t` is used to set the data on that boundary. If different boundaries are stored with different time intervals, then it is advised to call [`bem_sim.enforce_intersections()`](bem_sim.md#method-enforce_intersections).

## Method: get_stages()
`stages = manager.get_stages()`

Returns a list (cell array) of stage names this manager has access to. Entries are what can be passed in as the `stage` argument for the other methods in this object.


## Method: save_time()
`manager.save_time(sim,stage,bd_mask)`<br>
`manager.save_time(sim,stage,bd_mask, force_flush)`

Saves a current snapshot of the `sim` [`bem_sim`](bem_sim.md) object into a file (but buffered). If a buffer is full, then all boundaries in `bd_mask` are flushed into disk storage.

The simulation is saved into the stage given by `stage`. Only the boundaries specified by `bd_mask` are saved. Since each boundary has its own time series, this is not a problem. One can save the free surface every time step, while other boundaries may be saved at a lower rate (or not at all).

`bd_mask` must be given as a $1\times n$ array of indices, in order for the for-loop syntax `for bd_index = bd_mask` to work.

The optional `force_flush` argument can be given to force a flush command if true, even when no buffer is full.