function [max_time,sim] = gather_run_bouncingsoliton(Nx, dt, a0, T, varargin)
%simulates and stores data from a linear standing wave simulation
%
%   Nx - number of points on the free surface
%   dt - time step
%   a0 - initial amplitude
%   T - how much time to simulate for (seconds)
%
%
%   workers - the number of worker threads in parallel tasks
%   do_gui - whether or not plots should be drawn every time step
%   save_plot - whether or not save plots to the figure directory. Saving
%               only occurs when do_gui is true
%   data_filename - file to save the simulation data to. 'null' for no saving.
%   rebalance_walls - after every time step, whether or not to shift the spacing of the walls
%   br_x - x coordinate of the bottom right double-node (if set,
%               rebalance_walls is forced true.) (default 'default')
%   br_z - z coordinate of the bottom right double-node (if set,
%               rebalance_walls is forced true.) (default 'default')
%   L - horizontal length of the simulation
%   sim_args - cell array for bem_sim constructor
%
% simulation data is saved into the runs directory in the format:
%   time, energy error, volume error, flux integral, mean water level, condition number\n
%
do_gui = 0;
save_plot = 0;
data_filename = 'null';
rebalance_walls = 0;
br_x = 'default';
br_z = 'default';
workers = 4;
plot_full = 0;
wall_resolution_factor = 1;
L = 40;
sim_args = {};

p = inputParser;
addOptional(p,'do_gui',do_gui,@(x) islogical(x) || (x==0) || (x==1));
addOptional(p,'save_plot',save_plot,@(x) islogical(x) || (x==0) || (x==1));
addOptional(p,'data_filename',data_filename,@(x) ischar(x) || isstring(x));
addOptional(p,'rebalance_walls',rebalance_walls,@(x) islogical(x) || (x==0) || (x==1));
addOptional(p,'workers',workers, @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (mod(x,1) == 0) );
addOptional(p,'plot_full',plot_full,@(x) islogical(x) || (x==0) || (x==1));
addOptional(p,'wall_resolution_factor',wall_resolution_factor, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'br_x',br_x, @(x) (isnumeric(x) && isscalar(x))|| strcmp(x,'default'));
addOptional(p,'br_z',br_z, @(x) (isnumeric(x) && isscalar(x))|| strcmp(x,'default'));
addOptional(p,'L',L, @(x) isnumeric(x) && isscalar(x) && (x >= 0) );
addOptional(p,'sim_args',sim_args,@(x) iscell(x))
parse(p, varargin{:});

do_gui = p.Results.do_gui;
save_plot = p.Results.save_plot;
data_filename = p.Results.data_filename;
rebalance_walls = p.Results.rebalance_walls;
br_x = p.Results.br_x;
br_z = p.Results.br_z;
workers = p.Results.workers;
plot_full = p.Results.plot_full;
wall_resolution_factor = p.Results.wall_resolution_factor;
L = p.Results.L;
sim_args = p.Results.sim_args;

% note here 'L' is called twice
s = bem_sim('type','soliton_reflect','L',40,'parallel_workers',workers,'Nx',Nx,'dt',dt,'a0',a0,'wall_resolution_factor',wall_resolution_factor,'L',L,sim_args{:});


s.stepping.courant_lock = 0;
s.meta.plot_ylim = [-0.01 2.2];

force_rebalance = 0;
if strcmp(br_x,'default')
    br_x = s.meta.L;
else
    force_rebalance = 1;
end
if strcmp(br_z,'default')
    br_z = min(s.boundary.z);
else
    force_rebalance = 1;
end

if force_rebalance
    rebalance_walls = 1;
%shift non-free surface to accomodate br_x/z
%br is the first double node because the bottom boundary is first
    l = find(s.boundary.doublenode,2,'first');
    r = find(s.boundary.doublenode,2,'last');
    z = s.boundary.z;
    x = s.boundary.x;
    s.boundary.x(l(2)+1:r(1)) = linspace(x(l(2)),br_x,r(1)-l(2));
    s.boundary.z(l(2)+1:r(1)) = linspace(z(l(2)),br_z,r(1)-l(2));
    s.boundary.x(r(1)+1:r(2)) = linspace(br_x,x(r(2)),r(2)-r(1));
    s.boundary.z(r(1)+1:r(2)) = linspace(br_z,z(r(2)),r(2)-r(1));
end

if exist('data_filename','var')
    fname = data_filename;
else
    fname = sprintf('../runs/Nx%d_dt%f_a0%f_bouncingsoliton_run.txt',Nx,dt,a0);
end

t_target = T;
do_nc_export = ~strcmp(data_filename,'null');
ncfile = data_filename;
on_loop_callback = @(s) on_loop(s,rebalance_walls);


do_fig = do_gui;
do_save = save_plot;

fignum = 0;%saving the figures

s = s.update_history();

base_energy = s.history.energy(1);
base_volume = s.history.volume(1);

x0_fit = 20;
a0_fit = 0.1;

run_sim;
    
max_time = s.stepping.t;
sim = s;

end


function s = on_loop(s,rebalance_walls)

        if rebalance_walls
            %set z values of left and right walls to be linearly spaced
            % last indices form the right wall
            z = s.boundary.z;
            x = s.boundary.x;
            l = find(s.boundary.doublenode,2,'first');
            r = find(s.boundary.doublenode,2,'last');
            s.boundary.z(r(1)+1:r(2)) = linspace(z(r(1)),z(r(2)),r(2)-r(1));
            s.boundary.z(l(1)+1:l(2)) = linspace(z(l(1)),z(l(2)),l(2)-l(1));
            s.boundary.x(r(1)+1:r(2)) = linspace(x(r(1)),x(r(2)),r(2)-r(1));
            s.boundary.x(l(1)+1:l(2)) = linspace(x(l(1)),x(l(2)),l(2)-l(1));
        end
end

