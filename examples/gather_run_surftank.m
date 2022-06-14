function gather_run_surftank(Nx_FS,courant_target, tmax,ncfile,do_log,varargin)
%does a run using the surftank_IC initial condition.
%   Runs the surftank IC with x0 = 20, a0 = 1.5, theta = 64.5 degrees
%   with a courant lock on courant_target, up to tmax.
%   logging is done in ncfile, and if do_log is true, the print
%   statements are appended to a file of the format
%       surftanksimlog__YYYY_MM_DD.txt
%   unless the keyword argument log_filename is given.
%
% do_log defaults to true
%
% varargin - keyword arguments, as specified:
%
%   interp - M value for interpolation (default = 1)
%   interp_FS - M value for interpolation_FS (default = 1)
%   interp_sliding - M value for interpolation_sliding (default = 4)
%
%   integral_spline - integral spline flag in meta (default = 0)
%   load_from_nc - filename of netcdf file to load the initial condition
%       from. If loat_from_nc_t is not specified, t=0 is defaulted to
%   load_from_nc_t - the time to load from
%   sim_setup_handle - function handle that takes a bem_sim object and
%       outputs a bem_sim object with parameters that the user wants set.
%       The handle is used to modify the simulation right before it is run.
%       (default = @(x) x, the identity map)
%   log_filename - name of the logfile if do_log is true (0 for default name).
%   loop_callback - callback function relating to run_sim on_loop_callback
%   append_nc_callback - function relating to run_sim should_append_nc
%

interp = 1;
interp_FS = 1;
interp_sliding = 4;
integral_spline = 0;
load_from_nc = 'no';
load_from_nc_t = 0;
sim_setup_handle = @(x) x;
log_filename = 0;
append_nc_callback = @(t_last,t_now,steps) 1;
loop_callback = @(s) s;
do_fig = 0;
do_save = 0;

p = inputParser;
addOptional(p,'integral_spline',integral_spline,@(x) islogical(x) || (x==0) || (x==1))
addOptional(p,'load_from_nc',load_from_nc,@(x) ischar(x) || isstring(x))
addOptional(p,'load_from_nc_t',load_from_nc_t, @(x) (isnumeric(x) && isscalar(x)));
addOptional(p,'interp',interp, @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (mod(x,1) == 0) );
addOptional(p,'interp_FS',interp_FS, @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (mod(x,1) == 0) );
addOptional(p,'interp_sliding',interp_sliding, @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (mod(x,1) == 0) );
addOptional(p,'sim_setup_handle',sim_setup_handle, @(x) isa(x,'function_handle') );
addOptional(p,'log_filename',log_filename, @(x) x==0 || (ischar(x) || isstr(x)) );
addOptional(p,'append_nc_callback',append_nc_callback,@(x) isa(x,'function_handle'))
addOptional(p,'loop_callback',loop_callback,@(x) isa(x,'function_handle'))
addOptional(p,'do_plot_figs',do_fig,@(x) islogical(x) || (x==0) || (x==1))
addOptional(p,'do_save_figs',do_save,@(x) islogical(x) || (x==0) || (x==1))
parse(p, varargin{:});

interp = p.Results.interp;
interp_FS = p.Results.interp_FS;
interp_sliding = p.Results.interp_sliding;
integral_spline = p.Results.integral_spline;
load_from_nc = p.Results.load_from_nc;
load_from_nc_t = p.Results.load_from_nc_t;
sim_setup_handle = p.Results.sim_setup_handle;
log_filename = p.Results.log_filename;
should_append_nc = p.Results.append_nc_callback;
on_loop_callback = p.Results.loop_callback;
do_fig = p.Results.do_plot_figs;
do_save = p.Results.do_save_figs;


if ~exist('do_log','var')
    do_log = 1;
end

if strcmp(load_from_nc,'no')
    s = surftank_IC(Nx_FS,20,1.5,64.5*pi/180);
else
    man_temp = sim_netcdf_manager(load_from_nc,0);
    s = man_temp.get_sim(load_from_nc_t);
end

s.meta.integral_spline = integral_spline;
s.interpolation = get_interpolation_struct(interp,0,1);
s.interpolation_FS = get_interpolation_struct(interp_FS,0,interp_FS);
s.interpolation_sliding = get_interpolation_struct(interp_sliding,0,1);


t_target = tmax;
s.meta.sim_parallel_workers = 4;
do_nc_export = exist("ncfile","var") && (~strcmp(ncfile,"null"));
s.stepping.courant_lock = 1;
s.stepping.courant_number = courant_target;
s.stepping.timestep_method = 'taylor2';

s = sim_setup_handle(s);
if do_log
    if log_filename == 0
        logfile = sprintf('surftanksimlog__%s.txt',datestr(now,'yyyy_mm_dd'));
    else
        logfile = log_filename;             
    end
end
run_sim;
end

