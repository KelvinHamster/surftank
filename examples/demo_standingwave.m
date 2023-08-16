%Flat bottom, start with standing wave
%=================================



%create domain
if ~exist("L","var")
    L = 10;
end
if ~exist("N_wall","var")
N_wall = 9;
end
if ~exist("N_floor","var")
N_floor = 101;
end
if ~exist("N_FS","var")
N_FS = 100;
end
if ~exist("MODE","var")
MODE = 2;
end
if ~exist("a0","var")
a0 = 0.2;
end
if ~exist("parallel_workers","var")
    parallel_workers = 10;
end
if ~exist("dt","var")
    dt = 0.01;
end
if ~exist("TMAX","var")
    TMAX = 35;
end
if ~exist("do_absorbing_piston","var")
do_absorbing_piston = 1;
end
if ~exist("do_absorbing_layer","var")
do_absorbing_layer = 1;
end


if do_absorbing_piston
    s = bem_sim({ ...
        solid_bdry('surface_nodes',[linspace(0,0,N_wall);linspace(a0,-1,N_wall)]')...vertical left wall
        solid_bdry('surface_nodes',[linspace(0,L,N_floor);linspace(-1,-1,N_floor)]')...bottom floor
        piston_absorber_bdry('surface_nodes',[linspace(L,L,N_wall);linspace(-1,a0*(-1)^MODE,N_wall)]')...vertical right wall
        free_surface_bdry('surface_nodes',[linspace(L,0,N_FS);a0*cos(linspace(L,0,N_FS).*(MODE*pi/L))]')...free surface
        });
    s.boundaries{3}.h = 1;
else
    s = bem_sim({ ...
        solid_bdry('surface_nodes',[linspace(0,0,N_wall);linspace(a0,-1,N_wall)]')...vertical left wall
        solid_bdry('surface_nodes',[linspace(0,L,N_floor);linspace(-1,-1,N_floor)]')...bottom floor
        solid_bdry('surface_nodes',[linspace(L,L,N_wall);linspace(-1,a0*(-1)^MODE,N_wall)]')...vertical right wall
        free_surface_bdry('surface_nodes',[linspace(L,0,N_FS);a0*cos(linspace(L,0,N_FS).*(MODE*pi/L))]')...free surface
        });
end

s.boundaries{1}.regridding.mode = 'linearize_uniform';
s.boundaries{end-1}.regridding.mode = 'linearize_uniform';


s.meta.plot_xlim = [0-L*0.02,L + L*0.02];
s.meta.plot_ylim = [-1.5 1.5]*a0;
s.meta.parallel_workers = parallel_workers;
s.stepping.courant_lock = 0;
s.stepping.dt = dt;
s.plot_full();

if do_absorbing_layer
    FS_add_absorbing_layer(s,4,0,L,'nu0',0.5,"h",1);
end


start_time = tic();
while s.stepping.t < TMAX
    s.full_step();
    s.plot_full();

    ener = calc_energy(s);

    walltime = toc(start_time);
    wallsecs = mod(walltime,60);
    wallmins_ = (walltime - wallsecs)/60;
    wallmins = mod(wallmins_,60);
    wallhrs = (wallmins_ - wallmins)/60;
    logstr = sprintf('systime:%10.1fs (%02d:%02d:%06.3f): t=%10.4f (dt=%.6f, C=%.4f); E = %f\n',...
        walltime, wallhrs, wallmins, wallsecs, s.stepping.t, s.stepping.dt, s.stepping.courant_number, ...
        ener);
    fprintf(logstr);
    title(sprintf('t = %.3f, E = %.3e', s.stepping.t, ener));
end
