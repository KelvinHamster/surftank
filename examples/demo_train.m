
if ~exist("do_absorbing_piston","var")
do_absorbing_piston = 1;
end
if ~exist("do_absorbing_layer","var")
do_absorbing_layer = 1;
end

if ~exist("s","var")
    %create domain
    dt = 0.01;

    N_FS = 200;
    L = 40;
    N_wall = 20;
    
    a0 = 0.05; lam0 = 6;
    
    if do_absorbing_piston
        bds = {...
            free_surface_bdry('surface_nodes',[linspace(L,0,N_FS);zeros(1,N_FS)]'),...
            vertical_wavemaker_bdry('surface_nodes',[zeros(1,N_wall);linspace(0,-1,N_wall)]'),...
            solid_bdry('surface_nodes',[linspace(0,L,N_FS);zeros(1,N_FS)]'-1),...
            piston_absorber_bdry('surface_nodes',[L+zeros(1,N_wall);linspace(-1,0,N_wall)]')
        };
        bds{4}.h = 1;
    else
        bds = {...
            free_surface_bdry('surface_nodes',[linspace(L,0,N_FS);zeros(1,N_FS)]'),...
            vertical_wavemaker_bdry('surface_nodes',[zeros(1,N_wall);linspace(0,-1,N_wall)]'),...
            solid_bdry('surface_nodes',[linspace(0,L,N_FS);zeros(1,N_FS)]'-1),...
            solid_bdry('surface_nodes',[L+zeros(1,N_wall);linspace(-1,0,N_wall)]')
        };
    end
    
    s = bem_sim(bds);
    bds{2}.default_gen_functions(1,lam0,a0, 0,10);
    

    s.boundaries{2}.regridding.mode = 'linearize_uniform';
    s.boundaries{3}.regridding.mode = 'linearize_uniform';
    s.boundaries{4}.regridding.mode = 'linearize_uniform';
    
    s.meta.plot_xlim = [0-L*0.02,L + L*0.02];
    s.meta.plot_ylim = [-0.15,0.15];
    s.meta.parallel_workers = 10;
    s.stepping.dt = dt;
    s.stepping.courant_lock = 0;
    % s.stepping.courant_lock = 1;
    % s.stepping.courant_number = CFL;

    %=========absorbing beach
    if do_absorbing_layer
        FS_add_absorbing_layer(s,1, L-10,L,"nu0",0.5);
    end
    %=========
else
    fprintf("Variable 's' already exists. Assuming this is a bem_sim" + ...
        " that we should be running. Terminate and clear if not.\n")
end

stats_step = 10;
T = [];
VOLUME = [];
ENERGY = [];


TMAX = 50;
start_time = tic();
fignum = 0;
while s.stepping.t < TMAX
    %bdry_handle_regrid(s.boundaries{1});

    s.full_step();
    s.plot_full();
    vol = bem_integrate(s,@(x,z) x, 2);
    energy = calc_energy(s);
    
    if mod(fignum,stats_step) == 0
        T = [T s.stepping.t];
        VOLUME = [VOLUME vol];
        ENERGY = [ENERGY energy];
    end

    walltime = toc(start_time);
    wallsecs = mod(walltime,60);
    wallmins_ = (walltime - wallsecs)/60;
    wallmins = mod(wallmins_,60);
    wallhrs = (wallmins_ - wallmins)/60;
    logstr = sprintf('systime:%10.1fs (%02d:%02d:%06.3f): t=%10.4f (dt=%.6f, C=%.4f); E = %f\n',...
        walltime, wallhrs, wallmins, wallsecs, s.stepping.t, s.stepping.dt, s.stepping.courant_number, energy);
    fprintf(logstr);
    %saveas(1,sprintf('../figdump/fig%d.png',fignum)); fignum = fignum + 1;
    % fprintf("RG checks: global %f, local %f ", ...
    %     s.boundaries{1}.regridding.threshold_check_vals(1),...
    %     s.boundaries{1}.regridding.threshold_check_vals(2))
    % if s.boundaries{1}.regridding.did_regrid
    %     fprintf("(Regrid done)")
    % end
    % fprintf("\n")
end

