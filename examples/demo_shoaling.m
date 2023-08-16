% Waves running up a beach
%=================================

if ~exist("do_FS_regrid","var")
    do_FS_regrid = 1;
end
if ~exist("TMAX","var")
    TMAX = 35;
end
if ~exist("do_plot","var")
    do_plot = 1;
end
if ~exist("do_plotsave","var")
    do_plotsave = 0;
end
if ~exist("log_regrids","var")
    log_regrids = 1;
end
if ~exist("parallel_workers","var")
    parallel_workers = 10;
end

if ~exist("s","var")
    %create domain
    if ~exist("N_FS","var")
        N_FS = 400;
    end
    
    if ~exist("a0","var")
        a0 = 1.5;
    end

    if ~exist("CFL","var")
        CFL = 0.1;
    end

    if ~exist("lam0","var")
        lam0 = 20;
    end

    %boundary wave generation/absorption
    if ~exist("use_wavemaker","var")
        use_wavemaker = 1;
    end
    if ~exist("use_piston_absorber","var")
        use_piston_absorber = 1;
    end
    if ~exist("use_absorbing_layer","var")
        use_absorbing_layer = 1;
    end
    
    %give soliton IC? or still water?
    if ~exist("starting_soliton","var")
        starting_soliton = 0;
    end

    if ~exist("do_fileman","var")
        do_fileman = 0;
    end
    if ~exist("dt_store","var")
        dt_store = 0.1;
    end
    if ~exist("dtmax","var")
        dtmax = 0.02;
    end


    
    % surftank_IC(Nx,x0,a0,theta,varargin)
    bds = surftank_IC(N_FS,20,a0 * starting_soliton,64.5*pi/180);

    %===replace boundaries for wavemaker
    if use_piston_absorber
        bd_n = bds{end}.boundary_nodes;
        bds{end} = piston_absorber_bdry('surface_nodes',bd_n);
    end
    if use_wavemaker
        bd_n = bds{2}.boundary_nodes;
        bds{2} = vertical_wavemaker_bdry('surface_nodes',bd_n);
    end
    %===
    
    
    s = bem_sim(bds);

    %regrid right/left walls, regardless of absorber/wavemaker
    bds{end}.regridding.mode = 'linearize_uniform';
    bds{2}.regridding.mode = 'linearize_uniform';
    if use_wavemaker
        bds{2}.default_gen_functions(1,lam0,a0, 0,10);
        %wavemaker moves left/right, so we may need this
        bds{3}.regridding.mode = 'linearize_uniform';
    end

    if use_piston_absorber
        %absorber moves left/right, so we may need this
        bds{end-1}.regridding.mode = 'linearize_uniform';
    end
    
    nodes = s.get_nodelist(); x = nodes(:,1);
    L = max(x) - min(x);
    %sets the xlim of the plot
    s.meta.plot_xlim = [0-L*0.02,L + L*0.02] + min(x);
    %number of processes in parallel pool
    s.meta.parallel_workers = parallel_workers;
    %fix time step to courant number
    s.stepping.courant_lock = 1;
    s.stepping.courant_target = CFL;

    s.stepping.dtmax = dtmax;
    

    % be careful with this, we may want to constantly update the absorbing
    % layer struct with the proper max(x), since the absorbing piston
    % drifts.
    if use_absorbing_layer
        h_ = abs(bds{end}.boundary_nodes(1,2));
        FS_add_absorbing_layer(s,1,max(x)-10, max(x), "h", h_, 'nu0',0.5);
    end
    if do_fileman
        fman = bem_sim_file_manager(fmandir);
        ncstore_add_autosave(s,fman,save_stagename,dt_store);
    end
else
    fprintf("Variable 's' already exists. Assuming this is a bem_sim" + ...
        " that we should be running. Terminate and clear if not.\n")
end


%s.meta.plot_xlim = [70,78];
%s.meta.plot_ylim = [0,2];

if do_FS_regrid
    s.boundaries{1}.regridding = struct( "mode","nodeshift",...
        "precede_update",0,...
        "order",0,...
        "regrid_type","shift_by_curve2d_normu", ...
        "a",16, ...
        "threshold_type","factor_global_local", ...
        "global_threshold", 10,...
        "local_threshold", 1.5,...
        "local_check_val",-1,...
        "global_check_val",-1,...
        "did_regrid",0, ...
        "startnode", 5,...
        "endpadding",5, ...
        "buffer_length", 10, ...
        "reweight_min",0 ...
    );
end

start_time = tic();
fignum = 0;
while s.stepping.t < TMAX
    s.full_step(); %only thing needed to step

    %UI
    if do_plot
        s.plot_full();
    end

    %LOGGING
    walltime = toc(start_time);
    wallsecs = mod(walltime,60);
    wallmins_ = (walltime - wallsecs)/60;
    wallmins = mod(wallmins_,60);
    wallhrs = (wallmins_ - wallmins)/60;
    logstr = sprintf('systime:%10.1fs (%02d:%02d:%06.3f): t=%10.4f (dt=%.6f, C=%.4f);\n',...
        walltime, wallhrs, wallmins, wallsecs, s.stepping.t, s.stepping.dt, s.stepping.courant_number);
    fprintf(logstr);

    %saving the plot, if we want to; make sure `savedir` is defined
    if do_plot && do_plotsave
        saveas(1,sprintf('%s/fig%d.png',savedir,fignum));
        fignum = fignum + 1;
    end
    
    %additional log info: regridding
    if log_regrids && do_FS_regrid
        fprintf("RG checks: global %f, local %f ", ...
            s.boundaries{1}.regridding.global_check_val,...
            s.boundaries{1}.regridding.local_check_val)
        if s.boundaries{1}.regridding.did_regrid
            fprintf("(Regrid done)")
        end
        fprintf("\n")
    end
end