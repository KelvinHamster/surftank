if ~exist('do_fig','var')
    fprintf('No specified do_fig. defaulting to 0 (no plots)\n');
    do_fig = 0;
end

if ~exist('do_save','var')
    fprintf('No specified do_save. defaulting to 0 (no saving plots)\n');
    do_save = 0;
end

if ~exist('do_nc_export','var')
    fprintf('No specified do_nc_export. defaulting to 0 (no netcdf file)\n')
    do_nc_export = 0;
end

if ~exist('on_loop_callback','var')
    fprintf('no on_loop_callback specified. defaulting to identity map\n')
    on_loop_callback = @(s) s;
end

if ~exist('fignum','var')
    fignum = 0;
else
    if do_fig && do_save && fignum ~= 0
        fprintf('fignum defined. Assuming we do not start from plot 0.\n');
    end
end

if do_nc_export && ~exist('ncfile','var')
    error('NC exporting without "ncfile" specified! Please set the variable to a valid file!')
end

if ~exist('t_target','var')
    error('No target time "t_target" defined.')
end

if ~exist('s','var')
    error('No simulation "s" defined.')
end

if ~exist('dtmin','var')
    fprintf('No specified dtmin. defaulting to 10^(-6).\n')
    dtmin = 1e-6;
end
fprintf('The simulation will be halted if dt drops below %.5e\n',dtmin)

start_time = tic;
s = s.update_characteristics();
if do_nc_export && ~exist('manager','var')
    fprintf('no sim_netcdf_manager object as "manager". Creating a new one')
    fprintf(' without overwrite.\n')
    manager = sim_netcdf_manager(ncfile,0);
    manager.save_time(s);
end



if ~exist('logfile','var')
    fprintf('"logfile" variable not set. defaulting to "null".\n')
    logfile = 'null';
end

fprintf('logfile = "%s"\n    this is the file that is appended to for the print\n',logfile)
fprintf('    statements that follow. "null" means no logging to a file is done.\n');

max_pos2_init = max(s.boundary.x.^2 + s.boundary.z.^2);

%print parameters
logstr = sprintf('N = %d, N_FS = %d, interpolation.M = %d, interpolation_FS.M = %d, interpolation_sliding.M = %d\n',s.boundary.N, s.boundary.FS_end - s.boundary.FS_start + 1, s.interpolation.M, s.interpolation_FS.M, s.interpolation_sliding.M);
logstr = sprintf('%sintegral_spline = %d, g = %.5f, t_start = %.5f, method = %s, courant_lock = %d\n',logstr,s.meta.integral_spline,s.meta.g,s.stepping.t,s.stepping.timestep_method,s.stepping.courant_lock);
logstr = sprintf('%sinitialized dt = %.5f, initialized Courant number = %.5f\n',logstr,s.stepping.dt, s.stepping.courant_number);

fprintf('%s',logstr);
if ~strcmp(logfile,'null')
    fid = fopen(logfile,'a+');
    fprintf(fid,'%s',logstr);
    fclose(fid);
end


while s.stepping.t < t_target
    %=========[check blow-up criterion]===========
    max_pos2 = max(s.boundary.x.^2 + s.boundary.z.^2);
    %position size increases by a factor of 100
    if max_pos2 > (max_pos2_init*10000)
        logstr = sprintf('stopping criterion reached, largest node distance from the origin: %f\n',sqrt(max_pos2));
        fprintf('%s',logstr);
        if ~strcmp(logfile,'null')
            fid = fopen(logfile,'a+');
            fprintf(fid,'%s',logstr);
            fclose(fid);
        end
        break
    end
    
    if s.stepping.dt < dtmin
        logstr = sprintf('stopping criterion reached, dt < dtmin: %.5e < %.5e\n',s.stepping.dt, dtmin);
        fprintf('%s',logstr);
        if ~strcmp(logfile,'null')
            fid = fopen(logfile,'a+');
            fprintf(fid,'%s',logstr);
            fclose(fid);
        end
        break
    end
    
    %TODO energy/volume conservation criteria
    
    
    %=========[end blow-up criterion]===============
    s = s.step();
    s = s.handle_regrid(0);
    s = on_loop_callback(s);
    s = s.update_characteristics();
    if do_nc_export
        manager.save_time(s);
    end
    if do_fig
        s.plot_FS();
        if do_save
            saveas(1,sprintf('../figure/fig%d.png',fignum)); fignum = fignum + 1;
        end
    end
    walltime = toc(start_time);
    wallsecs = mod(walltime,60);
    wallmins_ = (walltime - wallsecs)/60;
    wallmins = mod(wallmins_,60);
    wallhrs = (wallmins_ - wallmins)/60;

    logstr = sprintf('systime:%10.1fs (%02d:%02d:%06.3f): t=%10.4f (dt=%.6f, C=%.4f, dxmin=%.4e, umax=%7f);\n',...
        walltime, wallhrs, wallmins, wallsecs, s.stepping.t, s.stepping.dt, s.stepping.courant_number, s.characteristics_.dxmin, s.characteristics_.umax);
    if s.regridding.did_regrid
        logstr = sprintf('%s    regrid performed at t=%.4f\n',logstr,s.stepping.t);
    end
    fprintf('%s',logstr);
    if ~strcmp(logfile,'null')
        fid = fopen(logfile,'a+');
        fprintf(fid,'%s',logstr);
        fclose(fid);
    end
end
if do_nc_export
    manager.flush_buffer();
end