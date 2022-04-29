classdef sim_netcdf_manager < handle
    %SIM_NETCDF_MANAGER Handles saving sim data to a netCDF file
    %   interfaces with a given netCDF file, buffering data and writing as
    %   necessary.
    
    properties
        filename
        buffer_size
        buffer_index
        buffer
        phase_info
        phase_bounds
    end
    
    methods
        function obj = sim_netcdf_manager(filename,overwrite)
            %SIM_NETCDF_MANAGER Construct an instance of this class
            %   Creates and initializes a netCDF file. If the file exists,
            %   then it is overwritten only if the overwrite flag is true.
            %   In any case, the file is populated with necessary
            %   attributes and variables.
            current_version = 1;
            
            obj.filename = filename;
            if overwrite && exist(filename, 'file')==2
              delete(filename);
            end
            
            if 0 == exist(filename, 'file')
                nccreate(filename,'phase_start','Dimensions',{'index',inf},...
                    'Format','netcdf4');
            end
            info = ncinfo(filename);
            if 0 == struct_arr_contains(info.Variables,'Name','phase_start')
                nccreate(filename,'phase_start','Dimensions',{'index',inf},...
                    'Format','netcdf4');
            end
            if 0 == struct_arr_contains(info.Variables,'Name','phase_end')
                nccreate(filename,'phase_end','Dimensions',{'index',inf},...
                    'Format','netcdf4');
            end
            if 0 == struct_arr_contains(info.Attributes,'Name','version')
                ncwriteatt(filename,'/','version',current_version)
            end
            
            obj.buffer_size = 10;
            obj.phase_info = [];
            obj.buffer = struct(...
                't',{},...
                'dt',{},...
                'x',{},...
                'z',{},...
                'c_num',{},...
                ...'K_n',{},...
                ...'K_d',{},...
                'phi',{},...
                'phi_n',{},...
                'energy',{},...
                'volume',{},...
                'mean_water_level',{},...
                'velocity_flux',{}...
            );
            obj.phase_info = struct(...
                'N',{},...
                'doublenode_indices',{},...
                'FS_start',{},...
                'FS_end',{},...
                'interp_M',{},...
                'interp_FS_M',{},...
                'interp_sliding_M',{},...
                'timestep_method',{},...
                'integral_spline',{}...
            );
            %append to phase_info relevant information
            i = 1; phasestr = sprintf('phase%d',i);
            while 0 ~= struct_arr_contains(info.Groups,'Name',phasestr)
                attr = info.Groups(i).Attributes;
                N_index = struct_arr_contains(info.Groups(i).Dimensions,'Name','N');
                obj.phase_info(i).N = info.Groups(i).Dimensions(N_index).Length;
                %obj.phase_info(i).doublenode_indices = ncread(filename,sprintf('%s/%s/doublenode_indices',filename,phasestr));
                obj.phase_info(i).doublenode_indices = attr(struct_arr_contains(attr,'Name','doublenode_indices')).Value;
                obj.phase_info(i).FS_start = attr(struct_arr_contains(attr,'Name','FS_start')).Value;
                obj.phase_info(i).FS_end = attr(struct_arr_contains(attr,'Name','FS_end')).Value;
                obj.phase_info(i).interp_M = attr(struct_arr_contains(attr,'Name','interp_M')).Value;
                obj.phase_info(i).interp_FS_M = attr(struct_arr_contains(attr,'Name','interp_FS_M')).Value;
                obj.phase_info(i).interp_sliding_M = attr(struct_arr_contains(attr,'Name','interp_sliding_M')).Value;
                obj.phase_info(i).timestep_method = attr(struct_arr_contains(attr,'Name','timestep_method')).Value;
                obj.phase_info(i).integral_spline = attr(struct_arr_contains(attr,'Name','integral_spline')).Value;
                i = i+1; phasestr = sprintf('phase%d',i);
            end
            obj.buffer_index = 1;

            ps = ncread(obj.filename,'phase_start');
            pe = ncread(obj.filename,'phase_end');
            obj.phase_bounds = zeros(length(ps),2);
            obj.phase_bounds(:,1) = ps;
            obj.phase_bounds(:,2) = pe;
        end
        
        function save_time(obj,sim)
            % Saves the current sim state to the file
            %Adds the current sim state as a snapshot to buffer for saving.
            %If the buffer is full, then it is flushed. If there is a phase
            %change, then a new phase is created and the buffer is flushed
            
            %TODO allow forcing a new phase
            sim = sim.update_characteristics();
            if ~obj.is_compatible_with_phase(sim,length(obj.phase_info))
                obj.flush_buffer();
                obj.append_phase(sim);
            end
            snapshot = struct(...
                't',sim.stepping.t,...
                'dt',sim.stepping.dt,...
                'x',sim.boundary.x,...
                'z',sim.boundary.z,...
                'c_num',sim.stepping.courant_number,...
                ...'K_n',sim.characteristics_.K_n,...
                ...'K_d',sim.characteristics_.K_d,...
                'phi',sim.characteristics_.phi,...
                'phi_n',sim.characteristics_.phi_n,...
                'energy',calc_energy(sim),...
                'volume',calc_volume(sim),...
                'mean_water_level',calc_mean_water_level(sim),...
                'velocity_flux',calc_total_flux(sim)...
            );
            obj.buffer(obj.buffer_index) = snapshot;
            obj.buffer_index = obj.buffer_index + 1;
            if obj.buffer_index > obj.buffer_size
                obj.flush_buffer();
            end
        end
        
        function append_phase(obj,sim)
            i = length(obj.phase_info)+1; phasestr = sprintf('phase%d',i);
            nccreate(obj.filename,sprintf('%s/t',phasestr),'Dimensions',{'t',inf},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/dt',phasestr),'Dimensions',{'t',inf},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/c_num',phasestr),'Dimensions',{'t',inf},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/x',phasestr),'Dimensions',{'t',inf,'N',sim.boundary.N},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/z',phasestr),'Dimensions',{'t',inf,'N',sim.boundary.N},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/phi',phasestr),'Dimensions',{'t',inf,'N',sim.boundary.N},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/phi_n',phasestr),'Dimensions',{'t',inf,'N',sim.boundary.N},'Format','netcdf4');
            %nccreate(obj.filename,sprintf('%s/K_n',phasestr),'Dimensions',{'t',inf,'N',sim.boundary.N,'N',sim.boundary.N},'Format','netcdf4');
            %nccreate(obj.filename,sprintf('%s/K_d',phasestr),'Dimensions',{'t',inf,'N',sim.boundary.N,'N',sim.boundary.N},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/energy',phasestr),'Dimensions',{'t',inf},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/volume',phasestr),'Dimensions',{'t',inf},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/mean_water_level',phasestr),'Dimensions',{'t',inf},'Format','netcdf4');
            nccreate(obj.filename,sprintf('%s/velocity_flux',phasestr),'Dimensions',{'t',inf},'Format','netcdf4');
            obj.phase_info(i).N = sim.boundary.N;
            ncwriteatt(obj.filename,phasestr,'doublenode_indices',find(sim.boundary.doublenode));
            obj.phase_info(i).doublenode_indices = find(sim.boundary.doublenode);
            ncwriteatt(obj.filename,phasestr,'FS_start',sim.boundary.FS_start);
            obj.phase_info(i).FS_start = sim.boundary.FS_start;
            ncwriteatt(obj.filename,phasestr,'FS_end',sim.boundary.FS_end);
            obj.phase_info(i).FS_end = sim.boundary.FS_end;
            ncwriteatt(obj.filename,phasestr,'interp_M',sim.interpolation.M);
            obj.phase_info(i).interp_M = sim.interpolation.M;
            ncwriteatt(obj.filename,phasestr,'interp_FS_M',sim.interpolation_FS.M);
            obj.phase_info(i).interp_FS_M = sim.interpolation_FS.M;
            ncwriteatt(obj.filename,phasestr,'interp_sliding_M',sim.interpolation_sliding.M);
            obj.phase_info(i).interp_sliding_M = sim.interpolation_sliding.M;
            ncwriteatt(obj.filename,phasestr,'timestep_method',sim.stepping.timestep_method);
            obj.phase_info(i).timestep_method = sim.stepping.timestep_method;
            ncwriteatt(obj.filename,phasestr,'integral_spline',sim.meta.integral_spline);
            obj.phase_info(i).integral_spline = sim.meta.integral_spline;
            
            ncwrite(obj.filename,'phase_start',sim.stepping.t,i);
            ncwrite(obj.filename,'phase_end',sim.stepping.t,i);
            obj.phase_bounds(i,:) = sim.stepping.t;
        end
        
        function compatibility = is_compatible_with_phase(obj,sim,phase_index)
            % returns if a sim is compatible with the phase (parameters same, and is at a later time)
            compatibility = 0;
            %phase index out of bounds
            if length(obj.phase_info) < phase_index || phase_index <= 0
                return
            end
            %time decreased? it must mean it was rolled back
            if any(size(obj.phase_bounds) < [phase_index,2]) || sim.stepping.t < obj.phase_bounds(phase_index,2)
                return
            end
            
            pinfo = obj.phase_info(phase_index);
            if pinfo.N ~= sim.boundary.N || pinfo.FS_start ~= sim.boundary.FS_start ||...
                    pinfo.FS_end ~= sim.boundary.FS_end ||...
                    pinfo.interp_M ~= sim.interpolation.M ||...
                    pinfo.interp_FS_M ~= sim.interpolation_FS.M ||...
                    pinfo.interp_sliding_M ~= sim.interpolation_sliding.M||...
                    (~strcmp(pinfo.timestep_method,sim.stepping.timestep_method))||...
                    pinfo.integral_spline ~= sim.meta.integral_spline
                return
            end
            
            %check doublenode equality
            dn1 = find(sim.boundary.doublenode);
            dn2 = pinfo.doublenode_indices;
            if length(dn1) ~= length(dn2)
                return
            end
            for i=1:length(dn1)
                if dn1(i) ~= dn2(i)
                    return
                end
            end
            compatibility = 1;
        end
        
        function flush_buffer(obj)
            buf = obj.buffer_index-1; %number of elements to flush
            if buf == 0
                return
            end
            phase = length(obj.phase_info); phasestr = sprintf('phase%d',phase);
            if phase == 0
                error('Illegal state. Attempting to flush with no phases')
            end
            N = obj.phase_info(phase).N;
            t = zeros(buf,1);
            dt = zeros(buf,1);
            c_num = zeros(buf,1);
            x = zeros(buf,N);
            z = zeros(buf,N);
            phi = zeros(buf,N);
            phi_n = zeros(buf,N);
            %K_n = zeros(buf,N,N);
            %K_d = zeros(buf,N,N);
            energy = zeros(buf,1);
            volume = zeros(buf,1);
            mean_water_level = zeros(buf,1);
            velocity_flux = zeros(buf,1);
            
            
            for i=1:buf
                t(i) = obj.buffer(i).t;
                dt(i) = obj.buffer(i).dt;
                c_num(i) = obj.buffer(i).c_num;
                x(i,:) = obj.buffer(i).x;
                z(i,:) = obj.buffer(i).z;
                phi_n(i,:) = obj.buffer(i).phi_n;
                phi(i,:) = obj.buffer(i).phi;
                %K_n(i,:,:) = obj.buffer(i).K_n;
                %K_d(i,:,:) = obj.buffer(i).K_d;
                energy(i) = obj.buffer(i).energy;
                volume(i) = obj.buffer(i).volume;
                mean_water_level(i) = obj.buffer(i).mean_water_level;
                velocity_flux(i) = obj.buffer(i).velocity_flux;
            end
            
            info = ncinfo(obj.filename);
            t_index = struct_arr_contains(info.Groups(phase).Dimensions,'Name','t');
            t_index = 1 + info.Groups(phase).Dimensions(t_index).Length;
            ncwrite(obj.filename,sprintf('%s/t',phasestr),t,t_index);
            ncwrite(obj.filename,sprintf('%s/dt',phasestr),dt,t_index);
            ncwrite(obj.filename,sprintf('%s/c_num',phasestr),c_num,t_index);
            ncwrite(obj.filename,sprintf('%s/x',phasestr),x,[t_index,1]);
            ncwrite(obj.filename,sprintf('%s/z',phasestr),z,[t_index,1]);
            ncwrite(obj.filename,sprintf('%s/phi',phasestr),phi,[t_index,1]);
            ncwrite(obj.filename,sprintf('%s/phi_n',phasestr),phi_n,[t_index,1]);
            %ncwrite(obj.filename,sprintf('%s/K_n',phasestr),K_n,[t_index,1,1]);
            %ncwrite(obj.filename,sprintf('%s/K_d',phasestr),K_d,[t_index,1,1]);
            ncwrite(obj.filename,sprintf('%s/energy',phasestr),energy,t_index);
            ncwrite(obj.filename,sprintf('%s/volume',phasestr),volume,t_index);
            ncwrite(obj.filename,sprintf('%s/mean_water_level',phasestr),mean_water_level,t_index);
            ncwrite(obj.filename,sprintf('%s/velocity_flux',phasestr),velocity_flux,t_index);
            
            %update phase bounds
            phase_times = ncread(obj.filename,sprintf('%s/t',phasestr));
            maxtime = max(phase_times);
            mintime = min(phase_times);
            ncwrite(obj.filename,'phase_start',mintime,phase);
            ncwrite(obj.filename,'phase_end',maxtime,phase);
            obj.phase_bounds(phase,1) = mintime;
            obj.phase_bounds(phase,2) = maxtime;
            obj.buffer_index = 1;
        end
        
        function [snapshot,phase] = get_snapshot(obj,t,phase)
            if ~exist('phase','var')
                phase = -1;
            end
            if phase <= 0 || phase > length(obj.phase_info)
                %get the phase with the closest time
                phase = 1; phasestr = sprintf('phase%d',phase);
                starts = ncread(obj.filename,'phase_start');
                ends = ncread(obj.filename,'phase_end');
                min_phase = 1;
                min_dist = inf;
                while phase <= length(obj.phase_info)
                    if (starts(phase) <= t && t <= ends(phase)) || abs(starts(phase) - t) < min_dist...
                            || abs(ends(phase) - t) < min_dist
                        %see if we have a closer time in this period
                        phaset = ncread(obj.filename,sprintf('%s/t',phasestr));
                        dist = min(abs(phaset - t));
                        if dist < min_dist
                            min_dist = dist;
                            min_phase = phase;
                        end
                    end
                    phase = phase + 1; phasestr = sprintf('phase%d',phase);
                end
                phase = min_phase;
            end
            %load the closest time of the phase
            phasestr = sprintf('phase%d',phase);
            phaset = ncread(obj.filename,sprintf('%s/t',phasestr));
            [min_dist,ind] = min(abs(phaset - t));
            if phase == length(obj.phase_info) && obj.buffer_index > 1
                %current phase, so also check the buffer
                buft = [obj.buffer(1:(obj.buffer_index-1)).t];
                [bufmin,bufind] = min(abs(buft - t));
                if isempty(phaset) || (bufmin < min_dist)
                    %closer time, so just return this one
                    snapshot = obj.buffer(bufind);
                    return
                end
            end
            
            
            snapshot = struct(...
                't',phaset(ind),...
                'dt',ncread(obj.filename,sprintf('%s/dt',phasestr),ind,1),...
                'x',ncread(obj.filename,sprintf('%s/x',phasestr),[ind,1],[1,inf]),...
                'z',ncread(obj.filename,sprintf('%s/z',phasestr),[ind,1],[1,inf]),...
                'c_num',ncread(obj.filename,sprintf('%s/c_num',phasestr),ind,1),...
                ...'K_n',ncread(obj.filename,sprintf('%s/K_n',phasestr),[ind,1,1],[1,inf,inf]),...
                ...'K_d',ncread(obj.filename,sprintf('%s/K_d',phasestr),[ind,1,1],[1,inf,inf]),...
                'phi',ncread(obj.filename,sprintf('%s/phi',phasestr),[ind,1],[1,inf]),...
                'phi_n',ncread(obj.filename,sprintf('%s/phi_n',phasestr),[ind,1],[1,inf]),...
                'energy',ncread(obj.filename,sprintf('%s/energy',phasestr),ind,1),...
                'volume',ncread(obj.filename,sprintf('%s/volume',phasestr),ind,1),...
                'mean_water_level',ncread(obj.filename,sprintf('%s/mean_water_level',phasestr),ind,1),...
                'velocity_flux',ncread(obj.filename,sprintf('%s/velocity_flux',phasestr),ind,1)...
            );
        end
        
        function [sim,phase] = get_sim(obj,t,phase)
            % loads a sim from the data contained in a snapshot
            if ~exist('phase','var')
                phase = -1;
            end
            [snapshot,phase] = get_snapshot(obj,t,phase);
            pinfo = obj.phase_info(phase);
            doublenode = zeros(1,pinfo.N); doublenode(pinfo.doublenode_indices) = 1;
            FS_i = pinfo.FS_start;
            FS_f = pinfo.FS_end;
            
            sim = sim_from_data(snapshot.x,snapshot.z,doublenode,snapshot.phi(FS_i:FS_f),...
                FS_i,FS_f);
            
            sim.meta.has_a0 = 0;
            sim.stepping.t = snapshot.t;
            sim.stepping.dt = snapshot.dt;
            sim.stepping.courant_number = snapshot.c_num;
            sim.stepping.timestep_method = pinfo.timestep_method;
            sim.interpolation = get_interpolation_struct(pinfo.interp_M,0,1);
            sim.interpolation_FS = get_interpolation_struct(pinfo.interp_M,0,pinfo.interp_FS_M);
            sim.interpolation_sliding = get_interpolation_struct(pinfo.interp_sliding_M,0,1);
        end
    end
end

function index = struct_arr_contains(st, fieldsearch, name)
    %searches st(i).(fieldsearch) to see if name is in there
    for i = 1:length(st)
        if strcmp(name,st(i).(fieldsearch))
            index = i;
            return
        end
    end
    index = 0;
end

