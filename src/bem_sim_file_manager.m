classdef bem_sim_file_manager < handle
    %BEM_SIM_FILE_MANAGER Handles saving sim data to a netCDF file
    %   interfaces with a given netCDF file, buffering data and writing as
    %   necessary.
    
    properties
        sim_directory
        version %keep this for backwards compatibility, if ever necessary

        buffer_size
        buffer %buffer.stage.bdry
    end
    
    methods
        function obj = bem_sim_file_manager(directory, buffer_size)
            %BEM_SIM_FILE_MANAGER Construct an instance of this class
            % sets up the directory 'directory' to manage a simulator.
            % If the directory does not exist, it is created.

            if ~exist("buffer_size","var")
                buffer_size = 50;
            end
            obj.buffer_size = buffer_size;

            obj.sim_directory = directory;
            if ~exist(directory,'dir')
                mkdir(directory);
                obj.version = 2;
            end
            % if future format versions exist, insert some filechecking
            % scheme here for version seeking.
            obj.version = 2;
        end

        function save_time(obj,sim,stage,bd_mask, force_flush)
            % Saves the current sim state to the given stage (string).
            % a snapshot is created only for the boundaries given by
            % bd_mask, which should be a 1xn array of indices.

            if ~exist("force_flush","var")
                force_flush = 0;
            end

            if obj.push_buffer(sim,stage,bd_mask) || force_flush
                obj.flush_buffer(stage,bd_mask);
            end
        end

        function sim = get_sim(obj,stage,t)
            load(sim_filename(obj,stage),'sim')
            ncf = nc_filename(obj,stage);
            info = ncinfo(ncf);
            
            sim.stepping.t = t;
            for g=info.Groups
                bdst = g.Name;
                bd = ncread(ncf,sprintf('%s/bdry_index',bdst));
                bdry = sim.boundaries{bd};

                %load closest time into bdry
                T = ncread(ncf,sprintf('%s/t',bdst));
                if isempty(T)
                    %just use defaults
                    continue
                end
                [~,t_ind] = min(abs(T - t));

                bdry.boundary_nodes(:,1) = ...
                    ncread(ncf,sprintf('%s/x',bdst),[t_ind,1],[1,inf]);
                bdry.boundary_nodes(:,2) = ...
                    ncread(ncf,sprintf('%s/z',bdst),[t_ind,1],[1,inf]);
                bdry.characteristics.phi = ...
                    ncread(ncf,sprintf('%s/phi',bdst),[t_ind,1],[1,inf]);
                bdry.characteristics.phi_n = ...
                    ncread(ncf,sprintf('%s/phi_n',bdst),[t_ind,1],[1,inf]);
                bdry.characteristics.phi_t = ...
                    ncread(ncf,sprintf('%s/phi_t',bdst),[t_ind,1],[1,inf]);
                bdry.characteristics.phi_tn = ...
                   ncread(ncf,sprintf('%s/phi_tn',bdst),[t_ind,1],[1,inf]);
            end

            %are there function handles?
            if ~isempty(sim.global_update_events)
                fprintf("Global update events nonempty. We cannot" + ...
                    " ensure that those function handles are still" + ...
                    " valid after file-loading. Unless you know what" + ...
                    " you are doing, please redefine them.\n");
            end
            if ~isempty(sim.global_step_events)
                fprintf("Global step events nonempty. We cannot" + ...
                    " ensure that those function handles are still" + ...
                    " valid after file-loading. Unless you know what" + ...
                    " you are doing, please redefine them.\n");
            end
            for bd=1:length(sim.boundaries)
                if sim.boundaries{bd}.has_function_handles
                fprintf("Boundary %d probably has function handles." + ...
                    " We cannot" + ...
                    " ensure that those function handles are still" + ...
                    " valid after file-loading. Unless you know what" + ...
                    " you are doing, please redefine them.\n",bd);
                end
            end

        end

        function stages = get_stages(obj)
            files = {dir(obj.sim_directory).name};
            %check for strings s where files contain strings
            % s + '_dat.nc' and s + '_obj.mat'
            dats=replace(files(endsWith(files,'_dat.nc')),'_dat.nc','');
            objs=replace(files(endsWith(files,'_obj.mat')),'_obj.mat','');
            
            stages = {};
            k = 0;
            for f=dats
                if any(strcmp(objs,f))
                    k = k + 1;
                    stages(k) = f;
                end
            end
        end

        function flush_buffer(obj,stage,bd_mask)

            ncf = nc_filename(obj,stage);
            obj.ensure_stage(stage);

            for b = bd_mask
                bdgp = obj.ensure_bdgroup(stage,b);
                bdst = bd_str(obj,b);
                if ~isfield(obj.buffer.(stage),bdst)
                    continue
                end

                %get the time index to place this frame in
                t_index = struct_arr_contains(bdgp.Dimensions,'Name','t');
                if t_index
                    t_index = 1 + bdgp.Dimensions(t_index).Length;
                else
                    %dim was not set yet; just populate the first slot.
                    t_index = 1;
                end
                tlen = obj.buffer.(stage).(bdst).buf_len;
                if tlen < 1
                    continue
                end
    
                ncwrite(ncf,sprintf('%s/t',bdst), ...
                    obj.buffer.(stage).(bdst).t(1:tlen),t_index);
                ncwrite(ncf,sprintf('%s/x',bdst), ...
                    obj.buffer.(stage).(bdst).x(1:tlen,:),[t_index,1]);
                ncwrite(ncf,sprintf('%s/z',bdst), ...
                    obj.buffer.(stage).(bdst).z(1:tlen,:),[t_index,1]);
                ncwrite(ncf,sprintf('%s/phi',bdst), ...
                    obj.buffer.(stage).(bdst).phi(1:tlen,:),[t_index,1]);
                ncwrite(ncf,sprintf('%s/phi_n',bdst), ...
                    obj.buffer.(stage).(bdst).phi_n(1:tlen,:),[t_index,1]);
                ncwrite(ncf,sprintf('%s/phi_t',bdst), ...
                    obj.buffer.(stage).(bdst).phi_t(1:tlen,:),[t_index,1]);
                ncwrite(ncf,sprintf('%s/phi_tn',bdst), ...
                    obj.buffer.(stage).(bdst).phi_tn(1:tlen,:),[t_index,1]);

                obj.buffer.(stage).(bdst).buf_len = 0;
            end
        end

        %TODO:
        % compatibility fcn (is this sim compatible with a given stage?)
    end

    methods (Access=protected)

        function ensure_stage(obj,stagename,sim)
            %makes sure the stage files are properly set/configured
            %a reference to the bem_sim object is needed, which can be
            %either passed through the second argument `sim`, or pulled
            %from the stage's buffer structure.
            %If the buffer is not configured, configure it as well. In such
            %a case, `sim` should be given as an argument.

            if ~exist("sim","var")
                sim = obj.buffer.(stagename).sim_obj_ref_;
            end
            
            sf = sim_filename(obj,stagename);
            if ~exist(sf,'file')
                save(sf,"sim");
            end
            %TODO if sim file exists: check compatibility
            
            ncf = nc_filename(obj,stagename);
            if ~exist(ncf,'file')
                nccreate(ncf,'version','Format','netcdf4');
                ncwrite(ncf,'version',obj.version);
            end

            if ~isfield(obj.buffer,stagename)
                obj.buffer.(stagename) = struct();
                obj.buffer.(stagename).sim_obj_ref_ = sim;
            end
        end

        function bdgp = ensure_bdgroup(obj,stage,bd,sim)
            %similar to ensure_stage, except, relevant to the .nc file,
            %prepares the group information.
            ncf = nc_filename(obj,stage);
            info = ncinfo(ncf);
            bdst = bd_str(obj,bd);

            if exist("sim","var")
                N = sim.boundaries{bd}.node_count;
            else
                N = length(obj.buffer.(stage).(bdst).x(1,:));
            end
            gid = struct_arr_contains(info.Groups,'Name',bdst);
            if gid %does the bdry group exist?
                %check compatibility
                bdgp = info.Groups(gid);

                %node count equal?
                n_ind = struct_arr_contains(bdgp.Dimensions,...
                    'Name','N');
                if n_ind > 0 && ...
                        bdgp.Dimensions(n_ind).Length ~= N
                    error("Node count for boundary %d changed in " + ...
                        "stage '%s'",bd,stage);
                end
            else
                nccreate(ncf,sprintf('%s/bdry_index',bdst), ...
                    'Format','netcdf4');
                ncwrite(ncf,sprintf('%s/bdry_index',bdst), ...
                    bd,1);
                nccreate(ncf,sprintf('%s/t',bdst), ...
                    'Dimensions',{'t',inf},'Format','netcdf4');
                nccreate(ncf,sprintf('%s/x',bdst),'Dimensions', ...
                    {'t',inf,'N',N},'Format','netcdf4');
                nccreate(ncf,sprintf('%s/z',bdst),'Dimensions', ...
                    {'t',inf,'N',N},'Format','netcdf4');
                nccreate(ncf,sprintf('%s/phi',bdst),'Dimensions', ...
                    {'t',inf,'N',N},'Format','netcdf4');
                nccreate(ncf,sprintf('%s/phi_n',bdst),'Dimensions', ...
                    {'t',inf,'N',N},'Format','netcdf4');
                nccreate(ncf,sprintf('%s/phi_t',bdst),'Dimensions', ...
                    {'t',inf,'N',N},'Format','netcdf4');
                nccreate(ncf,sprintf('%s/phi_tn',bdst),'Dimensions', ...
                    {'t',inf,'N',N},'Format','netcdf4');

                info = ncinfo(ncf);
                gid = struct_arr_contains(info.Groups,'Name',bdst);
                bdgp = info.Groups(gid);
            end
        end

        

        function should_flush = push_buffer(obj,sim,stage,bd_mask)
            %helper function: when saving, we want to push stuff onto the
            %buffer. This does that. Once complete, a value is returned for
            %whether or not a flush should be called, which should be
            %whenever a buffer gets full.

            sim.update_characteristics();
            obj.ensure_stage(stage,sim);
            
            should_flush = 0;
            for b = bd_mask
                bdry = sim.boundaries{b};
                bdst = bd_str(obj,b);
                if ~isfield(obj.buffer.(stage),bdst)
                    obj.gen_empty_bd_buffer(stage, b, bdry.node_count);
                end
                t_ind = min(obj.buffer.(stage).(bdst).buf_len + 1,...
                    obj.buffer_size);
                
                obj.buffer.(stage).(bdst).t(t_ind) = sim.stepping.t;
                obj.buffer.(stage).(bdst).x(t_ind,:) = ...
                    bdry.boundary_nodes(:,1)';
                obj.buffer.(stage).(bdst).z(t_ind,:) = ...
                    bdry.boundary_nodes(:,2)';
                obj.buffer.(stage).(bdst).phi(t_ind,:) = ...
                    bdry.characteristics.phi;
                obj.buffer.(stage).(bdst).phi_n(t_ind,:) = ...
                    bdry.characteristics.phi_n;
                obj.buffer.(stage).(bdst).phi_t(t_ind,:) = ...
                    bdry.characteristics.phi_t;
                obj.buffer.(stage).(bdst).phi_tn(t_ind,:) = ...
                    bdry.characteristics.phi_tn;
    
                obj.buffer.(stage).(bdst).buf_len = t_ind;
                if t_ind == obj.buffer_size
                    should_flush = 1;
                end
            end
        end



        function gen_empty_bd_buffer(obj,stage,boundary,num_nodes)
            obj.buffer.(stage).(bd_str(obj,boundary)) = struct( ...
                "buf_len",0, ...
                "t",zeros(obj.buffer_size,1),...
                "x",zeros(obj.buffer_size,num_nodes),...
                "z",zeros(obj.buffer_size,num_nodes),...
                "phi",zeros(obj.buffer_size,num_nodes),...
                "phi_n",zeros(obj.buffer_size,num_nodes),...
                "phi_t",zeros(obj.buffer_size,num_nodes),...
                "phi_tn",zeros(obj.buffer_size,num_nodes)...
                );
        end
    end
end

function s = bd_str(obj,bdry)
    if obj.version == 2
    s = sprintf("bd%d",bdry);
    return
    end
    s = sprintf("bd%d",bdry);
end

function fname = nc_filename(obj,stagename)
    fname = fullfile(obj.sim_directory,...
                sprintf("%s_dat.nc",stagename));
end
function fname = sim_filename(obj,stagename)
    fname = fullfile(obj.sim_directory,...
                sprintf("%s_obj.mat",stagename));
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
