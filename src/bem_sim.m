classdef bem_sim
    %Manages a simulator
    %   Detailed explanation goes here
    
    properties
        boundary
        interpolation
        interpolation_FS
        interpolation_sliding
        stepping
        history
        meta
        characteristics_
        regridding
    end
    
    methods
        function obj = bem_sim(varargin)
            sim_types = {'standing_wave', 'soliton_reflect','empty'};
            step_types = {'taylor2', 'abm4','euler'};
            
            %default quantities
            type = 'standing_wave';
            L = 'default';
            a0 = 'default';
            Nx = 100;
            h = 1;
            dt = 0.01;
            store_dt = 0.25;
            g = 9.81;
            parallel_workers = 0;
            timestep_method = 'taylor2';
            wall_resolution_factor = 1;
            history_store_full_boundary = 1;
            history_store_angle_FS = 0;
            x0 = 20;
            plot_xlim = 'default';
            plot_ylim = 'default';
            integral_spline = 0;
            %matfile = 'null';
            
            p = inputParser;
            addOptional(p,'type',type, @(x) any(validatestring(x,sim_types)));
            addOptional(p,'timestep_method',timestep_method, @(x) any(validatestring(x,step_types)));
            addOptional(p,'L',L, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
            addOptional(p,'a0',a0, @(x) (isstring(x) && strcmp(x,'default')) || ( isnumeric(x) && isscalar(x) && (x > 0) ) );
            addOptional(p,'Nx',Nx, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1) == 0) );
            addOptional(p,'h',h, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
            addOptional(p,'dt',dt, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
            addOptional(p,'store_dt',store_dt, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
            addOptional(p,'g',g, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
            addOptional(p,'parallel_workers',parallel_workers, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1) == 0) );
            addOptional(p,'wall_resolution_factor',wall_resolution_factor, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
            addOptional(p,'history_store_full_boundary',history_store_full_boundary,@(x) islogical(x) || (x==0) || (x==1));
            addOptional(p,'history_store_angle_FS',history_store_angle_FS,@(x) islogical(x) || (x==0) || (x==1));
            addOptional(p,'x0',x0, @(x) isnumeric(x) && isscalar(x) );
            addOptional(p,'plot_xlim',plot_xlim, @(x) (isnumeric(x) && length(x)==2 && (x(1) < x(2))) || (isstring(x) && strcmp(x,'default')) );
            addOptional(p,'plot_ylim',plot_ylim, @(x) (isnumeric(x) && length(x)==2 && (x(1) < x(2))) || (isstring(x) && strcmp(x,'default')) );
            %addOptional(p,'matfile',matfile,@(x) ischar(x) || isstring(x));
            
            parse(p, varargin{:});
            
            type = p.Results.type;
            timestep_method = p.Results.timestep_method;
            L = p.Results.L;
            a0 = p.Results.a0;
            Nx = p.Results.Nx;
            h = p.Results.h;
            dt = p.Results.dt;
            store_dt = p.Results.store_dt;
            g = p.Results.g;
            parallel_workers = p.Results.parallel_workers;
            wall_res_factor = p.Results.wall_resolution_factor;
            history_store_full_boundary = p.Results.history_store_full_boundary;
            history_store_angle_FS = p.Results.history_store_angle_FS;
            x0 = p.Results.x0;
            plot_xlim = p.Results.plot_xlim;
            plot_ylim = p.Results.plot_ylim;
            %matfile = p.Results.matfile;
            
            has_a0 =  ~(isstring(a0) && strcmp(a0,'default'));
            
            %Construct an instance of this class
            %   Detailed explanation goes here
            obj.boundary = struct('N',0,'x',[],'z',[],'doublenode',[],...
                'FS_start',-1,'FS_end',-1,'phi_FS',[]);
            
            %TODO delegate interpolation params to a helper function
            
            M = 1;
            Lmat = [1, -1; 0, 1];     
            Lpmat = [-1;1];      
            Lppmat = zeros(2,0);
            obj.interpolation = struct('M',M, 'Lagrange', Lmat, 'Lagrange_prime', Lpmat, 'Lagrange_primeprime',Lppmat);
            obj.interpolation_FS = struct('M',M, 'Lagrange', Lmat, 'Lagrange_prime', Lpmat, 'Lagrange_primeprime',Lppmat);
            obj.interpolation_sliding = get_interpolation_struct(4, 0, 1);

            obj.stepping = struct('t',0, 'dt',dt, 'store_dt', store_dt, 'timestep_method', timestep_method, ...
                'courant_lock',1,'courant_number',0.5);
            obj.history = struct('time',[],'energy',[],'volume',[],'KE',[],'PE',[],...
                'net_flux',[],'mean_water_level',[], 'condition_number',[]);
            obj.history.x_FS = {};
            obj.history.z_FS = {};
            obj.history.x_full = {};
            obj.history.z_full = {};
            obj.history.FS_start = {};
            obj.history.FS_end = {};
            obj.history.phi_FS = {};
            obj.history.beta_FS = {};
            switch type
                case 'standing_wave'
                    %standing_wave
                    init_standing_wave
                    
                    N = length(x) - 1;
                    obj.boundary.N = N;
                    obj.boundary.x = x(1:N);
                    obj.boundary.z = z(1:N);
                    obj.boundary.doublenode = doublenode;
                    obj.boundary.FS_start = FS_i;
                    obj.boundary.FS_end = FS_f;
                    obj.boundary.phi_FS = phi_FS;
                    obj.meta.has_a0 = 1;
                    obj.meta.a0 = a0;
                    obj.meta.omega = omega;
                case 'soliton_reflect'
                    init_soliton_reflect
                    
                    N = length(x) - 1;
                    obj.boundary.N = N;
                    obj.boundary.x = x(1:N);
                    obj.boundary.z = z(1:N);
                    obj.boundary.doublenode = doublenode;
                    obj.boundary.FS_start = FS_i;
                    obj.boundary.FS_end = FS_f;
                    obj.boundary.phi_FS = phi_FS;
                    obj.meta.has_a0 = 1;
                    obj.meta.a0 = a0;
                case 'empty'
                    obj.meta.has_a0 = has_a0;
                    obj.meta.a0 = a0;
            end
            if strcmp(plot_ylim,'default')
                obj.meta.plot_ylim = [-1.5,1.5];
            else
                obj.meta.plot_ylim = plot_ylim;
            end
            if isempty(obj.boundary.x)
                minx = 0;
                maxx = 1;
            else
                minx = min(obj.boundary.x);
                maxx = max(obj.boundary.x);
            end
            if strcmp(plot_xlim,'default')
                obj.meta.plot_xlim = [minx-(maxx-minx)/20,maxx+(maxx-minx)/20];
            else
                obj.meta.plot_xlim = plot_xlim;
            end
            obj.meta.sim_parallel_workers = parallel_workers;
            obj.meta.g = g;
            %obj.meta.matfile = matfile;
            if strcmp(L,'default')
                obj.meta.L = maxx - minx;
            else
                obj.meta.L = L;
            end
            obj.meta.wall_resolution_factor = wall_res_factor;
            obj.meta.history_store_full_boundary = history_store_full_boundary;
            obj.meta.history_store_angle_FS = history_store_angle_FS;
            obj.meta.integral_spline = integral_spline;
            obj.characteristics_ = struct('K_n',[],'K_d',[],'COEFS_n',[],'cosbeta',[],'sinbeta',[],'ctime','null','phi',[],'phi_n',[],'mat_A',[],'mat_b',[],'FS_spline_x',0,'FS_spline_z',0);
            obj.regridding = struct('regrid_type','none','threshold_type','manual','did_regrid',0);
        end
        
        function obj = update_characteristics(obj)
            %Ensures that the variables in characteristics_ match the current boundary values
            if (ischar(obj.characteristics_.ctime) && strcmp(obj.characteristics_.ctime,'null'))||(obj.characteristics_.ctime ~= obj.stepping.t)
                obj = eval_Kn_Kd(obj);
            end
        end
        
        function obj = step(obj)
            %Steps the simulation forward by dt, updating the free surface.
            % defers to step_sim.m
            
            obj = step_sim(obj);
        end
        
        function obj = update_history(obj)
            %Stores values in history when enough time has passed.
            
            
            
            %if we do not have a starting value in history, or enough time
            %has passed, append to the history.
            histlength = length(obj.history.time);
            if histlength == 0 || obj.stepping.t - obj.history.time(histlength) >= obj.stepping.store_dt - (10^-12)
                obj = obj.update_characteristics();
                
                x = [obj.boundary.x obj.boundary.x(1)];
                z = [obj.boundary.z obj.boundary.x(1)];
                
                FS_i = obj.boundary.FS_start;
                FS_f = obj.boundary.FS_end;
                
                beta = atan2(obj.characteristics_.sinbeta(FS_i:FS_f),obj.characteristics_.cosbeta(FS_i:FS_f));
                
                [E,ke,pe] = calc_energy(obj);
                obj = obj.append_history(obj.stepping.t,E,ke,pe,calc_volume(obj),x(FS_i:FS_f),z(FS_i:FS_f),...
                    obj.boundary.phi_FS,calc_total_flux(obj),calc_mean_water_level(obj),cond(obj.characteristics_.mat_A),...
                    obj.boundary.x,obj.boundary.z, FS_i, FS_f, beta);
                
            end
        end
        
        function obj = clear_history(obj)
            obj.history = struct('time',[],'energy',[],'volume',[],'KE',[],'PE',[],...
                'net_flux',[],'mean_water_level',[], 'condition_number',[]);
            obj.history.x_FS = {};
            obj.history.z_FS = {};
            obj.history.x_full = {};
            obj.history.z_full = {};
            obj.history.FS_start = {};
            obj.history.FS_end = {};
            obj.history.phi_FS = {};
            obj.history.beta_FS = {}; 
        end
        
        function mat_b = gen_mat_b(obj, bt, bv)
            N = obj.boundary.N;
            gD = find(~bt(1:N));
            gN = find(bt(1:N));
            
            mat_b = sum(bv(gN).*obj.characteristics_.K_d(:,gN), 2) - sum(bv(gD).*obj.characteristics_.COEFS_n(:,gD), 2);

            %the following is copy-paste of the end of the BEM code====================
            for i=find(obj.boundary.doublenode)
                ip1 = i + 1 - N*(i >= N);

                %if a node has a neumann condition, use phi_i = phi_{i+1}
                if bt(i)
                    if bt(ip1)                                     %neumann-neumann
                        mat_b(i) = 0;
                    else                                           %neumann-dirichlet
                        mat_b(i) = bv(ip1);
                    end
                else
                    %no dirichlet-dirichlet boundaries so default to
                                                                   %dirichlet-neumann
                        mat_b(ip1) = bv(i);
                end
            end
        end
        
        function obj = append_history(obj,t,e,ke,pe,volume,x_FS,z_FS,phi_FS,net_flux,mean_water_level,condition_number,...
                x_full,z_full, fs_start, fs_end,beta)
            %Appends the given values to their corresponding arrays in history.
            store_full = obj.meta.history_store_full_boundary;
            obj.history.time = [obj.history.time t];
            obj.history.energy = [obj.history.energy e];
            obj.history.KE = [obj.history.KE ke];
            obj.history.PE = [obj.history.PE pe];
            obj.history.volume = [obj.history.volume volume];
            index = length(obj.history.time);
            if store_full
                obj.history.x_full{index} = x_full;
                obj.history.z_full{index} = z_full;
            else
                obj.history.x_FS{index} = x_FS;
                obj.history.z_FS{index} = z_FS;
            end
            
            if obj.meta.history_store_angle_FS
                obj.history.beta_FS{index} = beta;
            end
            obj.history.FS_start{index} = fs_start;
            obj.history.FS_end{index} = fs_end;
            obj.history.phi_FS{index} = phi_FS;
            obj.history.net_flux = [obj.history.net_flux net_flux];
            obj.history.mean_water_level = [obj.history.mean_water_level mean_water_level];
            obj.history.condition_number = [obj.history.condition_number condition_number];
        end
        
        function phi = interior_eval(obj,x,z)
            %Evaluates phi on an interior point (x,z)
            obj = obj.update_characteristics();
            
            phi = bem_integrate(obj,obj.characteristics_.phi_n,3,x,z) - bem_integrate(obj,obj.characteristics_.phi,4,x,z);
        end
        
        function result = is_point_inside(obj,xp,zp)
            %returns whether or not (xp,zp) is on the interior of the boundary
            x = obj.boundary.x; z = obj.boundary.z;
            N = obj.boundary.N;
            
            %even odd rule
            %(https://en.wikipedia.org/wiki/Even-odd_rule)
            j = N;
            result = false;
            
            for i=1:N
                
                if (z(i) > zp) ~= (z(j) > zp)
                    %crossed z = zp line
                    dotnorm = (xp-x(i))*(z(j)-z(i)) - (zp-z(i))*(x(j)-x(i));
                    
                    
                    if (dotnorm < 0) ~= (z(j) < z(i))
                        %z(j) < z(i) iff norm is pointing to +x
                        %so we are here when (xp,zp) is on the left side of
                        %the line segment
                        result = ~result;
                    end
                end
                
                j = i;
            end
            
        end
        
        function plot_FS(obj,omega)
            %Plots the free surface in figure 1.
            %   If omega is given as an optional argument, the title will
            %   show the fraction of a period the simulation is currently at.
            %
            have_a0 = obj.meta.has_a0;
            if have_a0
                a0 = obj.meta.a0;
            else
                a0 = 1;
            end
            
            FS_i = obj.boundary.FS_start;
            FS_f = obj.boundary.FS_end;
            t = obj.stepping.t;
            
            %plot free surface boundary and nodes
            set(0,'CurrentFigure',1);
            plot(obj.boundary.x(FS_i:FS_f),obj.boundary.z(FS_i:FS_f)./a0)
            hold on
            scatter(obj.boundary.x(FS_i:FS_f),obj.boundary.z(FS_i:FS_f)./a0,'o')
            
            
            ylim(obj.meta.plot_ylim)
            if have_a0
                ylabel('eta/a0')
            else
                ylabel('eta')
            end
            xlim(obj.meta.plot_xlim)
            xlabel('x')
            
            if ~exist('omega','var') && isfield(obj.meta,'omega')
                omega = obj.meta.omega;
            end
            
            if exist('omega','var')
                title(sprintf('(w*t)/(2pi) = %.2f', omega*t/(2*pi)))
            else
                title(sprintf('t = %.4f', t))
            end
            hold off
            drawnow
            
        end
        
        function plot_full(obj,omega)
            %Plots the full boundary in figure 1.
            %   If omega is given as an optional argument, the title will
            %   show the fraction of a period the simulation is currently at.
            %
            
            t = obj.stepping.t;
            
            %plot free surface boundary and nodes
            set(0,'CurrentFigure',1);            
            plot(obj.boundary.x,obj.boundary.z)
            hold on
            scatter(obj.boundary.x,obj.boundary.z,'o')
            
            
            ylabel('z')
            xlim(obj.meta.plot_xlim)
            xlabel('x')
            
            if ~exist('omega','var') && isfield(obj.meta,'omega')
                omega = obj.meta.omega;
            end
            
            if exist('omega','var')
                title(sprintf('(w*t)/(2pi) = %.2f', omega*t/(2*pi)))
            else
                title(sprintf('t = %.4f', t))
            end
            hold off
            drawnow
            
        end
        
        function new_obj = rollback(obj,t)
            %returns a simulation rolled back to a history closest to time t
            %
            % non-free-surface boundaries are expected to not change, and
            % double nodes on the free surface are expected to only be on
            % the ends.
            if isempty(obj.history.time)
                new_obj = obj;
                return
            end
            [~,t_i] = min(abs(obj.history.time - t));
            new_obj = bem_sim('type', 'empty');
            new_obj.boundary = obj.boundary;
            new_obj.interpolation = obj.interpolation;
            new_obj.interpolation_FS = obj.interpolation_FS;
            new_obj.interpolation_sliding = obj.interpolation_sliding;
            new_obj.stepping = obj.stepping;
            new_obj.history = obj.history;
            new_obj.meta = obj.meta;
            
            new_obj.stepping.t = obj.history.time(t_i);
            
            %=== update boundary ===
            if obj.meta.history_store_full_boundary
                new_obj.boundary.x = obj.history.x_full{t_i};
                new_obj.boundary.z = obj.history.z_full{t_i};
            else
                new_obj.boundary.x = [obj.boundary.x(1:obj.boundary.FS_start - 1) obj.history.x_FS{t_i} obj.boundary.x(obj.boundary.FS_end + 1:obj.boundary.N)];
                new_obj.boundary.z = [obj.boundary.z(1:obj.boundary.FS_start - 1) obj.history.z_FS{t_i} obj.boundary.z(obj.boundary.FS_end + 1:obj.boundary.N)];
                FS_im1 = obj.boundary.FS_start - 1;
                FS_fp1 = obj.boundary.FS_end + 1;
                if FS_im1 == 0
                    FS_im1 = length(new_obj.boundary.x);
                elseif FS_fp1 > length(new_obj.boundary.x)
                    FS_fp1 = 1;
                end
              
                new_obj.boundary.x(FS_im1) = new_obj.boundary.x(obj.boundary.FS_start);
                new_obj.boundary.x(FS_im1) = new_obj.boundary.x(obj.boundary.FS_start);
                new_obj.boundary.x(FS_fp1) = new_obj.boundary.x(obj.boundary.FS_end);
                new_obj.boundary.z(FS_fp1) = new_obj.boundary.z(obj.boundary.FS_end);
            end
            new_obj.boundary.FS_start = obj.history.FS_start{t_i};
            new_obj.boundary.FS_end = obj.history.FS_end{t_i};
            new_obj.boundary.phi_FS = obj.history.phi_FS{t_i};
            new_obj.boundary.doublenode = (new_obj.boundary.x == circshift(new_obj.boundary.x,-1)) ...
                & (new_obj.boundary.z == circshift(new_obj.boundary.z,-1));
            %=== update history ===
            
            if obj.meta.history_store_full_boundary
                new_obj.history.x_full = new_obj.history.x_full(1:t_i);
                new_obj.history.z_full = new_obj.history.z_full(1:t_i);
            else
                new_obj.history.x_FS = new_obj.history.x_FS(1:t_i);
                new_obj.history.z_FS = new_obj.history.z_FS(1:t_i);
            end
            new_obj.history.time = new_obj.history.time(1:t_i);
            new_obj.history.energy = new_obj.history.energy(1:t_i);
            new_obj.history.KE = new_obj.history.KE(1:t_i);
            new_obj.history.PE = new_obj.history.PE(1:t_i);
            new_obj.history.volume = new_obj.history.volume(1:t_i);
            new_obj.history.FS_start = new_obj.history.FS_start(1:t_i);
            new_obj.history.FS_end = new_obj.history.FS_end(1:t_i);
            new_obj.history.phi_FS = new_obj.history.phi_FS(1:t_i);
            new_obj.history.net_flux = new_obj.history.net_flux(1:t_i);
            new_obj.history.mean_water_level = new_obj.history.mean_water_level(1:t_i);
            new_obj.history.condition_number = new_obj.history.condition_number(1:t_i);
            
        end
        
        function obj = remove_free_surface_nodes(obj, nodes)
            %removes free surface nodes at the specified indices
            %
            %   nodes is an array of integers, where k = nodes(i) requests
            %   the node indexed by FS_start + k to be removed. So nodes
            %   should take values between 1 and FS_end-FS_start-1
            %
            nodes = sort(nodes(:),1,'descend')';
            
            FS_i = obj.boundary.FS_start;
            FS_f = obj.boundary.FS_end;
            
            nodes = nodes + FS_i - 1;
            %remove invalid indices
            nodes = nodes((nodes < FS_f) & (nodes > FS_i));
            
            
            obj.boundary.x(nodes) = [];
            obj.boundary.z(nodes) = [];
            obj.boundary.doublenode(nodes) = [];
            obj.boundary.phi_FS(nodes - FS_i + 1) = [];

            obj.boundary.FS_end = obj.boundary.FS_start + length(obj.boundary.phi_FS) - 1;
            obj.boundary.N = length(obj.boundary.x);
            
            obj.characteristics_.ctime = 'null'; %force update next time we need characteristics_
            if isfield(obj.characteristics_,'time_derivs')
                obj.characteristics_ = rmfield(obj.characteristics_,'time_derivs');
            end
        end
        
        function obj = add_free_surface_nodes(obj, xi, tol)
            %adds new free surface nodes at the specified parameterization indices, interpolating as necessary.
            %
            %   a node is placed at each value of xi. For each xi(i), set
            %   k+s = xi(i) so that k is an integer and 0 <= s < 1. The
            %   node is added in segment k of the free surface (between
            %   nodes FS_start+k and FS_start+k+1) offset by s, which
            %   uniformly distributes the node onto the segment by arc
            %   length (s = 0 corresponds to node FS_start+k and s = 0.5
            %   corresponds to the midpoint of segment k).
            %
            %   tol can be specified as an optional argument (default 10^-3)
            %   and is the allowed deviation s from its actual value in node
            %   placement. The bisection method is used, so tol specifies
            %   the maximum domain size before termination.
            %
            %   this method updates characteristics_ prior to addition and
            %   sets characteristics_.ctime to 'null' so the next
            %   characteristics update will automatically call eval_Kn_Kd
            
            if ~exist('tol','var')
                tol = 0.001;
            end
            
            tol = tol * 2;%we will use the midpoint, so the search interval can be shrunk to [c - tol,c + tol], length of 2tol
            
            QUAD_T = [0.0052995325041750333469603 0.0277124884633837102743126 0.0671843988060841224019271 0.1222977958224984867952045 0.1910618777986781147149031 0.2709916111713863151599924 0.3591982246103705422868302 0.4524937450811812866824368 0.5475062549188187688287144 0.6408017753896294577131698 0.7290083888286137403511589 0.8089381222013218852850969 0.8777022041775015548381589 0.9328156011939158220869217 0.9722875115366163001340283 0.9947004674958249692551249]';
            QUAD_W = [0.0135762297058770482066636 0.0311267619693239468159351 0.0475792558412463928441127 0.0623144856277669384470030 0.0747979944082883679845608 0.0845782596975012679330064 0.0913017075224617918882686 0.0947253052275342510846201 0.0947253052275342510846201 0.0913017075224617918882686 0.0845782596975012679330064 0.0747979944082883679845608 0.0623144856277669384470030 0.0475792558412463928441127 0.0311267619693239468159351 0.0135762297058770482066636]';

            
            obj = obj.update_characteristics();
            
            obj.characteristics_.ctime = 'null'; %force update next time we need characteristics_
            if isfield(obj.characteristics_,'time_derivs')
                obj.characteristics_ = rmfield(obj.characteristics_,'time_derivs');
            end

            FS_i = obj.boundary.FS_start;
            FS_f = obj.boundary.FS_end;
            
            M_FS = obj.interpolation_FS.M;
            Lmat_FS = obj.interpolation_FS.Lagrange;
            
            phi_FS = obj.boundary.phi_FS;
            
            %work in descending order to avoid splicing conflicts
            xi = sort(xi(:),1,'descend')';
            
            for t = xi
                s = mod(t,1);
                k = round(t - s) + 1;
                ax = obj.characteristics_.FS_spline_x(k,:);
                az = obj.characteristics_.FS_spline_z(k,:);
                % r_x(t) = dot(ax,[t^3 t^2 t 1])
                
                axp = ax(1:3) .* (3:-1:1);
                azp = az(1:3) .* (3:-1:1);
                % r_x'(t) = dot(ax(1:3),[3t^2 2t 1])
                
                len = dot(vecnorm([sum(axp.* QUAD_T.^(2:-1:0),2) sum(azp.* QUAD_T.^(2:-1:0),2)],2,2) , QUAD_W);
                slen = s * len;
                %f(t) = int(|r'|)_0^t - s*len
                f = @(t) t*dot(vecnorm([sum(axp.* (QUAD_T.*t).^(2:-1:0),2) sum(azp.* (QUAD_T.*t).^(2:-1:0),2)],2,2) , QUAD_W) - slen;
                
                
                a = 0; b = 1; c = 0.5;
                %we do not need to check if a zero is between a and b, there is one, and its unique because f is strictly increasing
                fa = -slen; %fb = (1-s)*len;
                %binary search until we are within x tolerance
                while b-a > tol
                    fc = f(c);
                    if fa*fc <= 0
                        %f(a),f(c) have different signs, zero is between a and c.
                        b = c; %fb = fc;
                    else
                        %f(a),f(c) have same signs so f(c),f(b) have different signs, zero is between a and c.
                        a = c; fa = fc;
                    end
                    c = (b+a)/2;
                end
                
                %add a new node at c
                k_full = FS_i + k - 1;
                
                %phi at c is done by interpolation, sliding segment is same
                %as in Kn, Kd evaluations
                sliding_seg = max(FS_i, min(FS_f - M_FS, k_full - floor(M_FS/2)));
                shape_offset = k_full - sliding_seg;
                
                phic = dot(((c+shape_offset).^(0:M_FS) * Lmat_FS'), phi_FS(k-shape_offset:k-shape_offset+M_FS));
                
                obj.boundary.x = [obj.boundary.x(1:k_full) dot(ax,c.^(3:-1:0)) obj.boundary.x(k_full+1:obj.boundary.N)];
                obj.boundary.z = [obj.boundary.z(1:k_full) dot(az,c.^(3:-1:0)) obj.boundary.z(k_full+1:obj.boundary.N)];
                obj.boundary.doublenode = [obj.boundary.doublenode(1:k_full) 0 obj.boundary.doublenode(k_full+1:obj.boundary.N)];
                obj.boundary.phi_FS = [obj.boundary.phi_FS(1:k) phic obj.boundary.phi_FS(k+1:length(obj.boundary.phi_FS))];
                
                obj.boundary.FS_end = obj.boundary.FS_end + 1;
                obj.boundary.N = obj.boundary.N + 1;
                
            end
        end
        
        function obj = handle_regrid(obj,force)
            % Checks if regridding should be done (or forces it if force=1)
            %
            % Regridding and checks are handled by the parameters in
            % the regridding struct.
            %   regrid_type:
            %     'none' - no regriding is done. force_regrid does not do
            %           anything.
            %     'shift_by_curve3d' - uses curve_renode_3d() to evenly
            %           distribute by curve length in (x,z,phi * a/Umax)
            %           coordinates, where a is specified by
            %           regrid_param{1}, and Umax is calculated as
            %           the maximum sqrt(phi_s^2 + phi_n^2) on the free
            %           surface. The same number of nodes are used, and
            %           regrid_param{2} specifies the interpolation scheme,
            %           which is passed as the interp_scheme argument into
            %           the curve_3d functions.
            %     'shift_by_curve2d_normu' -uses curve_renode_by_integral()
            %           using the integrand sqrt(1 + (|u|*a/umax)^2), which
            %           is similar to shift_by_curve3d, except
            %           |nabla(phi)|^2 is used instead of |phi_s|^2. Like
            %           shift_by_curve3d, regrid_param{1} specifies a.
            %   threshold_type:
            %     'manual' - regridding is never performed unless force is
            %           set to 1.
            %     'factor_global' - regridding is performed when the
            %           metric differs from the desired value of the metric
            %           by a factor of threshold_param{1}. for
            %           shift_by_curve3d, this metric is the segment
            %           length/total length in 3d. The factor must be
            %           greater than 1.
            %     'factor_local' - regridding is performed when the
            %           metric differs from the metric of adjacent segments
            %           by a factor of threshold_param{1}. for
            %           shift_by_curve3d, this metric is the segment
            %           length/total length in 3d. The factor must be
            %           greater than 1.
            %     'factor_global_local' - uses both modes 'factor_desired'
            %           and 'factor_neighbors' with factors
            %           threshold_param{1} and threshold_param{2},
            %           respectively.
            valid_threshold_types = {'manual','factor_global','factor_local','factor_global_local'};
            valid_regrid_types = {'none','shift_by_curve3d','shift_by_curve2d_normu'};
            %check valid types
            if ~any(strcmp(valid_regrid_types,obj.regridding.regrid_type))
                fprintf('regrid_type "%s" is invalid. Change it in the regridding struct!\n',obj.regridding.regrid_type);
                return;
            end
            if ~any(strcmp(valid_threshold_types,obj.regridding.threshold_type))
                fprintf('threshold_type "%s" is invalid. Change it in the regridding struct!\n',obj.regridding.threshold_type);
                return;
            end
            obj.regridding.did_regrid = 0;
            if strcmp(obj.regridding.regrid_type,'none')
                return;
            end
            %exit conditions: various threshold types fail
            if ~force
                if strcmp(obj.regridding.regrid_type,'manual')
                end
                if strcmp(obj.regridding.regrid_type,'shift_by_curve3d')
                    obj = obj.update_characteristics();
                    FS_i = obj.boundary.FS_start;
                    FS_f = obj.boundary.FS_end;
                    a = obj.regridding.regrid_param{1};
                    V = [obj.boundary.x(FS_i:FS_f); obj.boundary.z(FS_i:FS_f);obj.boundary.phi_FS *a/obj.characteristics_.umax];
                    seglens = curve_lengths_3d(V,obj.regridding.regrid_param{2});
                    seglens = (seglens / sum(seglens)) * (FS_f-FS_i); %divide by average length
                    if strcmp(obj.regridding.threshold_type,'factor_global')
                        factor = obj.regridding.threshold_param{1};
                        if max(seglens) < factor && min(seglens) > 1/factor
                            %within expected range, so no regrid
                            return;
                        end
                    elseif strcmp(obj.regridding.threshold_type,'factor_local')
                        neighbor_ratios = seglens(1:end-1) ./ seglens(2:end);
                        factor = obj.regridding.threshold_param{1};
                        if max(neighbor_ratios) < factor && min(neighbor_ratios) > 1/factor
                            %within expected range, so no regrid
                            return;
                        end
                    elseif strcmp(obj.regridding.threshold_type,'factor_global_local')
                        neighbor_ratios = seglens(1:end-1) ./ seglens(2:end);
                        factor_global = obj.regridding.threshold_param{1};
                        factor_local = obj.regridding.threshold_param{2};
                        if max(neighbor_ratios) < factor_local && min(neighbor_ratios) > 1/factor_local...
                                && max(seglens) < factor_global && min(seglens) > 1/factor_global
                            %within expected range, so no regrid
                            return;
                        end
                    end
                elseif strcmp(obj.regridding.regrid_type,'shift_by_curve2d_normu')
                    obj = obj.update_characteristics();
                    FS_i = obj.boundary.FS_start;
                    FS_f = obj.boundary.FS_end;
                    a = obj.regridding.regrid_param{1};
                    V = [obj.boundary.x(FS_i:FS_f); obj.boundary.z(FS_i:FS_f); obj.boundary.phi_FS];
                    vnorm2 = obj.characteristics_.vnorm2_FS;
                    seglens = curve_lengths_by_integral(V(1:2,:),sqrt(1 + vnorm2*(a/obj.characteristics_.umax)^2));
                    seglens = (seglens / sum(seglens)) * (FS_f-FS_i); %divide by average length
                    if strcmp(obj.regridding.threshold_type,'factor_global')
                        factor = obj.regridding.threshold_param{1};
                        if max(seglens) < factor && min(seglens) > 1/factor
                            %within expected range, so no regrid
                            return;
                        end
                    elseif strcmp(obj.regridding.threshold_type,'factor_local')
                        neighbor_ratios = seglens(1:end-1) ./ seglens(2:end);
                        factor = obj.regridding.threshold_param{1};
                        if max(neighbor_ratios) < factor && min(neighbor_ratios) > 1/factor
                            %within expected range, so no regrid
                            return;
                        end
                    elseif strcmp(obj.regridding.threshold_type,'factor_global_local')
                        neighbor_ratios = seglens(1:end-1) ./ seglens(2:end);
                        factor_global = obj.regridding.threshold_param{1};
                        factor_local = obj.regridding.threshold_param{2};
                        fprintf('local [%g,%g]; global [%g,%g]',min(neighbor_ratios),max(neighbor_ratios),min(seglens),max(seglens));
                        if max(neighbor_ratios) < factor_local && min(neighbor_ratios) > 1/factor_local...
                                && max(seglens) < factor_global && min(seglens) > 1/factor_global
                            %within expected range, so no regrid
                            return;
                        end
                    end
                end
            end

            %do the regrid
            if strcmp(obj.regridding.regrid_type,'shift_by_curve3d')
                if force || ~(...
                        strcmp(obj.regridding.threshold_type,'factor_local') ||...
                        strcmp(obj.regridding.threshold_type,'factor_global') ||...
                        strcmp(obj.regridding.threshold_type,'factor_global_local'))
                    %have not calculated these already, so do that
                    obj = obj.update_characteristics();
                    FS_i = obj.boundary.FS_start;
                    FS_f = obj.boundary.FS_end;
                    a = obj.regridding.regrid_param{1};
                    V = [obj.boundary.x(FS_i:FS_f); obj.boundary.z(FS_i:FS_f);obj.boundary.phi_FS *(a/obj.characteristics_.umax)];
                end
                Vnew = curve_renode_3d(V,FS_f-FS_i+1,obj.regridding.regrid_param{2});
                obj.boundary.x(FS_i:FS_f) = Vnew(1,:);
                obj.boundary.z(FS_i:FS_f) = Vnew(2,:);
                obj.boundary.phi_FS = Vnew(3,:) * (obj.characteristics_.umax/a);
            elseif strcmp(obj.regridding.regrid_type,'shift_by_curve2d_normu')
                if force || ~(...
                        strcmp(obj.regridding.threshold_type,'factor_local') ||...
                        strcmp(obj.regridding.threshold_type,'factor_global') ||...
                        strcmp(obj.regridding.threshold_type,'factor_global_local'))
                    %have not calculated these already, so do that
                    obj = obj.update_characteristics();
                    FS_i = obj.boundary.FS_start;
                    FS_f = obj.boundary.FS_end;
                    a = obj.regridding.regrid_param{1};
                    V = [obj.boundary.x(FS_i:FS_f); obj.boundary.z(FS_i:FS_f); obj.boundary.phi_FS];
                    vnorm2 = obj.characteristics_.vnorm2_FS;
                end
                Vnew = curve_renode_by_integral(V,sqrt(1 + vnorm2*(a/obj.characteristics_.umax)^2),FS_f-FS_i+1,[1;1;0]);
                obj.boundary.x(FS_i:FS_f) = Vnew(1,:);
                obj.boundary.z(FS_i:FS_f) = Vnew(2,:);
                obj.boundary.phi_FS = Vnew(3,:);
            end
            %force a reload of the characteristics
            obj.regridding.did_regrid = 1;
            obj.characteristics_.ctime = 'null';
        end
    end
end

