classdef bem_sim < handle
    %Manages a simulator
    %   Detailed explanation goes here
    
    properties
        boundaries
        meta
        stepping
        global_step_events
        global_update_events
        global_step_event_orders
        global_update_event_orders
        regridding
    end
    
    methods
        function obj = bem_sim(boundaries,varargin)
            obj.global_step_events = {};
            obj.global_update_events = {};
            obj.global_step_event_orders = [];
            obj.global_update_event_orders = [];
            
            obj.boundaries = boundaries;
            %link boundaries
            for i=2:length(boundaries)
                intersection_build_default(boundaries{i-1},boundaries{i});
                obj.boundaries{i}.parent_sim = obj;
            end
            obj.boundaries{1}.parent_sim = obj;
            intersection_build_default(boundaries{end},boundaries{1});

            % args
            parallel_workers = 0;

            p = inputParser;
            addOptional(p,'parallel_workers',parallel_workers, @(x) ...
               isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1) == 0) );

            parse(p, varargin{:});

            parallel_workers = p.Results.parallel_workers;
            dt = 0.1;
            store_dt = 0.1;
            timestep_method = 'taylor2';
            courant_lock = 1;
            courant_number = 0.5;
            g = 9.81;
            Uwave = 0;

            %meta
            obj.meta = struct('parallel_workers',parallel_workers,...
                'char_t','null','g',g);
            
            
            obj.stepping = struct('t',0, 'dt',dt, 'store_dt', store_dt,...
                                  'Uwave',Uwave,...
                'timestep_method', timestep_method, ...
                'courant_lock',courant_lock,...
                'courant_number',courant_number, 'dtmax', 0.1);
        end

        function [nodes,starts] = get_nodelist(obj)
            % Generates a Nx2 array of nodes, across all of the boundaries.
            % Ignores the last node of each boundary in order to prevent
            % repetitions on double-nodes.
            %
            % @returns [nodes,starts] where nodes is the Nx2 array, and
            % starts is an array that gives the index in 'nodes'
            % corresponding to the first node in each boundary, that is,
            % nodes(starts(i),:) == boundaries{i}.boundary_nodes(1,:)

            N = 0;
            starts = zeros(1,length(obj.boundaries));
            for i= 1:length(obj.boundaries)
                starts(i) = N+1;
                N = N + obj.boundaries{i}.node_count - 1;
            end
            
            nodes = zeros(N,2);
            k = 1;
            
            for i= 1:length(obj.boundaries)
                n = obj.boundaries{i}.node_count;
                %skip last node of bdry, since it will be added
                %in the next one. (double nodes are only observed
                %from the CW side of the CCW boundary for our purposes.)
                nodes(k:k+n-2,:)=obj.boundaries{i}.boundary_nodes(1:n-1,:);
            
                k = k + n - 1;
            end
        end
        
        function update_characteristics(obj)
            %Checks to ensure that the BEM has been solved for this step
            if (ischar(obj.meta.char_t) &&...
                    strcmp(obj.meta.char_t,'null'))...
                    ||(obj.meta.char_t ~= obj.stepping.t)
                NUM_BDRYS = length(obj.boundaries);
                %=======================pre_solve
                %get list of priorities
                prios = [0,1,obj.global_update_event_orders];
                for i=1:NUM_BDRYS
                    bdry = obj.boundaries{i};
                    prios = [prios,...
                        bdry.update_orders];

                    if isempty(bdry.boundary_nodes)
                        bdry.reassemble_nodes();
                    end
                end
                %go up the priority list
                while ~isempty(prios)
                    prio = min(prios);
                    prios = prios(prios > prio);

                    if prio == 0
                        %solve system at order 0
                        BIE_integ_calc(obj);
                        [A,b]= matrix_build(obj,'build_A',1,'build_b',1,...
                            'b_fieldsrc_u','phi','b_fieldsrc_q','phi_n');
                        x = A\b;
                        j = 1;
                        for i=1:NUM_BDRYS
                            bdry = obj.boundaries{i};
                            n = bdry.node_count;
                            if strcmp(bdry.cond_type,"dirichlet")
                                bdry.characteristics.phi_n = x(j:j+n-1)';
                            elseif strcmp(bdry.cond_type,"neumann")
                                bdry.characteristics.phi = x(j:j+n-1)';
                            end
                            j = j + n;
                        end
                    elseif prio == 1
                        %solve system at order 1
                        [A,b]= matrix_build(obj,'build_A',1,'build_b',1,...
                           'b_fieldsrc_u','phi_t','b_fieldsrc_q','phi_tn');
                        x = A\b;
                        j = 1;
                        for i=1:NUM_BDRYS
                            bdry = obj.boundaries{i};
                            n = bdry.node_count;
                            if strcmp(bdry.cond_type,"dirichlet")
                                bdry.characteristics.phi_tn = x(j:j+n-1)';
                            elseif strcmp(bdry.cond_type,"neumann")
                                bdry.characteristics.phi_t = x(j:j+n-1)';
                            end
                            j = j + n;
                        end
                    end

                    for i=1:length(obj.global_update_event_orders)
                        if prio == obj.global_update_event_orders(i)
                            obj.global_update_events{i}(obj);
                        end
                    end
                    
                    for i=1:NUM_BDRYS
                        bdry = obj.boundaries{i};
                        if ismember(prio,bdry.update_orders)
                            bdry.on_update(prio);
                        end
                    end
                end
                %mark this time so that we don't have to call this again
                obj.meta.char_t = obj.stepping.t;
            end
        end

        function obj = step(obj,dt_overwrite)
            %Steps the simulation forward by dt, updating the free surface.
            % if dt_overwrite is specified, then steps by that; otherwise
            % defaults to bem_sim settings (bem_sim.stepping)
            % note that dt_overwrite will affect dt and courant_number in
            % the stepping struct.


            NUM_BDRYS = length(obj.boundaries);
            %=======================calculate the time step size
            dxmin = inf;
            umax = 0;
            for i=1:NUM_BDRYS
                bdry = obj.boundaries{i};
                if isfield(bdry.meta,'cfl_contribute') &&...
                        bdry.meta.cfl_contribute
                    dxmin = min(dxmin, bdry.characteristics.dxmin);
                    umax = max(umax, bdry.characteristics.umax);
                end
            end
            if umax<obj.stepping.Uwave
                umax = obj.stepping.Uwave;
            end
            if exist("dt_overwrite","var")
                %overwritten from step(dt) call.
                dt = dt_overwrite;
            elseif obj.stepping.courant_lock == 1
                %set dt to courant number
                dt = obj.stepping.courant_target * dxmin / umax;
            elseif obj.stepping.courant_lock == 0
                %fixed dt, just update courant number
                dt = obj.stepping.dt;
            elseif isfield(obj.stepping,"get_step_size") &&...
                    isa(obj.stepping.get_step_size,"function_handle")
                %manual stepsize calculation
                dt = obj.stepping.get_step_size(obj);

            else
                error("Step size choice unknown. Set courant_lock" + ...
                    "to 0 or 1, or set function handle get_step_size!")
            end
            
            %cap dt
            if isfield(obj.stepping,"dtmax")
                dt = min(dt,obj.stepping.dtmax);
            end
            if isfield(obj.stepping,"calc_courant")
                obj.stepping.courant_number =...
                    obj.stepping.calc_courant(obj);
            else
                obj.stepping.courant_number = umax * dt / dxmin;
            end
            obj.stepping.dt = dt;


            if isfield(obj.stepping,"dtmin") && ...
                    dt < obj.stepping.dtmin
                error("Step size (%.4e) dropped below specified " + ...
                    "'dtmin' value of (%.4e)",dt,obj.stepping.dtmin)
            end



            %=======================do the step
            %get list of priorities
            prios = obj.global_step_event_orders;
            for i=1:NUM_BDRYS
                prios = [prios,...
                    obj.boundaries{i}.step_orders];
            end
            %go up the priority list
            while ~isempty(prios)
                prio = min(prios);
                prios = prios(prios > prio);

                for i=1:length(obj.global_step_event_orders)
                    if prio == obj.global_step_event_orders(i)
                        obj.global_step_events{i}(obj);
                    end
                end
                
                for i=1:NUM_BDRYS
                    bdry = obj.boundaries{i};
                    if ismember(prio,bdry.step_orders)
                        bdry.on_step(dt,prio);
                    end
                end
            end
            obj.stepping.t = obj.stepping.t + dt;

            %=======================resolve bdry intersections
            obj.enforce_intersections();
        end

        function enforce_intersections(obj)
            for i=1:length(obj.boundaries)
                intersection_enforce(obj.boundaries{i});
            end
        end
        
        function plot_full(obj,omega)
            %Plots the full boundary in figure 1.
            % If omega is given as an optional argument, the title will
            % show the fraction of a period the simulation is currently at.
            %
            
            t = obj.stepping.t;
            
            %plot free surface boundary and nodes
            %set(0,'CurrentFigure',1);
            if ~ishandle(1)
                figure(1);
            else
                set(0,'CurrentFigure',1);
            end
            for i=1:length(obj.boundaries)
                bdry = obj.boundaries{i};
                plot(bdry.boundary_nodes(:,1),bdry.boundary_nodes(:,2))
                hold on
                scatter(bdry.boundary_nodes(:,1),...
                    bdry.boundary_nodes(:,2),'o')
            end
            
            if isfield(obj.meta,'plot_ylim')
                ylim(obj.meta.plot_ylim)
            end
            ylabel('z')
            if isfield(obj.meta,'plot_xlim')
                xlim(obj.meta.plot_xlim)
            end
            xlabel('x')
            
            if (~exist('omega','var')) && isfield(obj.meta,'omega')
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

        function regrid(obj)
            is_updated = (ischar(obj.meta.char_t) &&...
                    strcmp(obj.meta.char_t,'null'))...
                    ||(obj.meta.char_t ~= obj.stepping.t);
            %get list of priorities
            prios = [];
            NUM_BDRYS = length(obj.boundaries);
            for i=1:NUM_BDRYS
                bdry = obj.boundaries{i};
                if (~strcmp(bdry.regridding.mode,'none')) &&...
                    (bdry.regridding.precede_update || is_updated)
                    prios = [prios,...
                        bdry.regridding.order];
                end
            end
            %go up the priority list
            while ~isempty(prios)
                prio = min(prios);
                prios = prios(prios > prio);
                
                %current priority, make sure that there is actually desired
                %regridding. Additionally, if we do not have updated
                %characteristics, do the non-precede_update parts in a
                %second pass.
                for i=1:NUM_BDRYS
                    bdry = obj.boundaries{i};
                    if prio == bdry.regridding.order &&...
                        (~strcmp(bdry.regridding.mode,'none')) &&...
                        (bdry.regridding.precede_update || is_updated)

                        if bdry_handle_regrid(bdry)
                            obj.meta.char_t = 'null';
                        end
                    end
                end
            end
            
            if ~is_updated
                obj.update_characteristics();
                %get list of priorities
                prios = [];
                for i=1:NUM_BDRYS
                    bdry = obj.boundaries{i};
                    if (~strcmp(bdry.regridding.mode,'none')) &&...
                        ~(bdry.regridding.precede_update || is_updated)
                        prios = [prios,...
                            bdry.regridding.order];
                    end
                end
                %go up the priority list
                while ~isempty(prios)
                    prio = min(prios);
                    prios = prios(prios > prio);
                    
                    %second pass, same logic
                    for i=1:NUM_BDRYS
                        bdry = obj.boundaries{i};
                        if prio == bdry.regridding.order &&...
                            (~strcmp(bdry.regridding.mode,'none')) &&...
                            ~(bdry.regridding.precede_update || is_updated)
    
                            if bdry_handle_regrid(bdry)
                                obj.meta.char_t = 'null';
                            end
                        end
                    end
                end
            end
        end
        
        function full_step(obj)
            obj.update_characteristics();
            obj.step();
            obj.regrid();
            obj.update_characteristics();
        end

        function append_update_event(obj,handle,order)
            n = length(obj.global_update_event_orders);
            if n ~= length(obj.global_update_events)
                error("Global update event arrays do not match!" + ...
                    " The list of orders do not line up with the" + ...
                    " list of events!")
            end
            obj.global_update_events{n+1} = handle;
            obj.global_update_event_orders(n+1) = order;
        end
        function append_step_event(obj,handle,order)
            n = length(obj.global_step_event_orders);
            if n ~= length(obj.global_step_events)
                error("Global update event arrays do not match!" + ...
                    " The list of orders do not line up with the" + ...
                    " list of events!")
            end
            obj.global_step_events{n+1} = handle;
            obj.global_step_event_orders(n+1) = order;
        end

        function phi = interior_eval(obj,x,z,solve_phi_t)
            %Evaluates phi on an interior point (x,z)
            obj.update_characteristics();

            if ~exist("solve_phi_t","var")
                solve_phi_t = 0;
            end

            phi = cell(1,length(obj.boundaries));
            phin = cell(1,length(obj.boundaries));
            
            if solve_phi_t
                for i=1:length(obj.boundaries)
                    phi{i} = obj.boundaries{i}.characteristics.phi_t;
                    phin{i} = obj.boundaries{i}.characteristics.phi_tn;
                end
            else
                for i=1:length(obj.boundaries)
                    phi{i} = obj.boundaries{i}.characteristics.phi;
                    phin{i} = obj.boundaries{i}.characteristics.phi_n;
                end
            end
            
            phi = bem_integrate(obj,phin,3,[x,z])...
                - bem_integrate(obj,phi,4,[x,z]);
        end
        
        function result = is_point_inside(obj,xp,zp)
            % returns whether or not (xp,zp)
            % is on the interior of the boundary

            nodes = obj.get_nodelist();
            x = nodes(:,1); z = nodes(:,2);
            N = length(x);
            
            %even odd rule
            %(https://en.wikipedia.org/wiki/Even-odd_rule)
            j = N;
            result = false;
            
            for i=1:N
                
                if (z(i) > zp) ~= (z(j) > zp)
                    %crossed z = zp line
                    dotnorm=(xp-x(i))*(z(j)-z(i)) - (zp-z(i))*(x(j)-x(i));
                    
                    
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
        
    end
end

