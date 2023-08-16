classdef vertical_wavemaker_bdry < boundary
    %VERTICAL_WAVEMAKER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cond_type = "neumann"
        update_orders = [-1]
        step_orders = [0.1]
        %u(t,x,z) - function of horizontal velocity of inflow wave
        %should expect z as a 1xn array.
        func_u
        func_D %D(t) - ramp function (0 -> 1)

        %direction wave travels in (1 for normal direction,
        % -1 for reverse normal)
        wavedir
        has_function_handles = 1
    end
    
    methods
        function obj = vertical_wavemaker_bdry(varargin)
            %SOLID_BDRY Construct an instance of this class
            %   Detailed explanation goes here
            obj.init_surface(varargin{:});
            obj.wavedir = NaN;
        end
        
        function on_update(obj,priority)
            if isnan(obj.wavedir)
                %make a guess on the wave direction based on FS side
                CCW = obj.intersection_CCW.CCW_link;
                CW = obj.intersect_CW_link;
                if isa(CCW,"free_surface_bdry") && ...
                        ~isa(CW,"free_surface_bdry")
                    obj.wavedir = 1;
                elseif ~isa(CCW,"free_surface_bdry") && ...
                        isa(CW,"free_surface_bdry")
                    obj.wavedir = -1;
                else
                    error("Vertical wavemaker boundaries must have " + ...
                        "exactly one adjacent free surface boundary");
                end
            end
            %x should all be the same. Just pull first
            x = obj.boundary_nodes(1,1);

            %set phi_n and phi_tn according to u and D
            EPS = 1e-4;
            t = obj.parent_sim.stepping.t;
            z = obj.boundary_nodes(:,2)';
            D_t = (obj.func_D(t+EPS) - obj.func_D(t-EPS))/(2*EPS);
            u_t = (obj.func_u(t+EPS,x,z) - obj.func_u(t-EPS,x,z))/(2*EPS);
            
            obj.characteristics.phi_n = zeros(1,obj.node_count);
            obj.characteristics.phi_tn = zeros(1,obj.node_count);

            obj.characteristics.phi_n(:) = ...
                (obj.wavedir * obj.func_D(t)) .* obj.func_u(t,x,z);
            obj.characteristics.phi_tn(:) = ...
                (obj.wavedir * D_t) .* obj.func_u(t,x,z)...
                + (obj.wavedir * obj.func_D(t)) .* u_t;
        end

        function on_step(obj, dt, order)
            %lock x position to match free surface
            CCW = obj.intersection_CCW.CCW_link;
            CW = obj.intersect_CW_link;
            if isa(CCW,"free_surface_bdry") && ...
                    ~isa(CW,"free_surface_bdry")
                %CCW side is FS; match that
                obj.boundary_nodes(:,1) = CCW.boundary_nodes(1,1);
                %expect to be on the right side of the domain
                if isnan(obj.wavedir)
                    obj.wavedir = 1;
                end
            elseif ~isa(CCW,"free_surface_bdry") && ...
                    isa(CW,"free_surface_bdry")
                %CW side is FS; match that
                obj.boundary_nodes(:,1) = CW.boundary_nodes(end,1);

                %expect to be on the left side of the domain
                if isnan(obj.wavedir)
                    obj.wavedir = -1;
                end
            else
                error("Vertical wavemaker boundaries must have " + ...
                    "exactly one adjacent free surface boundary");
            end
        end


        function default_gen_functions(obj,h,lambda,amplitude,ramp_start...
                ,ramp_end)
            g = obj.parent_sim.meta.g;
            k = (2*pi)/lambda;
            omega = sqrt( g*k*tanh(k*h) ); %% FF added full dispersion relationship
            obj.func_u = @(t,x,z)...
                amplitude.*omega.* (cosh(k.*(z+h))./sinh(k*h)).*cos(k.*x - omega.*t); %amplitude.*omega./(k*h).*cos(k.*x - omega.*t);
            obj.func_D = @(t) 1 .* (t >= ramp_end)...
                + (1-cos((t-ramp_start)*pi/(ramp_end-ramp_start)))./2 ...
                .* (t > ramp_start & t < ramp_end);
            disp(sprintf('WM: default_gen_func: depth=%2.f (m), a0=%.3f (m), lambda=%.2f (m), omega=%.2f',h,amplitude,lambda,omega))
        end

        function ff_gen_functions(obj,h,freq,amplitude,phase,ramp_end)
            g = obj.parent_sim.meta.g;
            omega = 2*pi*freq;
            k = get_wavenum(omega, h);
            kh = k.*h
                        %w = sqrt(g*h)*k;
            obj.func_u = @(t,x,z)...
                sum( amplitude.*omega.* (cosh(k.*(z+h))./sinh(kh)).*cos(k.*x - omega.*t + phase));
            obj.func_D = @(t) 1 .* (t >= ramp_end)...
                + (1-cos((t)*pi/(ramp_end)))./2 ...
            .* (t >= 0 & t < ramp_end);
        end

        function packet_gen_functions(obj,h,freq,amplitude,packet_width,ramp_end)
            g = obj.parent_sim.meta.g;
            omega = 2*pi*freq;
            k = get_wavenum(omega, h);
            kh = k.*h
                        %w = sqrt(g*h)*k;
            obj.func_u = @(t,x,z)...
                exp( -((t-2*packet_width)/packet_width)^2).*amplitude.*omega.*  (cosh(k.*(z+h))./sinh(kh)).*cos(k.*x - omega.*t );
            obj.func_D = @(t) 1 .* (t >= ramp_end)...
                + (1-cos((t)*pi/(ramp_end)))./2 ...
            .* (t >= 0 & t < ramp_end);
        end
        

    end
end

