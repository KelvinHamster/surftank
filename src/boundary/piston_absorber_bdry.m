classdef piston_absorber_bdry < boundary
    %piston_absorber_bdry Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cond_type = "neumann"

        update_orders = [-0.5,0.5]
        step_orders = [0.1]

        %interp to calculate phi_ss
        char_interp

        %stationary water depth
        h

        %data/estimates used for updates
        u_upd
        du_upd
        phit_upd

        u_prev
        phit_prev
        z_prev
        t_prev

        normdir %direction of normal (+1 for +x, -1 for -x)

        has_function_handles = 0
    end
    
    methods
        function obj = piston_absorber_bdry(varargin)
            %Construct an instance of this class
            %   Detailed explanation goes here

            obj.char_interp = get_interpolation_struct(4,0,1);
            obj.init_surface(varargin{:});
            
            obj.u_upd = 0;
            obj.du_upd = 0;
            obj.phit_upd = [];

            obj.u_prev = 0;
            obj.phit_prev = [];
            obj.z_prev = [];
            obj.t_prev = NaN;

            obj.h = NaN;
        end
        
        function on_update(obj,order)

            if isnan(obj.h)
                obj.h = abs(obj.boundary_nodes(end,2)...
                          - obj.boundary_nodes(1,2));
                fprintf("stationary depth 'h' not set on absorber." + ...
                    " Defaulting to %f.\n",obj.h);
            end

            obj.normdir = sign(obj.boundary_nodes(end,2)...
                             - obj.boundary_nodes(1,2));
            if isnan(obj.t_prev) ...
                    || isempty(obj.z_prev)
                obj.u_upd = 0;
                obj.du_upd = 0;
                %take role of stationary solid if data is missing
                % only pre-BEM1 step needed
                if order == obj.update_orders(1)
                    obj.characteristics.phi_n = ...
                        zeros(1,obj.node_count);
                    obj.characteristics.phi_tn = ...
                        zeros(1,obj.node_count);
                end
            else
                if order == obj.update_orders(1)
                    obj.pre_BEM1();
                end
            end
        end
        function pre_BEM1(obj)

            % data is here, so treat as piston: need phi_n
            % which requires u, which needs phi_t. The troublesome part is
            % that estimating phi_t for this step is hard without
            % knowing phi for this step.

            %integrate for hydrodynamic force, divided by rho
            F = -0.5 * sum((obj.phit_upd(1:end-1) + obj.phit_upd(2:end))...
                .* diff(obj.z_prev)) * obj.normdir;
            u = (F/(obj.h * sqrt(obj.h * obj.parent_sim.meta.g)) ...
                + obj.u_prev) * 0.5;
            
            obj.du_upd = (u - obj.u_prev)/...
                (obj.parent_sim.stepping.t - obj.t_prev);
            obj.u_upd = u;
            
            obj.characteristics.phi_n = ...
                        zeros(1,obj.node_count) + u * obj.normdir;
        end
        function pre_BEM2(obj)
            
            n = obj.node_count;
            M = obj.char_interp.M;
            Lpmat = obj.char_interp.Lagrange_prime;
            Lppmat = obj.char_interp.Lagrange_primeprime;
            mnode = linspace(0,1,M+1);

            % calculate phi_tn, which requires Du/Dt, u, and phi_ss
            % (Du/Dt, u calculated on pre_BEM1)
            phi_ss = zeros(1,n);
            phi = obj.characteristics.phi;

            for j = 1:n
                %place the sliding element,
                %preferably so j is the midpoint
                k = j - floor(M/2);
                %shift so it lies exclusively on this boundary
                if k < 1
                    k = 1;
                elseif k+M > n
                    k = n-M;
                end
                
                t_interp = mnode(j - k + 1);
                seg_z = obj.boundary_nodes(k:k+M,2);
                
                zp = (t_interp.^(0:(M-1))*Lpmat') * seg_z;
                zpp = (t_interp.^(0:(M-2)) * Lppmat') * seg_z;
                jacobian = abs(zp);
                phi_s = dot((t_interp.^(0:(M-1))*Lpmat'),...
                    phi(k:k+M))./jacobian;


                phi_ss(j)=(dot((t_interp.^(0:(M-2))*Lppmat'),phi(k:k+M))...
                    - zp*zpp*phi_s)./(jacobian.^3);
                
            end
            obj.characteristics.phi_tn = ...
                        obj.du_upd * obj.normdir...
                        - (obj.u_upd * obj.normdir) * phi_ss;

        end

        function on_step(obj,dt,order)
            %case: missing data:
            if isnan(obj.t_prev) ...
                    || isempty(obj.z_prev)
                %solid bdry, but initialize some stuff:

                obj.phit_upd = zeros(1,obj.node_count);
    
                obj.u_prev = obj.u_upd;
                obj.phit_prev = obj.characteristics.phi_t;
                obj.z_prev = obj.boundary_nodes(:,2)';
                obj.t_prev = obj.parent_sim.stepping.t;
                return
            end

            %update x pos
            x = obj.boundary_nodes(1,1) + obj.u_upd * dt...
                + 0.5 * obj.du_upd * dt^2;
            obj.u_prev = obj.u_upd;

            CCW = obj.intersection_CCW.CCW_link;
            CW = obj.intersect_CW_link;
            if isa(CCW,"free_surface_bdry") && ...
                    ~isa(CW,"free_surface_bdry")
                %CCW side is FS; match that
                obj.boundary_nodes(:,1) = x;
            elseif ~isa(CCW,"free_surface_bdry") && ...
                    isa(CW,"free_surface_bdry")
                %CW side is FS; match that
                obj.boundary_nodes(:,1) = x;
            else
                error("Vertical wavemaker boundaries must have " + ...
                    "exactly one adjacent free surface boundary");
            end

            %in case z changed, we need to shift phi_prev and phit_upd
            %to these values of z. Utilize sliding segments.
            n = length(obj.z_prev);
            M = obj.char_interp.M;
            Lmat = obj.char_interp.Lagrange;
            mnode = linspace(0,1,M+1);

            phit_prev_ = zeros(1,obj.node_count);

            for j = 1:(n-1)
                %place the sliding element,
                %preferably so j is the midpoint
                k = j - floor(M/2);
                %shift so it lies exclusively on this boundary
                if k < 1
                    k = 1;
                elseif k+M > n
                    k = n-M;
                end
                
                %parameterize linearly
                param = (obj.boundary_nodes(:,2) - obj.z_prev(j))/...
                    (obj.z_prev(j+1) - obj.z_prev(j));

                %find nodes inside segment
                within = (0 <= param) & (param < 1);
                %for first and last segments, include exteriors
                if j==1
                    within = within | (param < 0);
                elseif j == n-1
                    within = within | (param >= 1);
                end
                
                phit_prev_(within) = (param(within).^(0:M) * Lmat')...
                    * obj.phit_prev(k:k+M)';
            end



            
            t = obj.parent_sim.stepping.t;
            dt1 = t - obj.t_prev;
            %TODO calc phi_t for next pre_BEM1
            %     dt0   dt1 T  dt2
            %   |-----|-----|------|
            %      phi_t---phi_t-->|
            %     utilize phi_t's from BEM2 to get to T+dt2


            %phi_t(T+dt2) ~ phi_t(T) + (Dphi_t/Dt)*dt2
            % where 
            % Dphi_t/Dt ~ (phi_t(T) - phi_t(T-dt1))/dt1
            coef = dt/dt1;
            obj.phit_upd = obj.characteristics.phi_t * (1 + coef)...
                - phit_prev_ * coef;
            obj.z_prev = obj.boundary_nodes(:,2)';

            %need this for next step
            obj.phit_prev = obj.characteristics.phi_t;
            obj.t_prev = t;
        end
    end
end

