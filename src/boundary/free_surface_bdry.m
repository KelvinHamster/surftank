classdef free_surface_bdry < boundary
    %FREE_SURFACE_BDRY Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties
        cond_type = "dirichlet"

        update_orders = [-1,0.5,1.5]
        step_orders = [0]
        char_interp %characteristic interpolation
        
        has_function_handles = 0
    end
    
    methods
        function obj = free_surface_bdry(varargin)
            %FREE_SURFACE_BDRY Construct an instance of this class
            %   Detailed explanation goes here
            obj.init_surface(varargin{:});

            obj.char_interp = get_interpolation_struct(4,0,1);
            obj.meta.cfl_contribute = 1;
        end
        
        function ensure_characteristics(obj)
            if ~isfield(obj.characteristics,'phi')
                fprintf("phi was not set on this free surface " + ...
                    "boundary. Defaulting to all zeros.\n");
                obj.characteristics.phi = zeros(1,obj.node_count);
            elseif obj.node_count ~= length(obj.characteristics.phi)
                error("The size of the phi array in this free " + ...
                    "surface boundary does not match the node count.\n");
            end
        end

        function on_update(obj, order)
            if order == obj.update_orders(1)
                obj.ensure_characteristics();
            elseif order == obj.update_orders(2)
                obj.post_BEM1();
            elseif order == obj.update_orders(3)
                obj.post_BEM2();
            end
        end

        function post_BEM1(obj)
            g = obj.parent_sim.meta.g;
            

            M = obj.char_interp.M;
            Lpmat = obj.char_interp.Lagrange_prime;
            Lppmat = obj.char_interp.Lagrange_primeprime;
            mnode = linspace(0,1,M+1);

            FS_size = obj.node_count;
            
            phi = obj.characteristics.phi;
            phi_n = obj.characteristics.phi_n;
            phi_ns = zeros(1,FS_size);
            phi_s = zeros(1,FS_size);
            phi_ss = zeros(1,FS_size);
            phi_t = zeros(1,FS_size);
            cosbeta = zeros(1,FS_size);
            sinbeta = zeros(1,FS_size);
            beta_s = zeros(1,FS_size);
            dxdt = zeros(1,FS_size);
            dzdt = zeros(1,FS_size);
            dphidt = zeros(1,FS_size);
            
            vnorm2_max = 0;
            for j = 1:FS_size
                %place the sliding element,
                %preferably so j is the midpoint
                k = j - floor(M/2);
                %shift so it lies exclusively on the free surface
                if k < 1
                    k = 1;
                elseif k+M > FS_size
                    k = FS_size-M;
                end
                
                t_interp = mnode(j - k + 1);
                seg_nodes = obj.boundary_nodes(k:k+M,:);
                
                %r'
                rp = (t_interp.^(0:(M-1))*Lpmat') * seg_nodes;
                %|r'|
                jacobian = norm(rp);

                %cos(beta) and sin(beta)
                cosbeta(j) = rp(1) / jacobian;
                sinbeta(j) = rp(2) / jacobian;

                
                %differentiate the interpolation of
                % phi on s to get phi_s
                phi_s(j) = dot((t_interp.^(0:(M-1))*Lpmat'),...
                    phi(k:k+M))./jacobian;
                %same thing with phi_n to phi_ns
                phi_ns(j) = dot((t_interp.^(0:(M-1))*Lpmat'),...
                    phi_n(k:k+M))./jacobian;
                %twice for phi_ss
                %      ((phi(r))'' = phi_s*|r'|' + phi_ss*|r'|)
                rpp = (t_interp.^(0:(M-2)) * Lppmat') * seg_nodes;
                %|r'|' = (r' dot r'')/jac
                shapeproj = dot(rp,rpp)./jacobian;
                phi_ss(j) = (dot((t_interp.^(0:(M-2))*Lppmat'),...
                    phi(k:k+M))   - shapeproj*phi_s(j))./(jacobian.^2);
                
               
                %beta_s: calculated using chain rule:
                %  beta(t) = atan(z'(t)/x'(t))
                % beta'(t) = atan'(z'(t)/x'(t)) * (z'(t)/x'(t))'
                %    = 1/(1 + (z'(t)/x'(t))^2)  * (x'(t)z''(t)
                %                     - z'(t)x''(t))/(x'(t))^2
                %using the substitution
                %    cosbeta = x'(t)/jacobian, sinbeta = z'(t)/jacobian
                %and converting from t space to s-space
                %beta_s = (cosbeta*z_ss - sinbeta * x_ss)/(jacobian^2);
                %derive it first, though; the signs may be wrong
                beta_s(j) = (rp(1)*rpp(2)-rp(2)*rpp(1))./(jacobian.^3);
                
                %use eulerian condition:
                % phi_t   = -1/2 |del phi|^2  - gz - P/rho   (P = 0)
                % Dphi/Dt = +1/2 |del phi|^2  - gz - P/rho   (P = 0)
                v_norm2 = phi_s(j)^2 + phi_n(j)^2;
                z = obj.boundary_nodes(j,2);
                phi_t(j)  = -v_norm2/2 - g * z;
                dphidt(j) =  v_norm2/2 - g * z;
                
                %space derivatives simply gradient of phi
                dxdt(j) = phi_s(j)*cosbeta(j) + phi_n(j)*sinbeta(j);
                dzdt(j) = phi_s(j)*sinbeta(j) - phi_n(j)*cosbeta(j);

                vnorm2_max = max(vnorm2_max,v_norm2);
            end
                
                
            obj.characteristics.phi_s = phi_s;
            obj.characteristics.phi_ns = phi_ns;
            obj.characteristics.phi_ss = phi_ss;
            obj.characteristics.phi_t = phi_t;
            obj.characteristics.cosbeta = cosbeta;
            obj.characteristics.sinbeta = sinbeta;
            obj.characteristics.beta_s = beta_s;
            obj.characteristics.dxdt = dxdt;
            obj.characteristics.dzdt = dzdt;
            obj.characteristics.dphidt = dphidt;
            obj.characteristics.dxmin = min(vecnorm( ...
                  obj.boundary_nodes(2:end,:) ...
                - obj.boundary_nodes(1:(end-1),:)...
                ,2,2));
            obj.characteristics.umax = sqrt(vnorm2_max);
        end


        function post_BEM2(obj)
            g = obj.parent_sim.meta.g;

            FS_size = obj.node_count;
            phi_t = obj.characteristics.phi_t;
            phi_tn = obj.characteristics.phi_tn;
            phi_ns = obj.characteristics.phi_ns;
            phi_s = obj.characteristics.phi_s;
            phi_n = obj.characteristics.phi_n;
            sinbeta = obj.characteristics.sinbeta;
            cosbeta = obj.characteristics.cosbeta;
            beta_s = obj.characteristics.beta_s;
            phi_ss = obj.characteristics.phi_ss;

            phi_ts = zeros(1,FS_size);
            d2xdt2 = zeros(1,FS_size);
            d2zdt2 = zeros(1,FS_size);
            d2phidt2 = zeros(1,FS_size);

            M = obj.char_interp.M;
            Lpmat = obj.char_interp.Lagrange_prime;
            mnode = linspace(0,1,M+1);

            for j = 1:FS_size
                %place the sliding element,
                %preferably so j is the midpoint
                k = j - floor(M/2);
                %shift so it lies exclusively on the free surface
                if k < 1
                    k = 1;
                elseif k+M > FS_size
                    k = FS_size-M;
                end
                
                t_interp = mnode(j - k + 1);
                seg_nodes = obj.boundary_nodes(k:k+M,:);

                %r'
                rp = (t_interp.^(0:(M-1))*Lpmat') * seg_nodes;
                %|r'|
                jacobian = norm(rp);
                
                %differentiate the interpolation phi_t ---> phi_ts
                phi_ts(j) = dot((t_interp.^(0:(M-1))*Lpmat'),...
                    phi_t(k:k+M))./jacobian;
                
                %D^2phi/Dt^2 = phi_s*phi_st + phi_n*phi_nt
                %            + phi_s(phi_s*phi_ss + phi_n*phi_ns)
                %            - phi_n(phi_n*phi_ss + phi_s*phi_ns)
                %            - phi_n*beta_s(phi_s^2 + phi_n^2)
                %            - g(phi_s*sinbeta - phi_n*cosbeta)
                %            - (1/rho)DP/Dt
                d2phidt2(j) = phi_s(j)*phi_ts(j)+phi_n(j)*phi_tn(j)...
                   +phi_s(j)*(phi_s(j)*phi_ss(j)+phi_n(j)*phi_ns(j))...
                   -phi_n(j)*(phi_n(j)*phi_ss(j)-phi_s(j)*phi_ns(j))...
                   -phi_n(j)*beta_s(j)*(phi_s(j)^2 + phi_n(j)^2)...
                   - g*(phi_s(j)*sinbeta(j) - phi_n(j)*cosbeta(j));
                
                %D^2x/Dt^2 = cosbeta(phi_ts+phi_s*phi_ss+phi_n*phi_ns)
                %    + sinbeta(
                %              phi_tn + phi_s*phi_ns
                %            - phi_n*phi_ss - beta_s(phi_n^2 + phi_s^2)
                %             )
                d2xdt2(j) = cosbeta(j)*(phi_ts(j)...
                    + phi_s(j)*phi_ss(j)...
                    + phi_n(j)*phi_ns(j))...
                    +sinbeta(j)*(phi_tn(j)...
                    + phi_s(j)*phi_ns(j)...
                    - phi_n(j)*phi_ss(j)...
                    - beta_s(j)*(phi_n(j)^2 + phi_s(j)^2));
                %D^2z/Dt^2 = sinbeta(phi_ts+phi_s*phi_ss+phi_n*phi_ns)
                %    + cos(
                %           - phi_tn - phi_s*phi_ns
                %           + phi_n*phi_ss + beta_s(phi_n^2 + phi_s^2)
                %         )
                d2zdt2(j)=sinbeta(j)*(phi_ts(j)+phi_s(j)*phi_ss(j) ...
                    + phi_n(j)*phi_ns(j))...
                    + cosbeta(j)*(-phi_tn(j) - phi_s(j)*phi_ns(j) ...
                    + phi_n(j)*phi_ss(j) ...
                    + beta_s(j)*(phi_n(j)^2 + phi_s(j)^2));
            end

            obj.characteristics.phi_ts = phi_ts;
            obj.characteristics.d2xdt2 = d2xdt2;
            obj.characteristics.d2zdt2 = d2zdt2;
            obj.characteristics.d2phidt2 = d2phidt2;
        end

        function on_step(obj,dt,order)
            
            %TODO different timestep methods; this is taylor2
            obj.boundary_nodes(:,1) = obj.boundary_nodes(:,1) +...
                dt.*obj.characteristics.dxdt'...
                + ((dt^2)/2).*obj.characteristics.d2xdt2';
            obj.boundary_nodes(:,2) = obj.boundary_nodes(:,2) +...
                dt.*obj.characteristics.dzdt'...
                + ((dt^2)/2).*obj.characteristics.d2zdt2';
            obj.characteristics.phi = obj.characteristics.phi +...
                dt.*obj.characteristics.dphidt...
                + ((dt^2)/2).*obj.characteristics.d2phidt2;
        end

        
    end
end

