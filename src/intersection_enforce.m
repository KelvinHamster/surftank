function intersection_enforce(CW)
%resolves the intersection requirement between bdry and its CCW neighbor

inter = CW.intersection_CCW;
if strcmp(inter.resolution,"ignore")
    return
end
CCW = inter.CCW_link;

if strcmp(inter.resolution,"midpoint")
    %place doublenode at the average between the two nodes
    if (~strcmp(CW.formulation.type,"nodes"))...
            || (~strcmp(CCW.formulation.type,"nodes"))
        error("Boundary intersection: type 'midpoint' presently" + ...
            " requires 'nodes' formulation type for both sides.")
    end
    r = (CW.boundary_nodes(end,:) + CCW.boundary_nodes(1,:))/2;
    CW.boundary_nodes(end,:) = r;
    CCW.boundary_nodes(1,:) = r;


elseif strcmp(inter.resolution,"match_CW")
    %place doublenode on the node on the clockwise side
    if (~strcmp(CCW.formulation.type,"nodes"))
        error("Boundary intersection: type 'match_CW' presently" + ...
            " requires 'nodes' formulation type for CCW side.")
    end
    CCW.boundary_nodes(1,:) = CW.boundary_nodes(end,:);


elseif strcmp(inter.resolution,"match_CCW")
    %place doublenode on the node on the counter-clockwise side
    if (~strcmp(CW.formulation.type,"nodes"))
        error("Boundary intersection: type 'match_CCW' presently" + ...
            " requires 'nodes' formulation type for CW side.")
    end
    CW.boundary_nodes(end,:) = CCW.boundary_nodes(1,:);

    
elseif strcmp(inter.resolution,"curve_intersect")
    fprintf("'curve_intersect' corner resolution is not known to work."+...
        " Probably avoid using it.");
    %place doublenode at the intersection between the interpolation of the
    %two curves; utilize Newton's method with a linear search to obtain a
    %solution. First, we want to get a multivar function to minimize
    [func_CW, deriv_CW, d2_CW, init_CW] = ...
        curve_intersect_form_pos_func(CW,0);
    [func_CCW, deriv_CCW, d2_CCW, init_CCW] = ...
        curve_intersect_form_pos_func(CCW,1);

    %minimize dist: f = sum((CCW - CW)^2)
    f = @(s,t) sum((func_CCW(t) - func_CW(s)).^2);
    f_s = @(s,t) (-2) * dot((func_CCW(t) - func_CW(s)),deriv_CW(s));
    f_t = @(s,t) 2 * dot((func_CCW(t) - func_CW(s)),deriv_CCW(t));

    
    f_ss = @(s,t) (-2) * (dot((func_CCW(t) - func_CW(s)),d2_CW(s))...
        -dot(deriv_CW(s),deriv_CW(s)));
    f_st = @(s,t) (-2) * dot(deriv_CW(s),deriv_CCW(t));
    f_tt = @(s,t) 2 * (dot((func_CCW(t) - func_CW(s)),d2_CCW(t)) ...
        + dot(deriv_CCW(t),deriv_CCW(t)));
    
    s0 = init_CW; t0 = init_CCW; fx = f(s0,t0);
    %modified newton's (min f): normally x1 = x0 - (D^2f(x))^(-1) * Df(x)
    %change so dir = - (D^2f(x))^(-1) * Df(x) (or -Df(x) if D^2f is close
    % to singular)
    %then use linear search
    ftol = 1e-10; % f=0 is the lowest; how close are we?
    dftol2 = 1e-6; % df = 0 is a fixed point; if |dir|^2 < dftol2, kill
    max_iters = 100;
    iter = 0;
    alpha = 1; %initial step size
    while iter < max_iters && fx > ftol
        iter = iter + 1;
        fss = f_ss(s0,t0); fst = f_st(s0,t0); ftt = f_tt(s0,t0);
        fs = f_s(s0,t0); ft = f_t(s0,t0);

        deter = fss * ftt - fst^2;
        if abs(deter) < 1e-4
            ds = fst * ft - ftt * fs;
            dt = fst * fs - fss * ft;

            %make sure we are heading towards a minimum,
            if ds*fs + dt*ft > dftol2
                ds = - fs;
                dt = - ft;
            end

        else
            ds = - fs;
            dt = - ft;
        end

        if ds*ds + dt*dt < dftol2
            break
        end
        %we have step dir. choose size
        alpha = alpha * 2;
        i = 1; MAX_LIN_ITERS = 100;
        while i < MAX_LIN_ITERS && f(s0+alpha*ds,t0+alpha*dt) >= fx
            alpha = alpha * 0.5;
            i = i+1;
        end
        
        s0 = s0+alpha*ds; t0 = t0+alpha*dt; fx = f(s0,t0);
    end

    %this should be very close to the true intersection, so just midpoint.
    r = 0.5 * (func_CW(s0) + func_CCW(t0));
    CW.boundary_nodes(end,:) = r;
    CCW.boundary_nodes(1,:) = r;
    if strcmp(CW.formulation.type,'parameterized')
        CW.formulation.param_CCW = s0;
    elseif strcmp(CW.formulation.type,'nodes')
        %update phi or phi_n (for now, only supported for 'nodes')
        M = CW.BIE_shape_interp.M;
        Lmat = CW.BIE_shape_interp.Lagrange;
        if strcmp(CW.cond_type,'dirichlet')
            if isfield(CW.characteristics,'phi')
                CW.characteristics.phi(end) = (s0.^(0:M) * Lmat') *...
                    CW.characteristics.phi((end-M):end)';
            end
            if isfield(CW.characteristics,'phi_t')
                CW.characteristics.phi_t(end) = (s0.^(0:M) * Lmat') *...
                    CW.characteristics.phi_t((end-M):end)';
            end
        else
            if isfield(CW.characteristics,'phi_n')
                CW.characteristics.phi_n(end) = (s0.^(0:M) * Lmat') *...
                    CW.characteristics.phi_n((end-M):end)';
            end
            if isfield(CW.characteristics,'phi_tn')
                CW.characteristics.phi_tn(end) = (s0.^(0:M) * Lmat') *...
                    CW.characteristics.phi_tn((end-M):end)';
            end
        end
    end
    if strcmp(CCW.formulation.type,'parameterized')
        CCW.formulation.param_CW = t0;
    elseif strcmp(CCW.formulation.type,'nodes')
        %update phi or phi_n (for now, only supported for 'nodes')
        M = CCW.BIE_shape_interp.M;
        Lmat = CCW.BIE_shape_interp.Lagrange;
        if strcmp(CCW.cond_type,'dirichlet')
            if isfield(CCW.characteristics,'phi')
                CCW.characteristics.phi(end) = (t0.^(0:M) * Lmat') *...
                    CCW.characteristics.phi(1:(M+1))';
            end
            if isfield(CCW.characteristics,'phi_t')
                CCW.characteristics.phi_t(end) = (t0.^(0:M) * Lmat') *...
                    CCW.characteristics.phi_t(1:(M+1))';
            end
        else
            if isfield(CCW.characteristics,'phi_n')
                CCW.characteristics.phi_n(end) = (t0.^(0:M) * Lmat') *...
                    CCW.characteristics.phi_n(1:(M+1))';
            end
            if isfield(CCW.characteristics,'phi_tn')
                CCW.characteristics.phi_tn(end) = (t0.^(0:M) * Lmat') *...
                    CCW.characteristics.phi_tn(1:(M+1))';
            end
        end
    end
    

else
    error("Boundary intersection: not a supported resolution" + ...
        " type (%s)",inter.resolution);
end

end

function [func,deriv,d2,init] = ...
            curve_intersect_form_pos_func(bdry,side)
    %gives a function used to solve the curve intersection problem
    % 'side' is 0 for CW side, and 1 for CCW side. Returned values:
    % func = parameterized boundary
    % deriv = derivative of func
    % d2 = 2nd derivative of func
    % init = start parameter, corner node right now is func(init)
    if side
        intersect = bdry.intersect_CW_link.intersection_CCW;
    else
        intersect = bdry.intersection_CCW;
    end
    if strcmp(bdry.formulation.type,'nodes')
        if side
            if isfield(intersect,'CCW_interp')
                interp = intersect.CCW_interp;
            else
                interp = get_interpolation_struct(3,0,1);
            end
        else
            if isfield(intersect,'CW_interp')
                interp = intersect.CW_interp;
            else
                interp = get_interpolation_struct(3,0,1);
            end
        end
        M = interp.M;
        Lmat = interp.Lagrange;
        Lpmat = interp.Lagrange_prime;
        Lppmat = interp.Lagrange_primeprime;
        
        if side
            seg_nodes = bdry.boundary_nodes(1:(M+1),:);
            init = 0;
        else
            seg_nodes = bdry.boundary_nodes((end-M):end,:);
            init = 1;
        end
        func = @(t) (t.^(0:M) * Lmat') * seg_nodes;
        deriv = @(t) (t.^(0:M-1) * Lpmat') * seg_nodes;
        d2 = @(t) (t.^(0:M-2) * Lppmat') * seg_nodes;
    elseif strcmp(bdry.formulation.type,'parameterized')
        func = bdry.formulation.surface_defn;
        EPS = 1e-5;
        if isfield(bdry.formulation,'surface_defn_deriv')
            deriv = bdry.formulation.surface_defn_deriv;
        else
            deriv = @(t) (func(t + EPS) - func(t - EPS))/(2*EPS);
        end
        if isfield(bdry.formulation,'surface_defn_deriv2')
            d2 = bdry.formulation.surface_defn_deriv2;
        else
            d2 = @(t) (func(t + EPS) + func(t - EPS) - 2*func(t))/(EPS^2);
        end
        if side %==1: CCW
            if isfield(bdry.formulation,'param_CCW')
                init = bdry.formulation.param_CCW;
            else
                init = 0;
            end
        else %side == 0: CW
            if isfield(bdry.formulation,'param_CW')
                init = bdry.formulation.param_CW;
            else
                init = 1;
            end
        end
    else
        error("'curve_intersect' doublenode resolution does not " + ...
            "support boundary formulation type '%s'.",CW.bd_type);
    end
end
