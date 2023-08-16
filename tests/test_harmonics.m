function maxerr = test_harmonics(sim,varargin)
%TEST_HARMONICS validates the integrator / poission solver to known fcns


% args
num_trials = 100;
test_KnKd = 1;
test_solver = 1;
verbosity = 1; % 0 - none; 1 - fails only; 2 - all
reltol = 1e-2;

p = inputParser;
addOptional(p,'num_trials',num_trials, @(x) ...
   isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1) == 0) );
addOptional(p,'test_KnKd',test_KnKd,@(x) islogical(x) || (x==0) || (x==1));
addOptional(p,'test_solver',test_solver,@(x) islogical(x) || (x==0) || (x==1));
addOptional(p,'verbosity',verbosity, @(x) ...
   isnumeric(x) && isscalar(x) && (mod(x,1) == 0) && (x >= 0) && (x <= 2));
addOptional(p,'reltol',reltol, @(x) ...
   isnumeric(x) && isscalar(x) && (x > 0));

parse(p, varargin{:});

num_trials = p.Results.num_trials;
test_KnKd = p.Results.test_KnKd;
test_solver = p.Results.test_solver;
verbosity = p.Results.verbosity;
reltol = p.Results.reltol;



BIE_integ_calc(sim);

alpha = 0;
for i=1:length(sim.boundaries)
    alpha = alpha - sum(sim.boundaries{i}.characteristics.K_n,2);
end
nodes = sim.get_nodelist();

minlen = max(vecnorm(nodes - circshift(nodes,1,1),2,2)); %largest dx

xmin = min(nodes(:,1));
xmax = max(nodes(:,1));
zmin = min(nodes(:,2));
zmax = max(nodes(:,2));
width = xmax - xmin;
height = zmax - zmin;
maxlen = min(width,height);

%harmonics:
%   linear: u = ax + bz + c
%   saddle: u = a(x^2 - z^2)
%   exp:    u = cos(ax+bz+c)*exp(az-bx+c)

maxerr = 0;

for trial_num=1:num_trials
    rng = rand();
    %center of u (in some sense)
    offx = xmin + width*rand();
    offz = zmin + height*rand();
    if rng < 0.1
        fcntype = "linear";
        a = (rand()-0.5)/width; b = (rand()-0.5) / height;
        u = @(x,z) a*(x-offx) + b*(z-offz);
    elseif rng < 0.2
        fcntype = "saddle";
        a = 0; b = 0;
        u = @(x,z) ((x-offx).^2 - (z-offz).^2)/(maxlen^2);
    else
        fcntype = "exp";
        %length scale in range [8*minlen,maxlen]
        r = 8*minlen + max(0, maxlen - 8*minlen)*rand();
        t = rand()*2*pi;
        a = cos(t)/r;
        b = sin(t)/r;
        u = @(x,z) cos(a*(x-offx)+b*(z-offz)) .*exp(a*(z-offz)-b*(x-offx));
    end
    arga = a; argb = b;

    err_KnKd = 0;
    err_solver = 0;
    solver_sol = [];
    u_inf = 0; %inf-norm of u
    q_inf = 0; %inf-norm of q
    for i=1:length(sim.boundaries)
        bdry = sim.boundaries{i};
        pts = bdry.boundary_nodes;
        gammap = pts(2:end,:) - pts(1:end-1,:);
        normal_ = flip(gammap,2) .* [1,-1] ./ vecnorm(gammap,2,2);
        normal = zeros(bdry.node_count,2);
        normal(2:end-1,:) = (normal_(1:end-1,:) + normal_(2:end,:))/2;
        normal(1,:) = normal_(1,:);
        normal(end,:) = normal_(end,:);
        
        %U = u(pts), Q = du/dn(pts)   by central difference approx
        EPS = 1e-6;
        U = u(pts(:,1),pts(:,2));
        Q = (u(pts(:,1)+normal(:,1)*EPS,pts(:,2)+normal(:,2)*EPS)...
            -u(pts(:,1)-normal(:,1)*EPS,pts(:,2)-normal(:,2)*EPS)...
            )/(2*EPS);
        bdry.characteristics.phi = U;
        bdry.characteristics.phi_t = U;
        bdry.characteristics.phi_n = Q;
        bdry.characteristics.phi_tn = Q;
        if test_KnKd
            err_KnKd = err_KnKd + bdry.characteristics.K_d*Q...
                -bdry.characteristics.K_n*U;
        end

        if strcmp(bdry.cond_type,'neumann')
            solver_sol = [solver_sol; U];
        elseif strcmp(bdry.cond_type,'dirichlet')
            solver_sol = [solver_sol; Q];
        end
        u_inf = max(u_inf, max(abs(U)));
        q_inf = max(q_inf, max(abs(Q)));
    end

    if test_KnKd
        %alpha terms; err is now the BIE error
        err_KnKd = err_KnKd - alpha.*u(nodes(:,1),nodes(:,2));
        err_KnKd = max(abs(err_KnKd));
    end
    if test_solver
        [A,b] = matrix_build(sim);
        err_solver = A\b - solver_sol;
        err_solver = max(abs(err_solver));
    end
    
    u_h1inf = max(u_inf,q_inf);
    maxerr = max(maxerr,max(err_KnKd,err_solver)/u_h1inf);
    if verbosity > 0 && (verbosity == 2 ||...
            (err_KnKd > reltol*u_h1inf ||...
            err_solver > reltol*u_h1inf))
        fprintf("Equation %d type: %s (a = %15.8e,b = %15.8e," + ...
            "x0 = %15.8f, z0 = %15.8f)\n",...
            trial_num,fcntype,arga,argb,offx,offz);
        fprintf("    max(u) = %15.8e\n",u_inf)
        fprintf("    max(q) = %15.8e\n",q_inf)
        if verbosity == 2 || err_KnKd > reltol*u_h1inf
            fprintf("     BIE error: %15.8e (rel = %15.8e)\n",...
                err_KnKd,err_KnKd/u_h1inf);
        end
        if verbosity == 2 || err_solver > reltol*u_h1inf
            fprintf("  solver error: %15.8e (rel = %15.8e)\n",...
                err_solver,err_solver/u_h1inf);
        end
    end

end
end

