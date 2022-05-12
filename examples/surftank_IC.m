function sim = surftank_IC(Nx,x0,a0,theta,varargin)
%SURFTANK_IC Sets up a bem_sim object with an initial condition for an overturn.
%   Initializes bathymetry and a soliton. The left-most wall is
%   positioned at x = 0. The returned object may need manual adjustment of:
%       meta.sim_parallel_workers (default = 0)
%       stepping (default: dt=0.01, store_dt=0.25, timestep_method='abm4')
%       interpolation (default: iso = 1, FS = 1, sliding = 4)
%   
%   Arguments:
%       Nx - integer, positive
%           number of free surface nodes, counting the endpoints. Nodes are
%           evenly distributed across the free surface, having a distance,
%           when projected on to the x-axis, of L/(Nx-1), where L is the
%           total length of the simulation, L1+L2+L3+L4 (see optional args)
%       x0 - float
%           Position for the center of the soliton initial condition along
%           the x-axis.
%
%       a0 - float, positive
%           initial height of the soliton. This is the maximum value in
%           height above the flat surface that the wave takes at t = 0.
%
%       theta - float
%           the angle at which the wave should be travelling. This scales
%           the surf ranch bathymetry but not the wave to simulate the wave
%           shoaling at a non-perpendicular angle to the beach. 0
%           represents the wave going directly to the beach, and pi/2
%           represents the wave going parallel (values should be between
%           -pi/2 and pi/2, not inclusive).
%
%
%   Optional Arguments:
%       h0 - float, positive (default = 2.56)
%           this is the initial depth, at the left-most section (deep).
%
%       L0 - float, positive (default = 20.0)
%           length of the deep section. The slope begins at x = L0
%       
%       L1 - float, positive (default = 10.0)
%           length of the slope section in x. The slope ends at x = L0+L1,
%           at the start of the sandbar.
%
%       h1 - float, positive (default = 0.95)
%           depth of the sandbar. The slope can be found as s = (h0-h1)/L1
%       
%       L2 - float, positive (default = 5.8)
%           width of the sandbar.
%
%       L3 - float, positive (default = 1)
%           length of the second slope section in x. This section spans
%           from x = L0+L1+L2 to x = L0+L1+L2+L3
%
%       h2 - float, positive (default = 1.3)
%           depth of the trough section.
%
%       L4 - float, positive (default = 4)
%           length of the trough section. At the end of this section, the
%           right wall is placed, ignoring the beach.
%
%       g - float, positive (default = 9.81)
%           acceleration due to gravity.
%       M - integer, positive (default = 1)
%           expected isoparametric order. The non-free surface is fixed to
%           contain (nM + 1) nodes for integers n on continuous segments of
%           the boundary to allow elements of order M (M+1 nodes) to be
%           used.
%       boundary_dx - float, positive (default = L/(Nx-1))
%           desired node spacing on the non-free-surface boundary. The
%           exact value is different in order to meet the (nM+1) - node
%           constraint, but n is chosen as to obtain the closest value to
%           boundary_dx. The default value is the same dx as implied by
%           choice of Nx.
%       wall_res_factor - float, positive (default = 2)
%           spatial resolution scaling on the left and right wall
%           boundaries. spatial resolution is roughly increased by this
%           factor, but the (nM+1) - node constraint is enforced,
%           preventing exact scaling.
h0 = 2.56;
L0 = 20.0;
L1 = 10.0;
h1 = 0.95;
L2 = 5.80;
L3 = 1.00;
h2 = 1.30;
L4 = 4.00;
g  = 9.81;
M = 1;
boundary_dx = 'default';
wall_res_factor = 2;

p = inputParser;
addOptional(p,'h0',h0, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'L0',L0, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'L1',L1, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'h1',h1, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'L2',L2, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'L3',L3, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'h2',h2, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'L4',L4, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p, 'g', g, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p, 'M', M, @(x) isnumeric(x) && isscalar(x) && (x > 0) && mod(x,1) == 0 );
addOptional(p, 'wall_res_factor', wall_res_factor, @(x) isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p, 'boundary_dx', boundary_dx,...
    @(x) (ischar(x) && strcmp(x,'default')) || (isnumeric(x) && isscalar(x) && (x > 0)) );

parse(p, varargin{:});
h0 = p.Results.h0;
L0 = p.Results.L0;
L1 = p.Results.L1;
h1 = p.Results.h1;
L2 = p.Results.L2;
L3 = p.Results.L3;
h2 = p.Results.h2;
L4 = p.Results.L4;
g = p.Results.g;
M = p.Results.M;
wall_res_factor = p.Results.wall_res_factor;
boundary_dx = p.Results.boundary_dx;

%rotate cross section: scale by 1/cos(theta)
rot_factor = sec(theta);
L0 = L0 * rot_factor;
L1 = L1 * rot_factor;
L2 = L2 * rot_factor;
L3 = L3 * rot_factor;
L4 = L4 * rot_factor;

%set up free surface
L = L0+L1+L2+L3+L4;
x_FS = linspace(L,0,Nx);

K0 = sqrt((3*a0)/(4*h0^3))*(1-5/8*(a0/h0)+71/128*(a0/h0)^2); %Wave Number
s2 = (sech(K0*((x_FS-x0)))).^2;

eta = h0*((a0/h0)*s2-3/4*(a0/h0)^2*(s2-s2.^2)+(a0/h0)^3*(5/8*s2-151/80*s2.^2+101/80*s2.^3));
ta = tanh(K0*((x_FS-x0)));
s2t = s2.*ta;
s4t = s2.^2.*ta;
phi_FS = (((a0/h0)*sqrt(g*h0)/(sqrt((3*a0)/(4*h0^3))))*(ta+(a0/h0)*(5/24*ta-1/3*s2t+3/4*(1+eta/h0).^2.*s2t)+(a0/h0)^2*(-1257/3200*ta+9/200*s2t+6/25*s4t+(1+eta/h0).^2.*(-9/32*s2t-3/2*s4t)+(1+eta/h0).^4.*(-3/16*s2t+9/16*s4t))));

%set up the rest

if (ischar(boundary_dx) && strcmp(boundary_dx,'default'))
    boundary_dx = L/(Nx - 1);
end

desired_nodes = wall_res_factor * h0/boundary_dx;
lwallnodes = max(M,M*round(desired_nodes/M))+1;

desired_nodes = L0/boundary_dx;
deepnodes = max(M,M*round(desired_nodes/M))+1;

desired_nodes = sqrt(L1^2 + (h1-h0)^2)/boundary_dx;
slope1nodes = max(M,M*round(desired_nodes/M))+1;

desired_nodes = L2/boundary_dx;
barnodes = max(M,M*round(desired_nodes/M))+1;

desired_nodes = sqrt(L3^2 + (h1-h2)^2)/boundary_dx;
slope2nodes = max(M,M*round(desired_nodes/M))+1;

desired_nodes = L4/boundary_dx;
troughnodes = max(M,M*round(desired_nodes/M))+1;

desired_nodes = wall_res_factor * h2/boundary_dx;
rwallnodes = max(M,M*round(desired_nodes/M))+1;


x = [x_FS, linspace(0,  0,lwallnodes), linspace(0  , L0,deepnodes), linspace(L0 ,L0+L1,slope1nodes), linspace(L0+L1,L0+L1+L2,barnodes),...
        linspace(L0+L1+L2,L0+L1+L2+L3,slope2nodes), linspace(L0+L1+L2+L3,L  ,troughnodes), linspace(L  ,L,rwallnodes)];
z = [eta,  linspace(0,-h0,lwallnodes), linspace(-h0,-h0,deepnodes), linspace(-h0,-h1  ,slope1nodes), linspace(-h1  ,     -h1,barnodes),...
        linspace(     -h1,  -h2      ,slope2nodes), linspace(-h2        ,-h2,troughnodes), linspace(-h2,0,rwallnodes)];

%enforce double node on ends of free surface
z(length(z)) = eta(1);
z(Nx+1) = eta(Nx);

doublenode = (x == circshift(x,-1)) & (z == circshift(z,-1));

%generate
sim = sim_from_data(x,z,doublenode,phi_FS,1,Nx);
sim.meta.g = g;
sim.meta.a0 = a0;
sim.meta.plot_ylim = [-0.1,1.5];
sim.meta.wall_resolution_factor = wall_res_factor;
end

