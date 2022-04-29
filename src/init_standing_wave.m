%initializes a standing wave simulation in the linear domain

%gravity
g = 9.81;

%5
if ~exist('wall_res_factor','var')
    wall_res_factor = 2.5;
end

%height of the wave
if (~exist('a0','var')) || strcmp(a0,'default')
    a0 = 0.0001;
end

%length of the pool
if ~exist('L','var') || strcmp(L,'default')
    L = 10;
end
%depth of the pool
if ~exist('h','var')
    h = 1;
end
%number of nodes in the standing wave
N = 2;

%element size
M = 4;


%surface points
if ~exist('Nx','var')
    Nx = 100;
end

%time step
if ~exist('dt','var')
    dt = 0.01;
end

%===========

k = pi*N/L; omega = sqrt(g*k*tanh(k*h));

x = linspace(0,L,Nx);

dx = x(2)-x(1);

eta = a0.*cos(k.*x);

%phi = -2a0*omega/k * (cosh(k(h+z))/sinh(kh)) * cos(kx)sin(wt) <- 0 at t=0

%===============
%BOTTOM: floor; no movement into or out of floor

floorpoints = M*round(Nx/M) + 1;

h0 = h;
h = h*ones(1,floorpoints);

xb = linspace(0,L,floorpoints);
yb = -h;
doublenode = zeros(1,floorpoints);
doublenode(floorpoints) = 1;

%RIGHT: wall

wallpoints = M*round((wall_res_factor*(h0-a0-dx)/dx)/M); %we include the last point later (to complete the element)

y0 = linspace(-h0,-a0-dx, wallpoints);
xb = [xb L*ones(1,wallpoints)];
yb = [yb y0];
doublenode = [doublenode zeros(1,wallpoints+1)];
doublenode(length(xb)+1) = 1;


%TOP: wave/eta, remembering to doublenode
FS_i = length(xb) + 2;

xb = [xb x(Nx) flip(x) x(1)];
yb = [yb eta(Nx) flip(eta) eta(1)];
doublenode = [doublenode zeros(1,Nx+1)];
doublenode(length(xb) - 1) = 1;


FS_f = length(xb) - 1;

%LEFT: wall; just go downward and use the zero-tangential-movement
%condition

y0 = linspace(-h0,-a0-dx, wallpoints);
xb = [xb zeros(1,length(y0))];
yb = [yb fliplr(y0)];

xb = [xb 0];
yb = [yb -h0];
doublenode = [doublenode zeros(1,length(y0))];
doublenode(length(doublenode)) = 1;


%transfer to the correct variables and pass it on
x = xb;
z = yb;

phi_FS = zeros(1,FS_f - FS_i + 1);


generated = 1;
t = 0;


pool_length = L;
wavenumber = k;