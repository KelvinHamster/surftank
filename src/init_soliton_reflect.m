%INIT_BATHYMETRY generates the boundary arrays for the BEM
%   Takes the following arguments (there are default values):
%   h - array of values for the bathymetry; spaced dx apart. The bathymetry
%       extends past the last value of h if length(h)*dx < L, keeping a
%       flat profile for that remaining length
%   L - length of the simulation
%   a0 - initial height of the soliton
%   x0 - position of the soliton, where x=0 is the left boundary
%   dx - x-spacing of both h and the boundary elements of the free surface
%   g - acceleration due to gravity
%



if ~exist('wall_res_factor','var')
    wall_res_factor = 2.5;
end

if ~exist('h','var')
    h=1;
end

if ~exist('L','var') || strcmp(L,'default')
    L = 100;
end

if ~exist('a0','var') || strcmp(a0,'default')
    a0 = 0.1;
end

if ~exist('x0','var')
    x0 = 20;
end

if ~exist('Nx','var')
    Nx = 200;
end

if ~exist('g','var')
    g = 9.81;
end

if ~exist('M','var')
    M = 4;
end

%x coordinates to deal with
x = linspace(0,L,Nx);

dx = x(2) - x(1);

%TOP: soliton; original code from Knowles and Yeh

h0 = h(1);
X0 = x0;
K0 = sqrt((3*a0)/(4*h0^3))*(1-5/8*(a0/h0)+71/128*(a0/h0)^2); %Wave Number

s2 = (sech(K0*((x-X0)))).^2;

eta = h0*((a0/h0)*s2-3/4*(a0/h0)^2*(s2-s2.^2)+(a0/h0)^3*(5/8*s2-151/80*s2.^2+101/80*s2.^3));
ta = tanh(K0*((x-X0)));
s2t = s2.*ta;
s4t = s2.^2.*ta;
ps = (((a0/h0)*sqrt(g*h0)/(sqrt((3*a0)/(4*h0^3))))*(ta+(a0/h0)*(5/24*ta-1/3*s2t+3/4*(1+eta/h0).^2.*s2t)+(a0/h0)^2*(-1257/3200*ta+9/200*s2t+6/25*s4t+(1+eta/h0).^2.*(-9/32*s2t-3/2*s4t)+(1+eta/h0).^4.*(-3/16*s2t+9/16*s4t))));

xb = fliplr(x);
yb = fliplr(eta);
phi_FS = fliplr(ps);
doublenode = zeros(1,length(xb));
doublenode(length(xb)) = 1;

%index of the first and last free surface points
FS_i = 1;
FS_f = length(x);

%LEFT: wall; just go downward and use the zero-tangential-movement
%condition
wall_1_size = M * round(wall_res_factor * length(-h(1):dx:0)./M) + 1;
wall_1_size = max(wall_1_size, M + 1);


y0 = linspace(-h(1), eta(1), wall_1_size);
y0 = y0(2:wall_1_size - 1);
xb = [xb xb(length(xb)) zeros(1,length(y0))];
yb = [yb yb(length(yb)) fliplr(y0)];
doublenode = [doublenode zeros(1,length(y0)+2)];
doublenode(length(doublenode)) = 1;


%BOTTOM: floor; no movement into or out of floor

floorsize = M*round(Nx./M)+1;

if floorsize < length(h)
    h = h(1:floorsize);
elseif floorsize > length(h)
    h = [h h(length(h))*ones(1,floorsize - length(h))];
end

xb = [xb xb(length(xb)) linspace(0,L,floorsize)];
yb = [yb -h(1) -h];
doublenode = [doublenode zeros(1,length(h))];
doublenode(length(doublenode)) = 1;


%RIGHT: wall

wall_2_size = M * round(wall_res_factor * length(-h(length(h)):dx:0)./M) + 1;
wall_2_size = max(wall_2_size, M + 1);


y0 = linspace(-h(length(h)),eta(length(eta)), wall_2_size);
y0 = y0(2:wall_2_size - 1);
xb = [xb xb(length(xb)) L*ones(1,length(y0))];
yb = [yb yb(length(yb)) y0];
doublenode = [doublenode zeros(1,length(y0)+2)];
doublenode(length(doublenode)) = 1;


%last point
xb = [xb L L];
yb = [yb eta(length(eta)) eta(length(eta))];



x = xb;
z = yb;

