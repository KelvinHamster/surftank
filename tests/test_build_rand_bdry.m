function bdry = test_build_rand_bdry(xf,zf,nodecount_target)
%TEST_BUILD_RAND_BDRY Builds a random boundary given the function handles
%   xf and zf should be functions defining the boundary on the domain [0,1]

if ~exist('dt_target','var')
    nodecount_target = 100;
end

M_MAX = 1; %CHANGE THIS LATER!!! ======================================
m = randi(M_MAX);
func_ord = randi(M_MAX);
fac = lcm(m,func_ord);
numpts = ceil(nodecount_target/fac)*fac + 1;

T = linspace(0,1,numpts);

x = xf(T);
z = zf(T);

%bdry interp rule
bdrule = rand();
if bdrule < 0.33
    bdrule = 'sliding';
elseif bdrule < 0.66
    bdrule = 'sequential';
else
    bdrule = 'spline';
end

frule = rand();
if frule < 0.5
    frule = 'sliding';
else
    frule = 'sequential';
end

bdry = solid_bdry('surface_nodes',[x;z]','bdry_interp_order',m,...
    'BIE_interp_order',func_ord,...
    'bdry_interp_rule',bdrule,'func_interp_rule',frule);
end

