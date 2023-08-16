function [sim] = test_build_random_polydomain(n,max_width,max_height)
%TEST_BUILD_RANDOM_POLYDOMAIN generates a random polygonal domain
% Builds the polygon using test_grow_random_poly.
% takes optional arguments for number of sides (default = 10) and
% maximum allowable width and height (default=1)



if ~exist("n","var")
    n = 10;
end
if ~exist("max_width","var")
    max_width = 1;
end
if ~exist("max_height","var")
    max_height = 1;
end

P = test_grow_random_poly(n,max_width,max_height);
bdrys = {};
for k=1:n
    a = P(k-1 + n*(k==1),:);
    b = P(k,:);
    

    M_MAX = 1; %CHANGE THIS LATER!!! ======================================
    m = randi(M_MAX);
    func_ord = randi(M_MAX);
    fac = lcm(m,func_ord);
    maxpts = 20; %softmax
    numpts = randi(ceil(maxpts/fac))*fac + 1;

    x = linspace(a(1),b(1),numpts);
    z = linspace(a(2),b(2),numpts);

    %keyboard

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

    bdrys{k} = solid_bdry('surface_nodes',[x;z]','bdry_interp_order',m,...
        'BIE_interp_order',func_ord,...
        'bdry_interp_rule',bdrule,'func_interp_rule',frule);
end

sim = bem_sim(bdrys);
end

