function sim = test_build_random_pieslicedomain()
%TEST_BUILD_RANDOM_PIESLICEDOMAIN Builds a sector domain

R = 10; %radius
t0 = rand() * 2*pi; %starting angle
phi = (0.1 + 0.8*rand()) * 2*pi; %how large (angle)

sim = bem_sim({...
    test_build_rand_bdry(@(t) cos(t0)*R*t,@(t) sin(t0)*R*t),...
    test_build_rand_bdry(@(t) cos(t0+phi*t)*R,@(t) sin(t0+phi*t)*R),...
    test_build_rand_bdry(@(t) cos(t0+phi)*R*(1-t),@(t) sin(t0+phi)*R*(1-t))});
end

