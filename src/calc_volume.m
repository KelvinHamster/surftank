function V = calc_volume(obj)
%Calculates the volume inside the boundary.


V = bem_integrate(obj,@(x,z) x, 2);

end

