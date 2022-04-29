function V = calc_volume(obj)
%Calculates the volume inside the boundary.


V = -bem_integrate(obj,obj.boundary.z,1);

end

