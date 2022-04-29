function [zbar] = calc_mean_water_level(obj)
%Calculates the mean water level by integrating z over the free surface
%   x,z represent the nodes of the boundary
%   FS_i,FS_f represent the indices of the first and last free surface
%   nodes, respectively.


FS_i = obj.boundary.FS_start;
FS_f = obj.boundary.FS_end;

zbar_int = obj.boundary.z;
zbar_int(1:FS_i - 1) = 0;
zbar_int(FS_f + 1: length(zbar_int)) = 0;

L = obj.boundary.x(FS_f) - obj.boundary.x(FS_i);

zbar = bem_integrate(obj,zbar_int,1)./L;

end