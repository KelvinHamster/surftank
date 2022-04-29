function obj = sim_from_data(x,z,doublenode,phi_FS, FS_start,FS_end, varargin)
%creates a bem_sim class from the input initial condition.
%
%   x - x values of the boundary
%   z - z values of the boundary
%   doublenode - double node flag array
%   phi_FS - phi on the free surface
%   FS_start - start index of the free surface
%   FS_end - end index of the free surface
%   varargin - additional keyword arguments for bem_sim initialization.
%
%
% L is set to the range of values in x, and a0 is set to max(z).
if exist('varargin','var')
    obj = bem_sim('type','empty',varargin{:});
else
    obj = bem_sim('type','empty');
end
obj.boundary.N = length(x);
obj.boundary.x = x;
obj.boundary.z = z;
obj.boundary.doublenode = doublenode;
obj.boundary.FS_start = FS_start;
obj.boundary.FS_end = FS_end;
obj.boundary.phi_FS = phi_FS;


minx = min(obj.boundary.x);
maxx = max(obj.boundary.x);
obj.meta.plot_xlim = [minx-(maxx-minx)/20,maxx+(maxx-minx)/20];
obj.meta.L = maxx-minx;
obj.meta.a0 = max(obj.boundary.z);
end

