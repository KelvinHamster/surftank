function [zbar] = calc_mean_water_level(obj)
%Calculates the mean water level by integrating z over the free surface

L = 0;

value = cell(1,length(obj.boundaries));
for i=1:length(value)
    bdry = obj.boundaries{i};
    if isa(bdry,"free_surface_bdry")
        value{i} = bdry.boundary_nodes(:,2);
        L = L + (max(bdry.boundary_nodes(:,1))...
            - min(bdry.boundary_nodes(:,1)));
    else
        value{i} = 0;
    end
end


zbar = bem_integrate(obj,value,1)./L;

end