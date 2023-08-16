function [flux] = calc_total_flux(obj)
%CALC_TOTAL_FLUX calculates dphi/dn integrated across the entire boundary

value = cell(1,length(obj.boundaries));
for i=1:length(value)
    bdry = obj.boundaries{i};
    value{i} = bdry.characteristics.phi_n;
end

flux = bem_integrate(obj,value,0);


end