function [E,KE,PE] = calc_energy(obj)
%Calculates the total, kinetic, and potential energy in the free surface


KE_int = cell(1,length(obj.boundaries));
PE_int = cell(1,length(obj.boundaries));
for i=1:length(KE_int)
    bdry = obj.boundaries{i};
    if isa(bdry,"free_surface_bdry")
        KE_int{i} = bdry.characteristics.phi .* bdry.characteristics.phi_n;
        PE_int{i} = bdry.boundary_nodes(:,2).^2;
    else
        KE_int{i} = 0;
        PE_int{i} = 0;
    end
end



KE = bem_integrate(obj,KE_int,0) / 2;
PE = bem_integrate(obj,PE_int,1) * (-obj.meta.g) / 2;
E = KE + PE;

end

