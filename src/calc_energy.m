function [E,KE,PE] = calc_energy(obj)
%Calculates the total, kinetic, and potential energy in the free surface

KE_int = obj.characteristics_.phi.*obj.characteristics_.phi_n;
PE_int = obj.boundary.z.^2;

FS_i = obj.boundary.FS_start;
FS_f = obj.boundary.FS_end;

KE_int(1:FS_i - 1) = 0;
PE_int(1:FS_i - 1) = 0;
KE_int(FS_f + 1: length(KE_int)) = 0;
PE_int(FS_f + 1: length(PE_int)) = 0;

if obj.meta.has_a0
    a0 = obj.meta.a0;
else
    a0 = 1;
end

L = obj.meta.L;
g = obj.meta.g;


KE = bem_integrate(obj,KE_int,0) * 1020/(2*a0*L);%1/2
PE = bem_integrate(obj,PE_int,1) * -1020*g/(2*a0*L);%g/2
E = KE + PE;

end

