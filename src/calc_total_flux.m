function [flux] = calc_total_flux(obj)
%CALC_TOTAL_FLUX calculates dphi/dn integrated across the entire boundary



flux_int = obj.characteristics_.phi_n;
flux = bem_integrate(obj,flux_int,0);


end