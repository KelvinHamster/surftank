function [Kd,Kn] = test_calc_reference_integral(sing_x,sing_y,gamma_x,gamma_y,gammap_x,gammap_y,N,a,b)
%TEST_CALC_REFERENCE_INTEGRAL compute Kn and Kd along a given segment
%   Use an expensive/more comprehensive quadrature to compute the integrals
%   (N G dgamma) and (N Gn dgamma) for given function handles.
%   given values (except sing,a,b) are function handles, which will be
%   integrated from a to b
G_dgam = @(t) log((gamma_x(t) - sing_x).^2 + (gamma_y(t) - sing_y).^2)...
    .*sqrt(gammap_x(t).^2 + gammap_y(t).^2)/(-4*pi);
Gn_dgam= @(t) ((gamma_x(t) - sing_x).*gammap_y(t) - (gamma_y(t) - sing_y).*gammap_x(t)) ...
    ./((-2*pi) .* ((gamma_x(t) - sing_x).^2 + (gamma_y(t) - sing_y).^2));
Kd = quadgk(@(t) N(t) .* G_dgam(t),a,b,'AbsTol',1e-6);
Kn = quadgk(@(t) N(t) .* Gn_dgam(t),a,b,'AbsTol',1e-6);

% if Kd == inf || Kn == inf || isnan(Kn) || isnan(Kd)
%     keyboard
% end
end