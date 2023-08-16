% Falk Feddersen (c) 2001
%
% function that takes the radian wave frequency and
% a vector of depths and returns the wavenumber at that
% depth by solving the dispersion relationship
%
% function k = get_wavenum(omega,h)

function k = get_wavenum(omega,h)

% returns the wavenumber of the gravity wave
% dispersion relation, by using newtons method

% the initial guess will be the shallow water wavenumber

g = 9.81;

k = omega./sqrt(g*h);


f = g*k.*tanh(k.*h) - omega.^2;

while max(abs(f))>1e-10
  dfdk = g*k.*h.*(sech(k.*h)).^2 + g*tanh(k.*h);
  k = k - f./dfdk;
  f = g*k.*tanh(k.*h) - omega.^2;
end


  
