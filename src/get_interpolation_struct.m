function s = get_interpolation_struct(M,x0,xm)
%Returns a structure as in bem_sim interpolation parameters given the order and bounds.
%   Takes the order M, and location of the first+last nodes (x0,xm respectively)
%   The returned struct has
%   M = order (equal to M passed in)
%   Lagrange = matrix of lagrange polynomial coefficients
%       L_j(x)  = sum(x^k * Lagrange(j+1,k+1)) from k=0 to k=M } where
%       L_j(x_i)= kronecker_delta(j,i), x_i = x0*(1 - i/M)+ xm*(i/M)
%
%   Lagrange_prime = matrix of the derivatives of lagrange polys
%       L_j'(x) = sum(x^k * Lagrange_prime(j+1,k+1)) from k=0 to k=M-1
%
%   Lagrange_primeprime = matrix of the 2nd derivatives of lagrange polys
%       L_j''(x)= sum(x^k * Lagrange_primeprime(j+1,k+1)) from k=0 to k=M-2
%
%
%   For bem_sim, x0 should always be zero, and xm should be 1 for
%   obj.interpolation and obj.interpolation_sliding. For
%   obj.interpolation_FS, though, xm should be equal to M
nodes = linspace(x0,xm,M+1);

Lmat = zeros(M+1);
Lpmat = zeros(M+1,M);
Lppmat = zeros(M+1,M-1);

%apparently there is a polynomial library
for j=1:M+1
    p = poly(nodes(1:(M+1) ~= j)); %poly with roots x_i for   i ~= (j-1)
    %normalize around p(x_i) = 1, and flip
    p = flip(p/polyval(p,nodes(j)));
    %L_j(x) = sum( x^(k-1) * p(k) ), k=1:M+1
    Lmat(j,:) = p;
    %L_j'(x) = sum( (k-1) * x^(k-2) * p(k  ) ), k=2:M+1
    %L_j'(x) = sum( (k  ) * x^(k-1) * p(k+1) ), k=1:M
    Lpmat(j,:) = p(2:M+1) .* (1:M);
    %L_j''(x) = sum( (k  )* (k-1) * x^(k-2) * p(k+1) ), k=2:M
    %L_j''(x) = sum( (k+1)* (k  ) * x^(k-1) * p(k+2) ), k=1:M-1
    Lppmat(j,:) = p(3:M+1) .* (1:M-1) .* (2:M);
end


s = struct('M',M, 'Lagrange', Lmat, 'Lagrange_prime', Lpmat, 'Lagrange_primeprime',Lppmat);
end

