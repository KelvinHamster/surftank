function [x0,a0,eta_err,phi_err] = solitonfit(x0,a0,x_FS,z_FS, phi_FS)
x_FS = squeeze(x_FS);
z_FS = squeeze(z_FS);

X = fminsearch(@(X) solitonerr(X(1),X(2), x_FS, z_FS) , [x0; a0]);

x0 = X(1);
a0 = X(2);
if exist('phi_FS','var')
    [eta_err,phi_err] = both_errors(X(1),X(2), x_FS, z_FS, phi_FS);
else
    eta_err = solitonerr(X(1),X(2), x_FS, z_FS);
    phi_err = 0;
end
end



function err = solitonerr(x0,a0,x_FS,z_FS)
%finds the error off of the given data to the soliton with the specified parameters

h = 1; g = 9.81;

x = x_FS;

h0 = h;
X0 = x0;
K0 = sqrt((3*a0)/(4*h0^3))*(1-5/8*(a0/h0)+71/128*(a0/h0)^2); %Wave Number

s2 = (sech(K0*((x-X0)))).^2;

eta = h0*((a0/h0)*s2-3/4*(a0/h0)^2*(s2-s2.^2)+(a0/h0)^3*(5/8*s2-151/80*s2.^2+101/80*s2.^3));

err = norm(z_FS - eta, 2);% + norm(phi_FS - ps, 2);

end

function [eta,phi] = both_errors(x0,a0,x_FS,z_FS, phi_FS)

h = 1; g = 9.81;

x = x_FS;

h0 = h;
X0 = x0;
K0 = sqrt((3*a0)/(4*h0^3))*(1-5/8*(a0/h0)+71/128*(a0/h0)^2); %Wave Number

s2 = (sech(K0*((x-X0)))).^2;

eta = h0*((a0/h0)*s2-3/4*(a0/h0)^2*(s2-s2.^2)+(a0/h0)^3*(5/8*s2-151/80*s2.^2+101/80*s2.^3));ta = tanh(K0*((x-X0)));
s2t = s2.*ta;
s4t = s2.^2.*ta;
ps = (((a0/h0)*sqrt(g*h0)/(sqrt((3*a0)/(4*h0^3))))*(ta+(a0/h0)*(5/24*ta-1/3*s2t+3/4*(1+eta/h0).^2.*s2t)+(a0/h0)^2*(-1257/3200*ta+9/200*s2t+6/25*s4t+(1+eta/h0).^2.*(-9/32*s2t-3/2*s4t)+(1+eta/h0).^4.*(-3/16*s2t+9/16*s4t))));


eta = norm(z_FS - eta, 2);
phi = norm(phi_FS - ps, 2);
end