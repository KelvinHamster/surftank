function [err,deriv_err] = curve_error_3d(V1,V2,interp_scheme, int_segs)
%Calculates the error of two curves, defined through interpolation of points
% This function attempts to evaluate how good of a reparameterization a
% curve is. This is used to find a curve that spaces nodes more evenly than
% the initial parameterization.
%
%   V1,V2 - 3xN matrices, where each column represents a point in 3d space
%
%   interp_scheme - string to represent how V1 and V2, respectively,
%   should be interpolated. The options are:
%       'linear' - does linear interpolation
%       'sliding3' - 3rd order (4 node) lagrange interpolation
%       'sliding5' - 5rd order (6 node) lagrange interpolation
%       'spline' - does spline interpolation
%   
%   int_segs - number of integral segments for the composite quadrature
%
% The error integral is calculated using composite Gauss-Legendre
% quadrature, where the integrals are evenly spaced by curve length.

v1_pts = size(V1,2);
v2_pts = size(V2,2);
dims = max(size(V1,1),size(V2,1));

%handle scheme
mode = 0; do_sld_seg = 0;
if strcmp(interp_scheme,'linear')
    mode = 1;
    sld_seg = get_interpolation_struct(1,0,1);
    do_sld_seg = 1;
elseif strcmp(interp_scheme,'sliding3')
    mode = 3;
    sld_seg = get_interpolation_struct(3,0,3);
    do_sld_seg = 1;
elseif strcmp(interp_scheme,'sliding5')
    mode = 5;
    sld_seg = get_interpolation_struct(5,0,5);
    do_sld_seg = 1;
end

if mode == 0
    v1_spline = zeros(dims,v1_pts-1,4);
    v2_spline = zeros(dims,v2_pts-1,4);
    for i = 1:dims
        v1_spline(i,:,:) = csape(1:v1_pts,V1(i,:)).coefs;
    end
    for i = 1:dims
        v2_spline(i,:,:) = csape(1:v2_pts,V2(i,:)).coefs;
    end
end



%integrate segment lengths
%Gauss Legendre 16: integrate on [0,1]
LEN_QUAD_T = [0.0052995325041750333469603 0.0277124884633837102743126 0.0671843988060841224019271 0.1222977958224984867952045 0.1910618777986781147149031 0.2709916111713863151599924 0.3591982246103705422868302 0.4524937450811812866824368 0.5475062549188187688287144 0.6408017753896294577131698 0.7290083888286137403511589 0.8089381222013218852850969 0.8777022041775015548381589 0.9328156011939158220869217 0.9722875115366163001340283 0.9947004674958249692551249]';
LEN_QUAD_W = [0.0135762297058770482066636 0.0311267619693239468159351 0.0475792558412463928441127 0.0623144856277669384470030 0.0747979944082883679845608 0.0845782596975012679330064 0.0913017075224617918882686 0.0947253052275342510846201 0.0947253052275342510846201 0.0913017075224617918882686 0.0845782596975012679330064 0.0747979944082883679845608 0.0623144856277669384470030 0.0475792558412463928441127 0.0311267619693239468159351 0.0135762297058770482066636]';


%how long is the segment between index k and k+1? This will give good
%initial guesses for where to match points to lengths
v1_seglens = zeros(v1_pts - 1, 1);
v2_seglens = zeros(v2_pts - 1, 1);

for k = 1:(v1_pts-1)
    if mode == 1 %trivial case: linear
        v1_seglens(k) = sqrt(sum((V1(:,k+1) - V1(:,k)).^2));
    elseif mode == 0 %spline; a bit more difficult; derivative of poly
        derivs = squeeze(v1_spline(:,k,:)) * ((LEN_QUAD_T.^(2:-1:-1))' .* (3:-1:0)');
        v1_seglens(k) = dot(LEN_QUAD_W, vecnorm(derivs,2,1));
    elseif do_sld_seg %general sliding segments
        sld_pos = max(1, min(v1_pts - sld_seg.M, k - floor(sld_seg.M/2)));
        sld_off = k - sld_pos;
        derivs = ((sld_off + LEN_QUAD_T).^(0:sld_seg.M-1) * sld_seg.Lagrange_prime') * V1(:,sld_pos:(sld_pos + sld_seg.M))';
        v1_seglens(k) = dot(LEN_QUAD_W, vecnorm(derivs,2,2));
    end
end
for k = 1:(v2_pts-1)
    if mode == 1 %trivial case: linear
        v2_seglens(k) = sqrt(sum((V2(:,k+1) - V2(:,k)).^2));
    elseif mode == 0 %spline; a bit more difficult; derivative of poly
        derivs = squeeze(v2_spline(:,k,:)) * ((LEN_QUAD_T.^(2:-1:-1))' .* (3:-1:0)');
        v2_seglens(k) = dot(LEN_QUAD_W, vecnorm(derivs,2,1));
    elseif do_sld_seg %general sliding segments
        sld_pos = max(1, min(v2_pts - sld_seg.M, k - floor(sld_seg.M/2)));
        sld_off = k - sld_pos;
        derivs = ((sld_off + LEN_QUAD_T).^(0:sld_seg.M-1) * sld_seg.Lagrange_prime') * V2(:,sld_pos:(sld_pos + sld_seg.M))';
        v2_seglens(k) = dot(LEN_QUAD_W, vecnorm(derivs,2,2));
    end
end

v1_cumlen = [0; cumsum(v1_seglens)];
v2_cumlen = [0; cumsum(v2_seglens)];


%positions at target times for integration points, do our integral off of a gauss-legendre that is not so overkill as the length integrator.
len_tol = 1e-4;
%Gauss Legendre 3: integrate on [0,1]
QUAD_N = 3;
QUAD_T = (1+[-sqrt(3/5);0;sqrt(3/5)])/2;
QUAD_W = [5/18; 8/18; 5/18];

err = 0; % accumulator for error integrals
deriv_err = 0;
for i = 1:int_segs
    %integrate error on [(i-1)/int_segs, i/int_segs] of the way across both curves
    int_acc = 0; % integral accumulator
    deriv_acc = 0;
    for j = 1:QUAD_N
        l1 = (i-1 + QUAD_T(j))*v1_cumlen(v1_pts)/int_segs; %length target for V1 curve
        l2 = (i-1 + QUAD_T(j))*v2_cumlen(v2_pts)/int_segs; %length target for V2 curve
        
        %target positions must lie on these segment indices
        k1 = find(v1_cumlen <= l1,1,"last");
        k2 = find(v2_cumlen <= l2,1,"last");
        
        if do_sld_seg %sliding segment
            %find point corresponding to position l1; use bisection method
            l1 = l1 - v1_cumlen(k1);
            sld_pos = max(1, min(v1_pts - sld_seg.M, k1 - floor(sld_seg.M/2)));
            sld_off = k1 - sld_pos;
            a = 0; b = 1; c = 0.5;
            while b - a > len_tol
                derivs = ((sld_off + c*LEN_QUAD_T).^(0:sld_seg.M-1) * sld_seg.Lagrange_prime') * V1(:,sld_pos:(sld_pos + sld_seg.M))';
                fc = c * dot(LEN_QUAD_W, vecnorm(derivs,2,2));
                if l1 < fc
                    b = c; %we are overshooting
                else
                    a = c; %we are undershooting
                end
                c = (b+a)/2;
            end
            V1t = ((sld_off + c).^(0:sld_seg.M) * sld_seg.Lagrange') * V1(:,sld_pos:(sld_pos + sld_seg.M))';
            V1pt = ((sld_off + c).^(0:sld_seg.M-1) * sld_seg.Lagrange_prime') * V1(:,sld_pos:(sld_pos + sld_seg.M))';


            %find point corresponding to position l2; use bisection method
            l2 = l2 - v2_cumlen(k2);
            sld_pos = max(1, min(v2_pts - sld_seg.M, k2 - floor(sld_seg.M/2)));
            sld_off = k2 - sld_pos;
            a = 0; b = 1; c = 0.5;
            while b - a > len_tol
                derivs = ((sld_off + c*LEN_QUAD_T).^(0:sld_seg.M-1) * sld_seg.Lagrange_prime') * V2(:,sld_pos:(sld_pos + sld_seg.M))';
                fc = c * dot(LEN_QUAD_W, vecnorm(derivs,2,2));
                if l2 < fc
                    b = c; %we are overshooting
                else
                    a = c; %we are undershooting
                end
                c = (b+a)/2;
            end
            V2t = ((sld_off + c).^(0:sld_seg.M) * sld_seg.Lagrange') * V2(:,sld_pos:(sld_pos + sld_seg.M))';
            V2pt = ((sld_off + c).^(0:sld_seg.M-1) * sld_seg.Lagrange_prime') * V2(:,sld_pos:(sld_pos + sld_seg.M))';


        elseif mode == 0 %spline
            %find point corresponding to position l1; use bisection method
            l1 = l1 - v1_cumlen(k1);
            a = 0; b = 1; c = 0.5;
            while b - a > len_tol
                derivs = squeeze(v1_spline(:,k1,:)) * (((LEN_QUAD_T*c).^(2:-1:-1))' .* (3:-1:0)');
                fc = c * dot(LEN_QUAD_W, vecnorm(derivs,2,1));
                if l1 < fc
                    b = c; %we are overshooting
                else
                    a = c; %we are undershooting
                end
                c = (b+a)/2;
            end
            V1t = squeeze(v1_spline(:,k1,:)) * (c.^(3:-1:0))';
            V1pt = squeeze(v1_spline(:,k1,:)) * (((c).^(2:-1:-1))' .* (3:-1:0)');


            %find point corresponding to position l2; use bisection method
            l2 = l2 - v2_cumlen(k2);
            a = 0; b = 1; c = 0.5;
            while b - a > len_tol
                derivs = squeeze(v2_spline(:,k2,:)) * (((LEN_QUAD_T*c).^(2:-1:-1))' .* (3:-1:0)');
                fc = c * dot(LEN_QUAD_W, vecnorm(derivs,2,1));
                if l2 < fc
                    b = c; %we are overshooting
                else
                    a = c; %we are undershooting
                end
                c = (b+a)/2;
            end
            V2t = squeeze(v2_spline(:,k2,:)) * (c.^(3:-1:0))';
            V2pt = squeeze(v2_spline(:,k2,:)) * (((c).^(2:-1:-1))' .* (3:-1:0)');


        end
        %normalize V1t,V2t according to ds (length of boundary curve)
        V1pt = V1pt / sqrt(V1pt(1).^2 + V1pt(2).^2);
        V2pt = V2pt / sqrt(V2pt(1).^2 + V2pt(2).^2);

        int_acc = int_acc + sum((V1t-V2t).^2)*QUAD_W(j);
        deriv_acc = deriv_acc + sum((V1pt-V2pt).^2)*QUAD_W(j);
    end
    err = err + int_acc/int_segs;
    deriv_err = deriv_err + deriv_acc/int_segs;
end

end