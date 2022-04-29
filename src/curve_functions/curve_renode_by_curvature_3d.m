function V1 = curve_renode_by_curvature_3d(V,numnodes,a,b)
%Creates a 3xnumnodes matrix representing a curve reweighed so that the nodes are concentrated at regions of high curvature
% the distance between two nodes should be chosen such that the integral
% between them of:
%       min(a + curvature,b) ds
% should be roughly the same.
%
%
%   V - 3xN matrices, where each column represents a point in 3d space
%
%   mindist - the smallest distance between two nodes (this should be
%   smaller than the initial distances)
%
%   interpolation is done using splines.
v_pts = size(V,2);
dims = size(V,1);

%generate spline
v_spline = zeros(dims,v_pts-1,4);
for i = 1:dims
    v_spline(i,:,:) = csape(1:v_pts,V(i,:)).coefs;
end


%integrate segment "lengths" (by a + curvature "norm")
%Gauss Legendre 16: integrate on [0,1]
LEN_QUAD_T = [0.0052995325041750333469603 0.0277124884633837102743126 0.0671843988060841224019271 0.1222977958224984867952045 0.1910618777986781147149031 0.2709916111713863151599924 0.3591982246103705422868302 0.4524937450811812866824368 0.5475062549188187688287144 0.6408017753896294577131698 0.7290083888286137403511589 0.8089381222013218852850969 0.8777022041775015548381589 0.9328156011939158220869217 0.9722875115366163001340283 0.9947004674958249692551249]';
LEN_QUAD_W = [0.0135762297058770482066636 0.0311267619693239468159351 0.0475792558412463928441127 0.0623144856277669384470030 0.0747979944082883679845608 0.0845782596975012679330064 0.0913017075224617918882686 0.0947253052275342510846201 0.0947253052275342510846201 0.0913017075224617918882686 0.0845782596975012679330064 0.0747979944082883679845608 0.0623144856277669384470030 0.0475792558412463928441127 0.0311267619693239468159351 0.0135762297058770482066636]';


%how "long" is the segment between index k and k+1? This will give good
%initial guesses for where to match points to lengths
v_seglens = zeros(v_pts - 1, 1);

for k = 1:(v_pts-1)
    dr = squeeze(v_spline(:,k,1:3)) * ((LEN_QUAD_T.^(2:-1:0))' .* (3:-1:1)');
    d2r = squeeze(v_spline(:,k,1:2)) * ((LEN_QUAD_T.^(1:-1:0))' .* (3:-1:2)' .* (2:-1:1)');
    ds = vecnorm(dr,2,1);
    v_seglens(k) = dot(LEN_QUAD_W, ds .* min(b,a + vecnorm(cross(dr,d2r),2,1)./(ds.^3)));
end

v_cumlen = [0; cumsum(v_seglens)];

v1_cumlen = 0:numnodes-1;
v1_cumlen = v1_cumlen * (v_cumlen(v_pts) / v1_cumlen(numnodes));

V1 = zeros(dims,numnodes);
len_tol = 1e-4;

V1(:,1) = V(:,1); V1(:,numnodes) = V(:,v_pts);
for j = 2:(numnodes-1)
    l = v1_cumlen(j); %length target for V curve
    
    %target positions must lie on these segment indices
    k = find(v_cumlen <= l,1,"last");
    if k == numnodes
        k = k - 1;
    end


    %find point corresponding to position l1; use bisection method
    l = l - v_cumlen(k);
    left = 0; right = 1; c = 0.5;
    while right - left > len_tol

        dr = squeeze(v_spline(:,k,1:3)) * (((LEN_QUAD_T*c).^(2:-1:0))' .* (3:-1:1)');
        d2r = squeeze(v_spline(:,k,1:2)) * (((LEN_QUAD_T*c).^(1:-1:0))' .* (3:-1:2)' .* (2:-1:1)');
        ds = vecnorm(dr,2,1);
        fc = c* dot(LEN_QUAD_W, ds .* min(b,a + vecnorm(cross(dr,d2r),2,1)./(ds.^3)));
        if l < fc
            right = c; %we are overshooting
        else
            left = c; %we are undershooting
        end
        c = (left+right)/2;
    end
    V1(:,j) = squeeze(v_spline(:,k,:)) * (c.^(3:-1:0))';
end

end

