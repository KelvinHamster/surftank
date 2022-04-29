function [In,Id] = segment_In_Id_integrate_FS(ax,az,t1,t2,r_sing,...
        shape_offset,M_FS,Lmat_FS,Lpmat_FS,...
        cos_alpha_max, subint_min, do_MCI)
%Adaptive integration scheme for an integral on the specified segment
%   Recursive function to integrate a segment with a large intercept
%   angle. Binary subdivisions are made until the intercept angle on
%   each segment is less than the specified alpha_max value
%
%   ax,az       if do_MCI is true, these are the coefficients to the
%            spline, where
%            r(t) = (ax4,az4) + (ax3,az3)t + (ax2,az2)t^2 + (ax1,az1)t^3
%               otherwise, these are the x and z locations of the elements
%
%   t1,t2       the segment to integrate is on (t1,t2).
%
%   r_sing      (x,z) location of the singularity
%
%   shape_offset     shift of the shape function to the parent segment
%
%   M_FS, Lmat_FS    data for free surface lagrange interpolation (shape
%                  function)
%
%   cos_alpha_max    cosine of the largest allowable angle, determines
%                 the threshold for finding if a subdivision is necessary.
%
%   subint_min       smallest subiterval size. if t2-t1 <= subint_min, then
%                 we do not subdivide anymore
%
%   returns vectors In, Id which are indexed by nodes i, for i=1,...,M_FS+1
%

%is a subdivision allowed by our min-length restriction?
if t2-t1 > subint_min
    %yes. alpha on this segment: heuristic from dot product
    if do_MCI
        v1 = (t1.^(3:-1:0)) * [ax;az]' - r_sing;%vec from singularity to r(t1)
        v2 = (t2.^(3:-1:0)) * [ax;az]' - r_sing;%vec from singularity to r(t2)
    else
        v1 = [ax(shape_offset+1),az(shape_offset+1)] - r_sing;%vec from singularity to r(t1)
        v2 = [ax(shape_offset+2),az(shape_offset+2)] - r_sing;%vec from singularity to r(t2)
    end
    if dot(v1./norm(v1),v2./norm(v2)) < cos_alpha_max
        %smaller cosine, larger angle; subdivide
        %return sum of integrals over subdivided segments
        [In,Id] = segment_In_Id_integrate_FS(ax,az,t1,t1+(t2-t1)/2,r_sing,...
                shape_offset,M_FS,Lmat_FS,Lpmat_FS,cos_alpha_max, subint_min, do_MCI);
        [In2,Id2] = segment_In_Id_integrate_FS(ax,az,t1+(t2-t1)/2,t2,r_sing,...
                shape_offset,M_FS,Lmat_FS,Lpmat_FS,cos_alpha_max, subint_min, do_MCI);
        In = In + In2;
        Id = Id + Id2;
        return
    end
end
%evaluate the integral here

%QUADRATURE_PTS = 16;
QUAD_T = [0.0052995325041750333469603 0.0277124884633837102743126 0.0671843988060841224019271 0.1222977958224984867952045 0.1910618777986781147149031 0.2709916111713863151599924 0.3591982246103705422868302 0.4524937450811812866824368 0.5475062549188187688287144 0.6408017753896294577131698 0.7290083888286137403511589 0.8089381222013218852850969 0.8777022041775015548381589 0.9328156011939158220869217 0.9722875115366163001340283 0.9947004674958249692551249]';
QUAD_W = [0.0135762297058770482066636 0.0311267619693239468159351 0.0475792558412463928441127 0.0623144856277669384470030 0.0747979944082883679845608 0.0845782596975012679330064 0.0913017075224617918882686 0.0947253052275342510846201 0.0947253052275342510846201 0.0913017075224617918882686 0.0845782596975012679330064 0.0747979944082883679845608 0.0623144856277669384470030 0.0475792558412463928441127 0.0311267619693239468159351 0.0135762297058770482066636]';

QUAD_T = t1 + QUAD_T*(t2-t1);
% Id = -1/(4*pi) * dot( shape_quad(:,i).*log(sum((r_quad - r_sing).^2, 2)).*norm_rp_quad , QUAD_W);

shape_quad = ((QUAD_T+shape_offset).^(0:M_FS) * Lmat_FS');
if do_MCI
    r_quad = (QUAD_T.^(3:-1:0)) * [ax;az]';
    rp_quad = ((3:-1:1).* QUAD_T.^(2:-1:0)) * [ax(1:3);az(1:3)]';
else
    seg_nodes = [ax', az'];
    r_quad = shape_quad * seg_nodes;
    rp_quad = ((QUAD_T+shape_offset).^(0:M_FS-1) * Lpmat_FS') * seg_nodes;
end




Id = zeros(1,M_FS+1);
In = zeros(1,M_FS+1);
for i=1:M_FS+1
    %I_d_FS(fs_elem,l,i) = -1/(4*pi) * dot( shape_quad(:,i).*log(sum((r_quad - r_sing).^2, 2)).*norm_rp_quad , QUAD_W);
    Id(i) = -(t2-t1)/(4*pi) * dot( shape_quad(:,i).*log(sum((r_quad - r_sing).^2, 2)).*vecnorm(rp_quad,2,2) , QUAD_W);
    
    %I_n_FS(fs_elem,l,i) = -1/(2*pi) * integral(@(t) shape_func(t) .* dot(r(t)-r_sing,norm_func(t))./sum((r(t) - r_sing).^2),0,1,'ArrayValued',1,'AbsTol',INTEGRAL_TOL);
    In(i) = -(t2-t1)/(2*pi) * dot(shape_quad(:,i) .* dot(r_quad-r_sing, flip(rp_quad,2).*[1,-1],2)./sum((r_quad - r_sing).^2,2) , QUAD_W );
end

end

