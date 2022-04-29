function result = bem_integrate(obj,value, mode, xl, zl)
%integrates the value given at each node
%   value is an array where value(i) is value at node i.
%   mode takes one of 4 values:
%       0 - default - integrate value d Gamma
%       1 - integrate value dx
%       2 - integrate value dz
%       3 - integrate value G((x,z),(x_l,z_l)) d Gamma
%       4 - integrate value G_n((x,z),(x_l,z_l)) d Gamma
%
% if mode = 3 or 4 is specified, xl and zl are used to locate the singularity.
%
% if only the free surface should be integrated, value can be set to 0 on
% non-free-surface nodes.

ADAPTIVE_THRESH = (10^-2);
%adaptive stopping conditions
cos_alpha_max = 0.7660;
subint_min = 0.01; %minimum length for subdivisions


if ~exist('mode','var')
    mode = 0;
end

obj = obj.update_characteristics();

if mode == 3 || mode == 4
    r_sing = [xl zl];
end

M = obj.interpolation.M;
Lmat = obj.interpolation.Lagrange;
Lpmat = obj.interpolation.Lagrange_prime;


M_FS = obj.interpolation_FS.M;
Lmat_FS = obj.interpolation_FS.Lagrange;
Lpmat_FS = obj.interpolation_FS.Lagrange_prime;

FS_i = obj.boundary.FS_start;
FS_f = obj.boundary.FS_end;
N = obj.boundary.N;
x = [obj.boundary.x obj.boundary.x(1)];
z = [obj.boundary.z obj.boundary.z(1)];
doublenode = obj.boundary.doublenode;

value = [value value(1)];


result = 0;


QUADRATURE_PTS = 16;
QUAD_T = [0.0052995325041750333469603 0.0277124884633837102743126 0.0671843988060841224019271 0.1222977958224984867952045 0.1910618777986781147149031 0.2709916111713863151599924 0.3591982246103705422868302 0.4524937450811812866824368 0.5475062549188187688287144 0.6408017753896294577131698 0.7290083888286137403511589 0.8089381222013218852850969 0.8777022041775015548381589 0.9328156011939158220869217 0.9722875115366163001340283 0.9947004674958249692551249]';
QUAD_W = [0.0135762297058770482066636 0.0311267619693239468159351 0.0475792558412463928441127 0.0623144856277669384470030 0.0747979944082883679845608 0.0845782596975012679330064 0.0913017075224617918882686 0.0947253052275342510846201 0.0947253052275342510846201 0.0913017075224617918882686 0.0845782596975012679330064 0.0747979944082883679845608 0.0623144856277669384470030 0.0475792558412463928441127 0.0311267619693239468159351 0.0135762297058770482066636]';

do_MCI = obj.meta.integral_spline;

r_quad = zeros(QUADRATURE_PTS,2);
rp_quad = zeros(QUADRATURE_PTS,2);
for k = FS_i:(FS_f-1)
    %free surface coefficients: a4 + a3*xi + a2*xi^2 + a1*xi^3
    ax = obj.characteristics_.FS_spline_x(k-FS_i+1,:);
    az = obj.characteristics_.FS_spline_z(k-FS_i+1,:);
    %r = @(t) [dot(ax,[t^3 t^2 t 1]) dot(az,[t^3 t^2 t 1])];%r
    %rp =@(t) [dot(ax,[3*t^2 2*t 1 0]) dot(az,[3*t^2 2*t 1 0])];%r'
    
    
    sliding_seg = max(FS_i, min(FS_f - M_FS, k - floor(M_FS/2)));
    shape_offset = k - sliding_seg;
    
    
    %check if the distance is small enough to warrant adaptive quadrature
    if mode == 3 || mode == 4
        r0 = [ax(4) az(4)];
        r1 = [dot(ax,ones(1,4)) dot(az,ones(1,4))];
        
        diff = r1 - r0;
        
        %project onto segment and see how far down the result is
        % seg_frac in [0,1] means that its on the segment
        length2 = sum(diff.^2);
        seg_frac = dot(diff, r_sing - r0)/length2;
        
        %capsule detection
        if (sum((r0-r_sing).^2) < ADAPTIVE_THRESH) || (sum((r1-r_sing).^2) < ADAPTIVE_THRESH) || (0 <= seg_frac && seg_frac <= 1 && (dot([1 -1].*flip(diff), r_sing - r0)/sqrt(length2))^2 < ADAPTIVE_THRESH)
            %integrals with shape functions - take linear combination
            if do_MCI
                [In,Id] = segment_In_Id_integrate_FS(ax,az,0,1,r_sing,shape_offset,M_FS,Lmat_FS,Lpmat_FS,...
                        cos_alpha_max, subint_min,1);
            else
                [In,Id] = segment_In_Id_integrate_FS(x(sliding_seg:sliding_seg+M_FS),...
                        z(sliding_seg:sliding_seg+M_FS),...
                        0,1,r_sing,shape_offset,M_FS,Lmat_FS,Lpmat_FS,cos_alpha_max, subint_min,0);
            end
            
            if mode == 3
                %Id
                result = result + dot(Id,value(sliding_seg:sliding_seg+M_FS));
                
            else
                %In
                result = result + dot(In,value(sliding_seg:sliding_seg+M_FS));
            end
            continue
        end
    end

    shape_quad = ((QUAD_T+shape_offset).^(0:M_FS) * Lmat_FS');
    shapep_quad = ((QUAD_T+shape_offset).^(0:M_FS-1) * Lpmat_FS');
    val_quad = shape_quad * value(sliding_seg:sliding_seg+M_FS)';
    
    
    
    
    %function at known sample points for quadrature expansion
    if do_MCI
        for i=1:QUADRATURE_PTS
            t = QUAD_T(i);
            r_quad(i,1) = dot(ax,[t^3 t^2 t 1]);
            r_quad(i,2) = dot(az,[t^3 t^2 t 1]);
            rp_quad(i,1) = dot(ax,[3*t^2 2*t 1 0]);
            rp_quad(i,2) = dot(az,[3*t^2 2*t 1 0]);
        end
    else
        seg_nodes = [x(sliding_seg:sliding_seg+M_FS)', z(sliding_seg:sliding_seg+M_FS)'];
        r_quad = shape_quad * seg_nodes;
        rp_quad = shapep_quad * seg_nodes;
    end
    norm_rp_quad = vecnorm(rp_quad,2,2);
    normal = fliplr(rp_quad).*[1 -1];
    
    switch(mode)
        case 0 %dGamma
            result = result + dot(val_quad .* norm_rp_quad, QUAD_W);
        case 1 %dx
            result = result + dot(val_quad .* rp_quad(:,1), QUAD_W);
        case 2 %dz
            result = result + dot(val_quad .* rp_quad(:,2), QUAD_W);
        case 3 %G dGamma
            result = result + -1/(2*pi) .* dot(val_quad .* norm_rp_quad .* log(vecnorm(r_quad - r_sing,2,2)), QUAD_W);
        case 4 %G_n dGamma
            result = result + -1/(2*pi) .* dot(val_quad .* (dot(r_quad - r_sing, normal,2)./sum((r_quad - r_sing).^2, 2)), QUAD_W);
    end
end

%non-free surface; use isoparametric elements
%for k = [1:M:(FS_i-1) FS_f:M:N]
k = 1;
while k <= N-M+1
    if doublenode(k)
        k = k+1;
        continue
    end
    if k >= FS_i && k < FS_f
        k = FS_f;
        continue
    end

    %Interpolate with lagrange: r = sum_{i=1}^{M+1} Li * r(k_i)
    x_elem = x(k:k+M)';
    z_elem = z(k:k+M)';
    val_elem = value(k:k+M)';
    %r = @(t) L(1:(M+1),t) * [x_elem z_elem];%r
    %rp =@(t) Lp(1:(M+1),t) * [x_elem z_elem];%r'
    
    %check if the distance is small enough to warrant adaptive quadrature
    if mode == 3 || mode == 4
        r0 = (0.^(0:M) * Lmat')*[x_elem z_elem];
        r1 = (1.^(0:M) * Lmat')*[x_elem z_elem];
        
        diff = r1 - r0;
        
        %project onto segment and see how far down the result is
        % seg_frac in [0,1] means that its on the segment
        length2 = sum(diff.^2);
        seg_frac = dot(diff, r_sing - r0)/length2;
        
        %capsule detection
        if (sum((r0-r_sing).^2) < ADAPTIVE_THRESH) || (sum((r1-r_sing).^2) < ADAPTIVE_THRESH) || (0 <= seg_frac && seg_frac <= 1 && (dot([1 -1].*flip(diff), r_sing - r0)/sqrt(length2))^2 < ADAPTIVE_THRESH)
            
            
            [In,Id] = segment_In_Id_integrate_nonFS(x_elem,z_elem,0,1,r_sing,...
                        M,Lmat,Lpmat,cos_alpha_max, subint_min);
            
            if mode == 3
                result = result + dot(Id,val_elem);
            else
                result = result + dot(In,val_elem);
            end
            k = k + M;
            continue
        end
    end
    
    
    %function at known sample points for quadrature expansion
    %r_quad = L(1:(M+1),QUAD_T)*[x_elem z_elem];
    r_quad = (QUAD_T.^(0:M) * Lmat')*[x_elem z_elem];
    %rp_quad = Lp(1:(M+1),QUAD_T)*[x_elem z_elem];
    rp_quad = (QUAD_T.^(0:(M-1)) * Lpmat')*[x_elem z_elem];
    norm_rp_quad = vecnorm(rp_quad, 2,2);
    
    normal = fliplr(rp_quad).*[1 -1]./norm_rp_quad;%n
    
    
   
    
    val_quad = (QUAD_T.^(0:M) * Lmat') * value(k:k+M)';
    
    switch(mode)
        case 0 %dGamma
            result = result + dot(val_quad .* norm_rp_quad, QUAD_W);
        case 1 %dx
            result = result + dot(val_quad .* rp_quad(:,1), QUAD_W);
        case 2 %dz
            result = result + dot(val_quad .* rp_quad(:,2), QUAD_W);
        case 3 %G dGamma
            result = result + -1/(2*pi) .* dot(val_quad .* norm_rp_quad .* log(vecnorm(r_quad - r_sing,2,2)), QUAD_W);
        case 4 %G_n dGamma
            result = result + -1/(2*pi) .* dot(val_quad .* norm_rp_quad .* (dot(r_quad - r_sing, normal,2)./sum((r_quad - r_sing).^2, 2)), QUAD_W);
    end
    k = k + M;
end




end

