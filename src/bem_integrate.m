function result = bem_integrate(obj,value, mode, r_sing)
%integrates the value given at each node
%   value is one of the following:
%       - a function handle value(x,z) returning the integrand at a given
%       point (x,z)
%       - a cell array describing the integrand along each boundary
%       (value{i} is one of the following to specify it along boundary i)
%           - array, where value{i}(j) is the integrand of node j
%           - function handle value{i}(x,z) giving the integrand at point
%             (x,z)
%           - a constant number for a constant integrand (good for
%           measuring length or excluding boundaries by using
%           constants 1 or 0 respectively)
%       function handles are only sampled at node points, so quickly
%       varying functions do not work well.
%   mode takes one of 4 values:
%       0 - default - integrate value d Gamma
%       1 - integrate value dx
%       2 - integrate value dz
%       3 - integrate value G((x,z),(x_l,z_l)) d Gamma
%       4 - integrate value G_n((x,z),(x_l,z_l)) d Gamma
%
% if mode = 3 or 4 is specified, r_sing (1x2 array: (x,z) ) is used to
% locate the singularity.
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


num_bdrys = length(obj.boundaries);

%========ensure value is a cell array to make code uniform
if isa(value,'function_handle') || (isnumeric(value) && length(value) == 1)
    f = value;
    value = cell(1,num_bdrys);
    for bd=1:num_bdrys
        value{bd} = f;
    end
end


QUADRATURE_PTS = 16;
QUAD_T = [0.0052995325041750333469603 0.0277124884633837102743126 0.0671843988060841224019271 0.1222977958224984867952045 0.1910618777986781147149031 0.2709916111713863151599924 0.3591982246103705422868302 0.4524937450811812866824368 0.5475062549188187688287144 0.6408017753896294577131698 0.7290083888286137403511589 0.8089381222013218852850969 0.8777022041775015548381589 0.9328156011939158220869217 0.9722875115366163001340283 0.9947004674958249692551249]';
QUAD_W = [0.0135762297058770482066636 0.0311267619693239468159351 0.0475792558412463928441127 0.0623144856277669384470030 0.0747979944082883679845608 0.0845782596975012679330064 0.0913017075224617918882686 0.0947253052275342510846201 0.0947253052275342510846201 0.0913017075224617918882686 0.0845782596975012679330064 0.0747979944082883679845608 0.0623144856277669384470030 0.0475792558412463928441127 0.0311267619693239468159351 0.0135762297058770482066636]';
%========integrate over each bdry
result = 0;
for bd=1:num_bdrys
    bdry = obj.boundaries{bd};
    n = bdry.node_count;
    x_b = bdry.boundary_nodes(:,1);
    z_b = bdry.boundary_nodes(:,2);

    f = value{bd};
    %========make f into an array
    if isnumeric(f) && all(size(f) == 1)
        %we have a constant
        if f == 0
            %constant zero; just skip this bdry
            continue;
        else
            f = zeros(1,n) + f;
        end
    elseif isa(f,'function_handle')
        f_ = f;
        f = zeros(1,n);
        for j=1:n
            f(j) = f_(x_b(j),z_b(j));
        end
    end

    %use similar integ scheme to BIE_integ_calc
    M = bdry.BIE_shape_interp.M;
    Lmat = bdry.BIE_shape_interp.Lagrange;
    Lpmat = bdry.BIE_shape_interp.Lagrange_prime;

    %boundary interpolation
    if isfield(bdry.formulation,'bdry_interp')
        M_bdry = bdry.formulation.bdry_interp.M;
        Lmat_bdry = bdry.formulation.bdry_interp.Lagrange;
        Lpmat_bdry = bdry.formulation.bdry_interp.Lagrange_prime;
    else
        M_bdry = 0;
        Lmat_bdry = 0;
        Lpmat_bdry = 0;
    end

    % boundary segments: we integrate individual pieces, where
    % each piece is the interval between two consecutive nodes
    bdry_sf_bdry_interp_rule = bdry.formulation.bdry_interp_rule;
    if strcmp(bdry_sf_bdry_interp_rule,'sliding')
        % sliding: center at elem (seg connecting two consecutive nodes),
        % one for each elem
        half_M = floor(M_bdry/2);
        seg_starts = max(1,min(n-M_bdry,(1:(n-1))-half_M));%start of interp

        %segment parameters; so we do not broadcast whole [x_b,z_b] arrs
        seg_params = zeros(length(seg_starts),M_bdry+1,2);
        for i=1:(M_bdry+1)
            seg_params(:,i,1) = x_b(seg_starts + i - 1);
            seg_params(:,i,2) = z_b(seg_starts + i - 1);
        end
    elseif strcmp(bdry_sf_bdry_interp_rule,'sequential')
        % sequential: line up segments consecutively
        if mod(n-1,M_bdry) ~= 0
            error("Boundary %d: sequential boundary rule (M=%d) "+ ...
                "is incompatible with %d nodes!",...
                bdry_index,M_bdry,n);
        end
        seg_starts = M_bdry * floor(((1:(n-1)) - 1)/M_bdry) + 1;

        %segment parameters; so we do not broadcast whole [x_b,z_b] arrs
        seg_params = zeros(length(seg_starts),M_bdry+1,2);
        for i=1:(M_bdry+1)
            seg_params(:,i,1) = x_b(seg_starts + i - 1);
            seg_params(:,i,2) = z_b(seg_starts + i - 1);
        end
    elseif strcmp(bdry_sf_bdry_interp_rule,'spline')
        splx = csape(1:n, x_b).coefs;
        splz = csape(1:n, z_b).coefs;
    end
    bdry_sf_func_interp_rule = bdry.formulation.func_interp_rule;

    for elem=1:(n-1)

%==========================================================================
        %shape functions replaced by value fcn
        if strcmp(bdry_sf_func_interp_rule,'sequential') ||...
                strcmp(bdry_sf_func_interp_rule,'sliding')
            %nonzero shape functions:
            %seg_starts(elem) to seg_starts(elem) + M
            if strcmp(bdry_sf_func_interp_rule,'sliding')
                half_M = floor(M/2);
                shape_start = max(1,min(n-M,elem-half_M));
            elseif strcmp(bdry_sf_func_interp_rule,'sequential')
                if mod(n-1,M)~=0
                    error("Boundary %d: sequential function rule (M=%d)"...
                    +" is incompatible with %d nodes!",...
                    bdry_index,M,n);
                end
                shape_start = M * floor((elem - 1)/M) + 1;
            end       
        else
            error("Boundary %d: Invalid function interpolation rule!",...
                bdry_index);
        end
%==========================================================================
        val_elem = zeros(M+1,1);
        val_elem(:) = f(shape_start:(shape_start + M));
        x_seg = zeros(1,M_bdry+1);
        z_seg = zeros(1,M_bdry+1);
        x_seg(:) = x_b(seg_starts(elem):(seg_starts(elem)+M_bdry));
        z_seg(:) = z_b(seg_starts(elem):(seg_starts(elem)+M_bdry));
    
        %first check (for mode 3,4) if we are near a singularity
        if mode == 3 || mode == 4
            r0 = [x_b(elem), z_b(elem)];
            r1 = [x_b(elem+1), z_b(elem+1)];
            
            diff = r1 - r0;
            
            %project onto segment and see how far down the result is
            % seg_frac in [0,1] means that its on the segment
            length2 = sum(diff.^2);
            seg_frac = dot(diff, r_sing - r0)/length2;
            
            %capsule detection
            if (sum((r0-r_sing).^2) < ADAPTIVE_THRESH) || ...
                    (sum((r1-r_sing).^2) < ADAPTIVE_THRESH) ||...
                    (0 <= seg_frac && seg_frac <= 1 &&...
                    (dot([1 -1].*flip(diff), r_sing - r0)...
                    /sqrt(length2))^2 < ADAPTIVE_THRESH)
                
                
                if strcmp(bdry_sf_bdry_interp_rule,'spline')
                    [In,Id] = segment_In_Id_integrate( ...
                        splx(elem,:),splz(elem,:),0,1,...
                        r_sing,elem-shape_start,0,0,0,M,...
                            Lmat,Lpmat,cos_alpha_max, subint_min,1);
                else
                    bdry_off = elem - seg_starts(elem);
                    [In,Id] = segment_In_Id_integrate( ...
                            x_seg, z_seg,...
                            bdry_off/M_bdry,(bdry_off+1)/M_bdry,...
                            r_sing,seg_starts(elem)-shape_start,...
                            M_bdry,Lmat_bdry,Lpmat_bdry,M,Lmat,Lpmat,...
                            cos_alpha_max, subint_min,0);
                end
                
                if mode == 3
                    result = result + dot(Id,val_elem);
                else
                    result = result + dot(In,val_elem);
                end
                continue
            end
        end

        % not adaptive, so just do regular quadrature:
        % we need interpolated f, boundary, and (for mode 3,4) Green's fcn
        % no special quadrature needed
        if strcmp(bdry_sf_bdry_interp_rule,'spline')
            a = [splx(elem,:)', splz(elem,:)'];
            r_quad = (QUAD_T.^(3:-1:0)) * a;
            rp_quad = ((3:-1:1).*QUAD_T.^(2:-1:0)) * a(1:3,:);
        else
            offset = elem - seg_starts(elem);
            seg_nodes = squeeze(seg_params(elem,:,:));
            QUAD_REPARAM = (offset + QUAD_T)./M_bdry;
            r_quad = (QUAD_REPARAM.^(0:M_bdry) * Lmat_bdry') * seg_nodes;
            rp_quad = (QUAD_REPARAM.^(0:M_bdry-1) * Lpmat_bdry')*seg_nodes;
        end
        val_quad = ((offset + QUAD_T/M).^(0:M) * Lmat') * val_elem;
        norm_rp_quad = vecnorm(rp_quad,2,2);

        switch(mode)
            case 0 %dGamma
                result = result + dot(val_quad .* norm_rp_quad, QUAD_W);
            case 1 %dx
                result = result + dot(val_quad .* rp_quad(:,1), QUAD_W);
            case 2 %dz
                result = result + dot(val_quad .* rp_quad(:,2), QUAD_W);
            case 3 %G dGamma
                result = result + -1/(2*pi) .* dot(val_quad .* norm_rp_quad...
                    .* log(vecnorm(r_quad - r_sing,2,2)), QUAD_W);
            case 4 %G_n dGamma
                normal = fliplr(rp_quad).*[1 -1]./norm_rp_quad;%n
                result = result + -1/(2*pi) .* dot(val_quad .* ...
                    (dot(r_quad - r_sing, normal,2)./...
                    sum((r_quad - r_sing).^2, 2)), QUAD_W);
        end
    end

end





end

