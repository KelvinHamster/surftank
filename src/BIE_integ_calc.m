function BIE_integ_calc(obj)
%BIE_INTEG_CALC Summary of this function goes here
%   Detailed explanation goes here

workers = obj.meta.parallel_workers;

%an interval is treated as quasisingular if one if the nodes is within
%sqrt(quasisingular_thresh) distance (in 2-norm) from the singularity
quasisingular_thresh = 5*10^-2;
cos_alpha_max = 0.7660;
subint_min = 0.0625; %minimum length for subdivisions
%tolerance for when to use L'Hospital's rule in the v-calculation
EPS = 10^-6;

%obtain x and z arrays; this will be used as the list of singularities
%we want to ignore the double nodes, since we are not examining the
%direction/angle of the path in this list.

%this is because we are forming a system of N equations (number of
%unique singularities) and 2*(N + num_doublenodes) variables. But
%the boundary conditions provide (var or equation for) half of these
%variables. The remaining num_doublenodes variables are accounted for by
%giving conditions between the two nodes that form each doublenode.
%The condition depends on BC (%dirichlet-dirichlet -> gradient matching,
%dirichlet-neumann or neumann-neumann -> phi equality)

[singularities,sing_indices] = obj.get_nodelist();
N = size(singularities,1);

%iterate each boundary, solving the integrals of the shape functions within
%the boundary
for bdry_index = 1:length(obj.boundaries)
    bdry = obj.boundaries{bdry_index};
    n = bdry.node_count;
    x_b = bdry.boundary_nodes(:,1);
    z_b = bdry.boundary_nodes(:,2);
    K_n = zeros(N,n); K_d = zeros(N,n);
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

    %not sure why failing to set these in this way throws an error,
    %but here we are.
    splx = zeros(n-1,1); splz = zeros(n-1,1);
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
    %each segment: integ over [elem, elem+1]
    %for elem = 1:(n-1)
    parfor (elem = 1:(n-1), workers)
        %constants for quadrature
        QUADRATURE_PTS = 16;
        LOG_QUAD_T = [0.0038978344871159159405749 0.0230289456168732385721309 0.0582803983062404121207045 0.1086783650910540383049963 0.1726094549098439456802367 0.2479370544705784829009332 0.3320945491299171492549647 0.4221839105819485959969484 0.5150824733814626243955104 0.6075561204477287757796944 0.6963756532282140421230565 0.7784325658732653696603165 0.8508502697153911276117810 0.9110868572222718952957621 0.9570255717035421882954438 0.9870478002479844414907006]';
        LOG_QUAD_W = [0.0607917100435912335920641 0.1029156775175821408874199 0.1223556620460092003721542 0.1275692469370159898289785 0.1230135746000709101588555 0.1118472448554857223701475 0.0965963851521243477282752 0.0793566643514731356878755 0.0618504945819652041105741 0.0454352465077266717830007 0.0310989747515818051870617 0.0194597659273608412922041 0.0107762549632055264908770 0.0049725428900876415469479 0.0016782011100511945445035 0.0002823537646684363164144]';
    
        QUAD_T = [0.0052995325041750333469603 0.0277124884633837102743126 0.0671843988060841224019271 0.1222977958224984867952045 0.1910618777986781147149031 0.2709916111713863151599924 0.3591982246103705422868302 0.4524937450811812866824368 0.5475062549188187688287144 0.6408017753896294577131698 0.7290083888286137403511589 0.8089381222013218852850969 0.8777022041775015548381589 0.9328156011939158220869217 0.9722875115366163001340283 0.9947004674958249692551249]';
        QUAD_W = [0.0135762297058770482066636 0.0311267619693239468159351 0.0475792558412463928441127 0.0623144856277669384470030 0.0747979944082883679845608 0.0845782596975012679330064 0.0913017075224617918882686 0.0947253052275342510846201 0.0947253052275342510846201 0.0913017075224617918882686 0.0845782596975012679330064 0.0747979944082883679845608 0.0623144856277669384470030 0.0475792558412463928441127 0.0311267619693239468159351 0.0135762297058770482066636]';
        
        %matlab forces reduction ops on arrays to be on the whole array,
        %unless there is syntax I don't know of. Our red-op is addition.
        K_n_accum = zeros(N,n);
        K_d_accum = zeros(N,n);

        %matlab is screaming about temp vars, so force set here.
        shape_start = 0; shape_end = 0; shape0 = 0; shape1 = 0;
        shape_quad = 0; shapep_quad = 0; shape_logquad = 0;
        shape_oneminuslogquad = 0; ax = 0; az = 0; seg_nodes = 0;
%==========================================================================
        %get shape functions
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
            shape_end = shape_start + M;

            offset = (elem - shape_start)/M;
            %integrate between [offset, offset + 1/M]
            %(node elem to node elem+1, where the shape starts at elem-off)

            shape0 = ((offset).^(0:M) * Lmat');
            shape1 = ((offset+1/M).^(0:M) * Lmat');
            shape_quad = ((offset + QUAD_T/M).^(0:M) * Lmat');
            shapep_quad = ((offset + QUAD_T/M).^(0:M-1) * Lpmat');
            shape_logquad = ((offset + LOG_QUAD_T/M)...
                .^(0:M) * Lmat');
            shape_oneminuslogquad = ((offset+(1-LOG_QUAD_T)/M)...
                .^(0:M) * Lmat');            
        else
            error("Boundary %d: Invalid function interpolation rule!",...
                bdry_index);
        end

%==========================================================================
        %interpolation segment
        if strcmp(bdry_sf_bdry_interp_rule,'spline')
            %cubic coefs, higher order first
            ax = splx(elem,:); az = splz(elem,:);
            %these are already interp on interval [0,1]
            r0 = [ax(4), az(4)];
            rp0 = [ax(3), az(3)];
            r1 = [sum(ax), sum(az)];
            rp1 = [dot(ax,[3,2,1,0]),dot(az,[3,2,1,0])];

            r_quad = (QUAD_T.^(3:-1:0)) * [ax',az'];
            rp_quad = ((3:-1:1).*QUAD_T.^(2:-1:0)) * [ax(1:3)',az(1:3)'];

%             r_logquad = (LOG_QUAD_T.^(3:-1:0)) * [ax',az'];
            rp_logquad = ((3:-1:1).*LOG_QUAD_T.^(2:-1:0))...
                * [ax(1:3)',az(1:3)'];

%             r_oneminuslogquad = ((1-LOG_QUAD_T).^(3:-1:0)) * [ax',az'];
            rp_oneminuslogquad = ((3:-1:1).*(1-LOG_QUAD_T).^(2:-1:0))...
                * [ax(1:3)',az(1:3)'];
        else
            seg_nodes = squeeze(seg_params(elem,:,:));
            % parameterize curve: 
            % [0,1] -> interp boundary(elem:elem+1)
            offset = elem - seg_starts(elem);
            r0 = seg_nodes(offset+1,:);
            r1 = seg_nodes(offset+2,:);
            rp0=((offset/M_bdry).^(0:M_bdry-1)*Lpmat_bdry')...
                *seg_nodes;
            rp1=(((offset+1)/M_bdry).^(0:M_bdry-1)*Lpmat_bdry')...
                *seg_nodes;

            QUAD_REPARAM = (offset + QUAD_T)./M_bdry;
            r_quad = (QUAD_REPARAM.^(0:M_bdry) * Lmat_bdry') * seg_nodes;
            rp_quad = (QUAD_REPARAM.^(0:M_bdry-1) * Lpmat_bdry')*seg_nodes;

            LOG_QUAD_REPARAM = (offset + LOG_QUAD_T)./M_bdry;
%             r_logquad = (LOG_QUAD_REPARAM.^(0:M_bdry) * Lmat')...
%                 * seg_nodes;
            rp_logquad = (LOG_QUAD_REPARAM.^(0:M_bdry-1) * Lpmat_bdry')...
                * seg_nodes;

            ONEMINUS_LOG_QUAD_REPARAM = (offset + (1-LOG_QUAD_T))/M_bdry;
%             r_oneminuslogquad = (ONEMINUS_LOG_QUAD_REPARAM.^(0:M_bdry)...
%                 * Lmat') * seg_nodes;
            rp_oneminuslogquad =(ONEMINUS_LOG_QUAD_REPARAM.^(0:M_bdry-1)...
                * Lpmat_bdry') * seg_nodes;
        end

%==========================================================================
        %speed of curve: used in integral
        dgam_quad = vecnorm(rp_quad,2,2);
        %different shape used for vectorized accumulation
        dgam_quad_ = zeros(1,1,QUADRATURE_PTS);
        dgam_quad_(1,1,:) = dgam_quad;

        %different shape of r_quad as well:
        r_quad_ = zeros(1,2,QUADRATURE_PTS);
        r_quad_(1,:,:) = r_quad';
        
        diff = r1 - r0;
        length2 = sum(diff.^2);
        invlength2 = 1/length2;
        invlength = sqrt(invlength2);
        %project onto segment and see how far down the result is
        % seg_frac in [0,1] means that its on the segment
        seg_frac = (singularities - r0) * (diff' .* invlength2);

        %which segments is a simple quadrature valid? we do this by
        %doing capsule collision (are we within a distance eps of a line
        %segment?
        far_away=sum((singularities - r0).^2,2) > quasisingular_thresh...
            ...             ^out of circle around r0
            & sum((singularities - r1).^2,2) > quasisingular_thresh...
            ...             ^and out of circle around r1
            & (seg_frac < 0 | seg_frac > 1 ...
            | abs((singularities - r0) * ([invlength -invlength]...
                    .*flip(diff))') > sqrt(quasisingular_thresh));
            % and outside rectangle aligned to line segment with height eps
        
        %autofail singularities:
        seg_sing_ind = sing_indices(bdry_index) + elem;
        seg_sing_indm1 = seg_sing_ind - 1;
        seg_sing_ind = seg_sing_ind - N*(seg_sing_ind > N);
        far_away(seg_sing_indm1) = 0;
        far_away(seg_sing_ind) = 0;
%==========================================================================
        %calculate the simple quadrature integrals
        K_d_accum(far_away,shape_start:shape_end) = (-1/(4*pi))*squeeze(...
         log(sum((r_quad_-singularities(far_away,:)).^2,2)).*dgam_quad_)...
         * (QUAD_W .* shape_quad);
        
        rp_quad_ = zeros(1,2,QUADRATURE_PTS);
        rp_quad_(1,:,:) = rp_quad';
        K_n_accum(far_away,shape_start:shape_end) = (-1/(2*pi))*squeeze(...
            sum((r_quad_-singularities(far_away,:)) .* ...
            (flip(rp_quad_,2).*[1,-1]),2)./...
            (sum((r_quad_-singularities(far_away,:)).^2,2)))...
            * (QUAD_W .* shape_quad);
%==========================================================================
        %calculate integrals with singularities

        %singularity at start (seg_sing_ind-1)
        %v alignment
        % first figure out which transformation to v we should use
        if abs(rp0(1)) >= abs(rp0(2))
            %inclined closer to x-axis: use delta z / delta x
            %v = @(t) atan((dot(az,[t^3 t^2 t 1]) - r0(2))/...
            %    (dot(ax,[t^3 t^2 t 1]) - r0(1)));
            v0 = atan(rp0(2)/rp0(1));
            v1 = atan((r1(2) - r0(2))./(r1(1) - r0(1)));
            v_quad = atan((r_quad(:,2) - r0(2))./(r_quad(:,1) - r0(1)));
        else
            %inclined closer to z-axis: use delta x / delta z
            %v = @(t) -atan((dot(ax,[t^3 t^2 t 1]) - r0(1))/...
            %    (dot(az,[t^3 t^2 t 1]) - r0(2)));
            v0 = -atan(rp0(1)/rp0(2));
            v1 = -atan((r1(1) - r0(1))./(r1(2) - r0(2)));
            v_quad = -atan((r_quad(:,1) - r0(1))./(r_quad(:,2) - r0(2)));
        end
        %K_d
        %integrate 1/(2*pi) N_j(t)J(t) [ -log(|r(t) - r0| / t) - log(t) ]
        K_d_accum(seg_sing_indm1,shape_start:shape_end) =...
            (-1/(2*pi)) .* ((dgam_quad...
            .* log(vecnorm(r_quad - singularities(seg_sing_indm1,:),...
            2,2)./QUAD_T))' * (QUAD_W .* shape_quad))...
            ...multiplication by -log(t) is handled by quadrature
            + (1/(2*pi)) .* (vecnorm(rp_logquad, 2,2)' *...
            (LOG_QUAD_W .* shape_logquad));

        %K_n
        %evaluate -1(2*pi) [ N_j(1)v(1) - N_j(0)v(0) 
        % - int_0^1( v(t)N_j'(t) , dt) ]
        K_n_accum(seg_sing_indm1,shape_start:shape_end) =...
            (-1/(2*pi)) .* (shape1.*v1 - shape0.*v0...
            - (v_quad)' * (QUAD_W .* shapep_quad));
        
        %singularity at end (seg_sing_ind)
        % first figure out which transformation to v we should use
        if abs(rp0(1)) >= abs(rp0(2))
            %inclined closer to x-axis: use delta z / delta x
            v0 = atan((r0(2) - r1(2))./(r0(1) - r1(1)));
            v1 = atan(rp1(2)/rp1(1));
            v_quad = atan((r_quad(:,2) - r1(2))./(r_quad(:,1) - r1(1)));
        else
            %inclined closer to z-axis: use delta x / delta z
            v0 = -atan((r0(1) - r1(1))./(r0(2) - r1(2)));
            v1 = -atan(rp1(1)/rp1(2));
            v_quad = -atan((r_quad(:,1) - r1(1))./(r_quad(:,2) - r1(2)));
        end
        
        %K_d
        %integrate 1/(2*pi) N_j(t)J(t) [ -log(|r(t) - r1| /(1-t))-log(1-t)]
        K_d_accum(seg_sing_ind,shape_start:shape_end) =...
            (-1/(2*pi)) .* ((dgam_quad...
            .* log(vecnorm(r_quad - singularities(seg_sing_ind,:),...
            2,2)./(1-QUAD_T)))' * (QUAD_W .* shape_quad))...
            ...multiplication by -log(t) is handled by quadrature
            + (1/(2*pi)) .* (vecnorm(rp_oneminuslogquad, 2,2)' *...
            (LOG_QUAD_W .* shape_oneminuslogquad));

        %K_n
        %evaluate -1(2*pi) [ N_j(1)v(1) - N_j(0)v(0) 
        % - int_0^1( v(t)N_j'(t) , dt)]
        K_n_accum(seg_sing_ind,shape_start:shape_end) =...
            (-1/(2*pi)) .* (shape1.*v1 - shape0.*v0...
            - (v_quad)' * (QUAD_W .* shapep_quad));
        
%==========================================================================
        %loop over other nodes; quasisingular
        %don't loop over the singularity nodes; mask them out.
        far_away(seg_sing_indm1) = 1;
        far_away(seg_sing_ind) = 1;
        close = find(~far_away);
        if ~isempty(close)
        for j = close'
            if strcmp(bdry_sf_bdry_interp_rule,'spline')
                [In,Id] = segment_In_Id_integrate(ax,az,0,1,...
                        singularities(j,:),elem-shape_start,0,0,0,M,...
                        Lmat,Lpmat,cos_alpha_max, subint_min,1);
            else
                bdry_off = elem - seg_starts(elem);
                [In,Id] = segment_In_Id_integrate( ...
                        seg_nodes(:,1)', seg_nodes(:,2)',...
                        bdry_off/M_bdry,(bdry_off+1)/M_bdry,...
                        singularities(j,:),seg_starts(elem)-shape_start,...
                        M_bdry,Lmat_bdry,Lpmat_bdry,M,Lmat,Lpmat,...
                        cos_alpha_max, subint_min,0);
            end
            K_d_accum(j,shape_start:shape_end) = Id;
            K_n_accum(j,shape_start:shape_end) = In;
        end
        end

        K_n = K_n + K_n_accum;
        K_d = K_d + K_d_accum;
    end
    
    bdry.characteristics.K_n = K_n;
    bdry.characteristics.K_d = K_d;
    
end

end