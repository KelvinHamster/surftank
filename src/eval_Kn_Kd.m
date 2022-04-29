function [obj, I_d_FS] = eval_Kn_Kd(obj)
%Performs the first step of the BEM, evaluating shape function integrals on the boundary and storing them in temp_.
%   Takes an optional argument workers, which determines the number of workers.
%   parallel for-loop to process the integrals.

workers = obj.meta.sim_parallel_workers;

%an interval is treated as quasisingular if one if the nodes is within
%sqrt(quasisingular_thresh) distance (in 2-norm) from the singularity
quasisingular_thresh = 5*10^-2;
cos_alpha_max = 0.7660;
subint_min = 0.0625; %minimum length for subdivisions

%required sim values

x = [obj.boundary.x obj.boundary.x(1)];
z = [obj.boundary.z obj.boundary.x(1)];
doublenode = obj.boundary.doublenode;

FS_i = obj.boundary.FS_start;
FS_f = obj.boundary.FS_end;

N = obj.boundary.N;

%characterize the interpolation scheme
M = obj.interpolation.M;
Lmat = obj.interpolation.Lagrange;
Lpmat = obj.interpolation.Lagrange_prime;
mnode = linspace(0,1,M+1);

do_MCI = obj.meta.integral_spline;

%TODO handle spline corner conditions
if FS_i<FS_f
    fs_x = csape(FS_i:FS_f, x(FS_i:FS_f)).coefs;
    fs_z = csape(FS_i:FS_f, z(FS_i:FS_f)).coefs;
    obj.characteristics_.FS_spline_z = fs_z;
    obj.characteristics_.FS_spline_x = fs_x;
end
M_FS = obj.interpolation_FS.M;
Lmat_FS = obj.interpolation_FS.Lagrange;
Lpmat_FS = obj.interpolation_FS.Lagrange_prime;
%tolerance for when to use L'Hospital's rule in the v-calculation
EPS = 10^-6;

FS_size = (FS_f-FS_i);

I_d_FS = zeros(FS_size,N,M_FS+1);%indexed by element ID, singularity location, index on element
I_n_FS = zeros(FS_size,N,M_FS+1);%indexed by element ID, singularity location, index on element
sliding_segs = zeros(FS_size,1);%t=0 for the sliding element at each segment

i_itermax = M_FS+1; %matlab does not like upper limits not like this

%for fs_elem = 1:FS_size%k = FS_i:(FS_f-1)
parfor (fs_elem = 1:FS_size, workers)%k = FS_i:(FS_f-1)
    k = fs_elem + FS_i - 1;
    sliding_seg = max(FS_i, min(FS_f - M_FS, k - floor(M_FS/2)));
    shape_offset = k - sliding_seg;
    sliding_segs(fs_elem) = sliding_seg;
    %   [--element--]
    % - o - o - o - o - o -
    %   ^   k  k+1
    % sliding_seg
    %
    %   ^---^
    % shape_offset
    %
    % elements: j = sliding_seg:sliding_seg+M_FS contribute with shape
    % function L_{j-sliding_seg+1}(t + shape_offset)
    
    %16 node integration, t-nodes, w-weights
    QUADRATURE_PTS = 16;
    LOG_QUAD_T = [0.0038978344871159159405749 0.0230289456168732385721309 0.0582803983062404121207045 0.1086783650910540383049963 0.1726094549098439456802367 0.2479370544705784829009332 0.3320945491299171492549647 0.4221839105819485959969484 0.5150824733814626243955104 0.6075561204477287757796944 0.6963756532282140421230565 0.7784325658732653696603165 0.8508502697153911276117810 0.9110868572222718952957621 0.9570255717035421882954438 0.9870478002479844414907006]';
    LOG_QUAD_W = [0.0607917100435912335920641 0.1029156775175821408874199 0.1223556620460092003721542 0.1275692469370159898289785 0.1230135746000709101588555 0.1118472448554857223701475 0.0965963851521243477282752 0.0793566643514731356878755 0.0618504945819652041105741 0.0454352465077266717830007 0.0310989747515818051870617 0.0194597659273608412922041 0.0107762549632055264908770 0.0049725428900876415469479 0.0016782011100511945445035 0.0002823537646684363164144]';

    QUAD_T = [0.0052995325041750333469603 0.0277124884633837102743126 0.0671843988060841224019271 0.1222977958224984867952045 0.1910618777986781147149031 0.2709916111713863151599924 0.3591982246103705422868302 0.4524937450811812866824368 0.5475062549188187688287144 0.6408017753896294577131698 0.7290083888286137403511589 0.8089381222013218852850969 0.8777022041775015548381589 0.9328156011939158220869217 0.9722875115366163001340283 0.9947004674958249692551249]';
    QUAD_W = [0.0135762297058770482066636 0.0311267619693239468159351 0.0475792558412463928441127 0.0623144856277669384470030 0.0747979944082883679845608 0.0845782596975012679330064 0.0913017075224617918882686 0.0947253052275342510846201 0.0947253052275342510846201 0.0913017075224617918882686 0.0845782596975012679330064 0.0747979944082883679845608 0.0623144856277669384470030 0.0475792558412463928441127 0.0311267619693239468159351 0.0135762297058770482066636]';

    %store values of parameterization r at points of the log quadrature
    r_logquad = zeros(QUADRATURE_PTS,2);
    rp_logquad = zeros(QUADRATURE_PTS,2);
    %store values of parameterization r at points of the log quadrature with
    %(1-t) because this quadrature is not symmetric
    r_oneminuslogquad = zeros(QUADRATURE_PTS,2);
    rp_oneminuslogquad = zeros(QUADRATURE_PTS,2);
    %store values of parameterization r at points of the G-L quadrature, which
    %is symmetric, so we can flip up-down for one-minus-quad
    r_quad = zeros(QUADRATURE_PTS,2);
    rp_quad = zeros(QUADRATURE_PTS,2);
    
    %store values of the shape function
    shape0 = ((shape_offset).^(0:M_FS) * Lmat_FS');
    shape1 = ((shape_offset+1).^(0:M_FS) * Lmat_FS');
    shape_quad = ((QUAD_T+shape_offset).^(0:M_FS) * Lmat_FS');
    shapep_quad = ((QUAD_T+shape_offset).^(0:M_FS-1) * Lpmat_FS');
    shape_logquad = ((LOG_QUAD_T+shape_offset).^(0:M_FS) * Lmat_FS');
    shape_oneminuslogquad = ((shape_offset+1-LOG_QUAD_T).^(0:M_FS) * Lmat_FS');
    
    %free surface coefficients: a4 + a3*xi + a2*xi^2 + a1*xi^3
    ax = fs_x(fs_elem,:);
    az = fs_z(fs_elem,:);
    %r = @(t) [dot(ax,[t^3 t^2 t 1]) dot(az,[t^3 t^2 t 1])];%r
    %rp =@(t) [dot(ax,[3*t^2 2*t 1 0]) dot(az,[3*t^2 2*t 1 0])];%r'
    
    
    
        if do_MCI
    
        %function at known sample points for quadrature expansion
        r0 = [ax(4) az(4)];
        rp0 = [ax(3) az(3)];
        r1 = [sum(ax) sum(az)];
        rp1 = [dot(ax,[3 2 1 0]) dot(az,[3 2 1 0])];

        for i=1:QUADRATURE_PTS
            t = QUAD_T(i);
            r_quad(i,1) = dot(ax,[t^3 t^2 t 1]);
            r_quad(i,2) = dot(az,[t^3 t^2 t 1]);
            rp_quad(i,1) = dot(ax,[3*t^2 2*t 1 0]);
            rp_quad(i,2) = dot(az,[3*t^2 2*t 1 0]);
            t = LOG_QUAD_T(i);
            r_logquad(i,1) = dot(ax,[t^3 t^2 t 1]);
            r_logquad(i,2) = dot(az,[t^3 t^2 t 1]);
            rp_logquad(i,1) = dot(ax,[3*t^2 2*t 1 0]);
            rp_logquad(i,2) = dot(az,[3*t^2 2*t 1 0]);
            t = 1-t;
            r_oneminuslogquad(i,1) = dot(ax,[t^3 t^2 t 1]);
            r_oneminuslogquad(i,2) = dot(az,[t^3 t^2 t 1]);
            rp_oneminuslogquad(i,1) = dot(ax,[3*t^2 2*t 1 0]);
            rp_oneminuslogquad(i,2) = dot(az,[3*t^2 2*t 1 0]);
        end
    
    else
        seg_nodes = [x(sliding_seg:sliding_seg+M_FS)', z(sliding_seg:sliding_seg+M_FS)'];
        
        r0 = [x(k), z(k)]; r1 = [x(k+1), z(k+1)];
        rp0 = ((shape_offset).^(0:M_FS-1) * Lpmat_FS') * seg_nodes;
        rp1 = ((shape_offset+1).^(0:M_FS-1) * Lpmat_FS') * seg_nodes;
        r_quad = shape_quad * seg_nodes;
        rp_quad = shapep_quad * seg_nodes;
        rp_logquad = ((LOG_QUAD_T+shape_offset).^(0:M_FS-1) * Lpmat_FS') * seg_nodes;
        rp_oneminuslogquad = ((shape_offset+1-LOG_QUAD_T).^(0:M_FS-1) * Lpmat_FS') * seg_nodes;
    
    end
    
    norm_rp_quad = vecnorm(rp_quad,2,2);
    normal_times_abs_rp_quad = fliplr(rp_quad).*[1 -1];
    
    %only two nodes with nonzero shape function j=k and j=k+1

    %for wrapping over; equate point N+1 and point 1
    kplus1 = k + 1;
    if kplus1 > N
        kplus1 = kplus1 - N;
    end
    
    diff = r1 - r0;
    length2 = sum(diff.^2);
    invlength2 = 1/length2;
    invlength = sqrt(invlength2);
    
    %scatter([r0(1) r1(1)], [r0(2) r1(2)])
    %plot(QUAD_T,shapep_quad(:,1))
    %hold on;
    
    for l = 1:N
        
        r_sing = [x(l) z(l)];%singularity point at l
        
        %project onto segment and see how far down the result is
        % seg_frac in [0,1] means that its on the segment
        seg_frac = dot(diff, r_sing - r0)/length2;
        

        if k == l || (k==(l+1 - N*(l==N)) && doublenode(l))
            %contains a singularity at the start of this segment
            
            %v alignment
            % first figure out which transformation to v we should use
            if abs(rp0(1)) >= abs(rp0(2))
                %inclined closer to x-axis: use delta z / delta x
                %v = @(t) atan((dot(az,[t^3 t^2 t 1]) - r0(2))/(dot(ax,[t^3 t^2 t 1]) - r0(1)));
                v0 = atan(rp0(2)/rp0(1));
                v1 = atan((r1(2) - r0(2))./(r1(1) - r0(1)));
                v_quad = atan((r_quad(:,2) - r0(2))./(r_quad(:,1) - r0(1)));
            else
                %inclined closer to z-axis: use delta x / delta z
                %v = @(t) -atan((dot(ax,[t^3 t^2 t 1]) - r0(1))/(dot(az,[t^3 t^2 t 1]) - r0(2)));
                v0 = -atan(rp0(1)/rp0(2));
                v1 = -atan((r1(1) - r0(1))./(r1(2) - r0(2)));
                v_quad = -atan((r_quad(:,1) - r0(1))./(r_quad(:,2) - r0(2)));
            end
            
            
            for i=1:i_itermax
                %K_d
                %integrate 1/(2*pi) N_j(t)J(t) [ -log(|r(t) - r0| / t) - log(t) ]
                I_j = -dot(shape_quad(:,i) .* norm_rp_quad .* log(vecnorm(r_quad - r_sing, 2,2)./QUAD_T) , QUAD_W);
                % multiplication by -log(t) is handled by quadrature
                I_j = I_j + dot(shape_logquad(:,i) .* vecnorm(rp_logquad, 2,2), LOG_QUAD_W);
                I_d_FS(fs_elem,l,i) = 1/(2*pi) * I_j;
                
                %K_n
                %evaluate -1(2*pi) [ N_j(1)v(1) - N_j(0)v(0) - int_0^1( v(t)N_j'(t) , dt) ]
                I_n_FS(fs_elem,l,i) = -1/(2*pi) * (shape1(i).*v1 - shape0(i).*v0 - dot( shapep_quad(:,i).*v_quad , QUAD_W));
                
            end

        elseif kplus1 == l || (kplus1==(l-1 + N*(l==1)) && doublenode(kplus1))
            %contains a singularity at the end of the segment; our
            %log_integral_solve code only works when the singularity is at
            %t=0, so swap the bounds for log integrals.
            
            
            % first figure out which transformation to v we should use
            if abs(rp0(1)) >= abs(rp0(2))
                %inclined closer to x-axis: use delta z / delta x
                %v = @(t) atan((dot(az,[t^3 t^2 t 1]) - r1(2))/(dot(ax,[t^3 t^2 t 1]) - r1(1)));
                %v_int = dot(atan((r_quad(:,2) - r1(2))./(r_quad(:,1) - r1(1))), QUAD_W);
                v0 = atan((r0(2) - r1(2))./(r0(1) - r1(1)));
                v1 = atan(rp1(2)/rp1(1));
                v_quad = atan((r_quad(:,2) - r1(2))./(r_quad(:,1) - r1(1)));
            else
                %inclined closer to z-axis: use delta x / delta z
                %v = @(t) -atan((dot(ax,[t^3 t^2 t 1]) - r1(1))/(dot(az,[t^3 t^2 t 1]) - r1(2)));
                %v_int = -dot(atan((r_quad(:,1) - r1(1))./(r_quad(:,2) - r1(2))), QUAD_W);
                v0 = -atan((r0(1) - r1(1))./(r0(2) - r1(2)));
                v1 = -atan(rp1(1)/rp1(2));
                v_quad = -atan((r_quad(:,1) - r1(1))./(r_quad(:,2) - r1(2)));
            end
            
            for i=1:i_itermax
                %K_d
                %integrate 1/(2*pi) N_j(t)J(t) [ -log(|r(t) - r1| /(1-t)) - log(1-t) ]
                I_j = -dot(shape_quad(:,i) .* norm_rp_quad .* log(vecnorm(r_quad - r_sing, 2,2)./(1-QUAD_T)) , QUAD_W);
                % multiplication by -log(t') is handled by quadrature (t' = 1-t)
                I_j = I_j + dot(shape_oneminuslogquad(:,i) .* vecnorm(rp_oneminuslogquad, 2,2), LOG_QUAD_W);
                I_d_FS(fs_elem,l,i) = 1/(2*pi) * I_j;
                
                %K_n
                %evaluate -1(2*pi) [ N_j(1)v(1) - N_j(0)v(0) - int_0^1( v(t)N_j'(t) , dt) ]
                I_n_FS(fs_elem,l,i) = -1/(2*pi) * (shape1(i).*v1 - shape0(i).*v0 - dot( shapep_quad(:,i).*v_quad , QUAD_W));
                
            end
            
        elseif sum((r_sing - r0).^2) < quasisingular_thresh || sum((r_sing - r1).^2) < quasisingular_thresh || (0 <= seg_frac && seg_frac <= 1 && (dot([1 -1].*flip(diff), r_sing - r0)*invlength)^2 < quasisingular_thresh)
            %no singularities but we are close to one
            
            if do_MCI
                [In,Id] = segment_In_Id_integrate_FS(ax,az,0,1,r_sing,shape_offset,M_FS,Lmat_FS,Lpmat_FS,...
                        cos_alpha_max, subint_min,1);
            else
                [In,Id] = segment_In_Id_integrate_FS(x(sliding_seg:sliding_seg+M_FS),...
                        z(sliding_seg:sliding_seg+M_FS),...
                        0,1,r_sing,shape_offset,M_FS,Lmat_FS,Lpmat_FS,cos_alpha_max, subint_min,0);
            end
                    
                    
            for i=1:i_itermax
                I_n_FS(fs_elem,l,i) = In(i);
                I_d_FS(fs_elem,l,i) = Id(i);
                
                %I_d_FS(fs_elem,l,i) = -1/(4*pi) * dot( shape_quad(:,i).*log(sum((r_quad - r_sing).^2, 2)).*norm_rp_quad , QUAD_W);
                
                %K_n
                %evaluate -1(2*pi) N_j(t)J(t) <r(t), n(t)>/|r(t)-r_sing|^2
                %TODO better integration scheme
                %I_n_FS_start(fs_elem,l) = -1/(2*pi) * dot( (1-QUAD_T).*(dot(r_quad - r_sing, normal,2)./sum((r_quad - r_sing).^2, 2)).*norm_rp_quad , QUAD_W);
                %I_n_FS_end(fs_elem,l) = -1/(2*pi) * dot( (QUAD_T).*(dot(r_quad - r_sing, normal, 2)./sum((r_quad - r_sing).^2, 2)).*norm_rp_quad , QUAD_W);
                %r = @(t) [dot(ax,[t^3 t^2 t 1]) dot(az,[t^3 t^2 t 1])];%r
                %rp =@(t) [dot(ax,[3*t^2 2*t 1 0]) dot(az,[3*t^2 2*t 1 0])];%r'
                %norm_func =@(t) flip(rp(t)).*[1 -1];
                %shape_func = @(t) (t+shape_offset).^(0:M_FS) * Lmat_FS(i,:)';
                %I_n_FS(fs_elem,l,i) = -1/(2*pi) * integral(@(t) shape_func(t) .* dot(r(t)-r_sing,norm_func(t))./sum((r(t) - r_sing).^2),0,1,'ArrayValued',1,'AbsTol',INTEGRAL_TOL);
                
            end
            
            
        else
            %contains no singularities; we can solve with basic quadrature
            %normal =@(t) flip(rp(t)).*[1 -1]./norm(rp(t));%n
            
            for i=1:i_itermax
                %K_d
                %integrate -1/(2*pi) N_j(t)J(t) log(|r(t) - r_sing|)
                I_d_FS(fs_elem,l,i) = -1/(4*pi) * dot( shape_quad(:,i).*log(sum((r_quad - r_sing).^2, 2)).*norm_rp_quad , QUAD_W);
                
                %K_n
                %evaluate -1(2*pi) N_j(t)J(t) <r(t), n(t)>/|r(t)-r_sing|^2
                I_n_FS(fs_elem,l,i) = -1/(2*pi) * dot( shape_quad(:,i).*(dot(r_quad - r_sing, normal_times_abs_rp_quad, 2)./sum((r_quad - r_sing).^2, 2)) , QUAD_W);
                
            end
        end
    end
end

%non-free surface; use isoparametric elements


valid_K = zeros(1,N);
k = 1;
total_k = 0;
while k <= N-M+1
    if doublenode(k)
        k = k+1;
        continue;
    end
    if k >= FS_i && k < FS_f
        k = FS_f;
        continue
    end
    total_k = total_k + 1;
    valid_K(total_k) = k;
    k = k + M;
end
valid_K = valid_K(1:total_k);


I_d_iso = zeros(total_k,N,M+1);
I_n_iso = zeros(total_k,N,M+1);
i_itermax = M+1; %matlab does not like upper limits not like this

%for iso_elem = 1:total_k
parfor (iso_elem = 1:total_k,workers)
    k = valid_K(iso_elem);
    
    
    %16 node integration, t-nodes, w-weights
    QUADRATURE_PTS = 16;
    LOG_QUAD_T = [0.0038978344871159159405749 0.0230289456168732385721309 0.0582803983062404121207045 0.1086783650910540383049963 0.1726094549098439456802367 0.2479370544705784829009332 0.3320945491299171492549647 0.4221839105819485959969484 0.5150824733814626243955104 0.6075561204477287757796944 0.6963756532282140421230565 0.7784325658732653696603165 0.8508502697153911276117810 0.9110868572222718952957621 0.9570255717035421882954438 0.9870478002479844414907006]';
    LOG_QUAD_W = [0.0607917100435912335920641 0.1029156775175821408874199 0.1223556620460092003721542 0.1275692469370159898289785 0.1230135746000709101588555 0.1118472448554857223701475 0.0965963851521243477282752 0.0793566643514731356878755 0.0618504945819652041105741 0.0454352465077266717830007 0.0310989747515818051870617 0.0194597659273608412922041 0.0107762549632055264908770 0.0049725428900876415469479 0.0016782011100511945445035 0.0002823537646684363164144]';

    QUAD_T = [0.0052995325041750333469603 0.0277124884633837102743126 0.0671843988060841224019271 0.1222977958224984867952045 0.1910618777986781147149031 0.2709916111713863151599924 0.3591982246103705422868302 0.4524937450811812866824368 0.5475062549188187688287144 0.6408017753896294577131698 0.7290083888286137403511589 0.8089381222013218852850969 0.8777022041775015548381589 0.9328156011939158220869217 0.9722875115366163001340283 0.9947004674958249692551249]';
    QUAD_W = [0.0135762297058770482066636 0.0311267619693239468159351 0.0475792558412463928441127 0.0623144856277669384470030 0.0747979944082883679845608 0.0845782596975012679330064 0.0913017075224617918882686 0.0947253052275342510846201 0.0947253052275342510846201 0.0913017075224617918882686 0.0845782596975012679330064 0.0747979944082883679845608 0.0623144856277669384470030 0.0475792558412463928441127 0.0311267619693239468159351 0.0135762297058770482066636]';


    %Interpolate with lagrange: r = sum_{i=1}^{M+1} Li * r(k_i)
    x_elem = x(k:k+M)';
    z_elem = z(k:k+M)';
    %r = @(t) L(1:(M+1),t) * [x_elem z_elem];%r
    %rp =@(t) Lp(1:(M+1),t) * [x_elem z_elem];%r'
    
    
    %function at known sample points for quadrature expansion
    r0 = [x_elem(1) z_elem(1)];
    %rp0 = Lp(1:(M+1),0)*[x_elem z_elem];
    %rp0 = (Lpmat(:,1)')*[x_elem z_elem];
    r1 = [x_elem(M+1) z_elem(M+1)];
    %rp1 = Lp(1:(M+1),1)*[x_elem z_elem];
    %rp1 = sum(Lpmat,2)'*[x_elem z_elem];
    %r_quad = L(1:(M+1),QUAD_T)*[x_elem z_elem];
    r_quad = (QUAD_T.^(0:M) * Lmat')*[x_elem z_elem];
    %rp_quad = Lp(1:(M+1),QUAD_T)*[x_elem z_elem];
    rp_quad = (QUAD_T.^(0:(M-1)) * Lpmat')*[x_elem z_elem];
    norm_rp_quad = vecnorm(rp_quad, 2,2);
    %r_logquad = L(1:(M+1),LOG_QUAD_T)*[x_elem z_elem];
    %r_logquad = (LOG_QUAD_T.^(0:M) * Lmat')*[x_elem z_elem];
    %rp_logquad = Lp(1:(M+1),LOG_QUAD_T)*[x_elem z_elem];
    %rp_logquad = (LOG_QUAD_T.^(0:(M-1)) * Lpmat')*[x_elem z_elem];
    %r_oneminuslogquad = L(1:(M+1),1-LOG_QUAD_T)*[x_elem z_elem];
    %r_oneminuslogquad = ((1-LOG_QUAD_T).^(0:M) * Lmat')*[x_elem z_elem];
    %rp_oneminuslogquad = Lp(1:(M+1),1-LOG_QUAD_T)*[x_elem z_elem];
    %rp_oneminuslogquad = ((1-LOG_QUAD_T).^(0:M-1) * Lpmat')*[x_elem z_elem];
    norm_quad = fliplr(rp_quad).*[1 -1]./norm_rp_quad;%n
    
    
    %for quasisingular checking- approx as line
    diff = r1 - r0;
    length2 = sum(diff.^2);
    invlength2 = 1/length2;
    invlength = sqrt(invlength2);
    
    for l = 1:N
        r_sing = [x(l) z(l)];%singularity point at l
        
        
        
        %project onto segment and see how far down the result is
        % seg_frac in [0,1] means that its on the segment
        seg_frac = dot(diff, r_sing - r0)/length2;

        %the only nodes with nonzero shape function are on this element;
        %also consider wrap-around
        lmk = l - k + N.*(k > l);
        %need to consider double nodes
        if lmk == N-1 && doublenode(l)
            lmk = 0;%l=k and l=k-1 are the same node
        elseif lmk == M+1 && doublenode(l - 1 + N*(l==1))
            lmk = M;%l=k+M and l=k+M+1 are the same node
        end
        if (0 <= lmk && lmk <= M)
            %contains a singularity in this segment
            tl = mnode(lmk + 1); %where on the parameterization the singularity is
            %rp_sing = rp(tl); %derivative at the singularity
            rp_sing = (tl.^(0:(M-1)) * Lpmat')*[x_elem z_elem];
            
            %v substitution at certain values
            v_quad = zeros(QUADRATURE_PTS,1);
            if abs(rp_sing(1)) >= abs(rp_sing(2))
                %inclined closer to x-axis: use delta z / delta x
                v_sing = atan(rp_sing(2)/rp_sing(1));
                if lmk == 0
                    v0 = v_sing;
                    v1 = atan((r1(2) - r_sing(2))/(r1(1) - r_sing(1)));
                elseif lmk == M
                    v0 = atan((r0(2) - r_sing(2))/(r0(1) - r_sing(1)));
                    v1 = v_sing;
                else
                    v0 = atan((r0(2) - r_sing(2))/(r0(1) - r_sing(1)));
                    v1 = atan((r1(2) - r_sing(2))/(r1(1) - r_sing(1)));
                end
                %v = @(t) getfield(struct('a',atan((dot(r(t),[0 1]) - r_sing(2))/(dot(r(t),[1 0]) - r_sing(1))), 'b',v0),char('a' + (abs(dot(r(t),[1 0]) - r_sing(1)) < EPS)));
                %I_j = integral_solve(@(t) v(t)*Lp(i,t));
                for pt = 1:QUADRATURE_PTS
                    t = QUAD_T(pt);
                    if abs(r_quad(pt,1) - r_sing(1)) < EPS
                        v_quad(pt) = v_sing;
                    else
                        v_quad(pt) = atan((r_quad(pt,2) - r_sing(2))/(r_quad(pt,1) - r_sing(1)));
                    end
                end
            else
                %inclined closer to z-axis: use delta x / delta z
                v_sing = -atan(rp_sing(1)/rp_sing(2));
                if lmk == 0
                    v0 = v_sing;
                    v1 = -atan((r1(1) - r_sing(1))/(r1(2) - r_sing(2)));
                elseif lmk == M
                    v0 = -atan((r0(1) - r_sing(1))/(r0(2) - r_sing(2)));
                    v1 = v_sing;
                else
                    v0 = -atan((r0(1) - r_sing(1))/(r0(2) - r_sing(2)));
                    v1 = -atan((r1(1) - r_sing(1))/(r1(2) - r_sing(2)));
                end
                %v = @(t) getfield(struct('a',-atan((dot(r(t),[1 0]) - r_sing(1))/(dot(r(t),[0 1]) - r_sing(2))), 'b', v_sing),char('a' + (abs(dot(r(t),[0 1]) - r_sing(2)) < EPS)));
                %I_j = integral_solve(@(t) v(t)*Lp(i,t));
                for pt = 1:QUADRATURE_PTS
                    t = QUAD_T(pt);
                    if abs(r_quad(pt,2) - r_sing(2)) < EPS
                        v_quad(pt) = v_sing;
                    else
                        v_quad(pt) = -atan((r_quad(pt,1) - r_sing(1))/(r_quad(pt,2) - r_sing(2)));
                    end
                end
            end
            

            for i = 1:i_itermax
                %==============node j; shape function L_i(x)========
                %j = k+i-1;
                %if j > N
                %    j = j - N;
                %end
                % I_{d_jl}^k = -1/2pi int(N_j|r'| * [log(|r - r_sing|/|t-tl|) + log|t-tl|], dt)
                %f = @(t) L(i,t).*norm(rp(t));
                %I_j = -integral_solve(@(t) f(t)*log(norm(r(t) - r_sing)/abs(t-tl)) );
                I_j = -dot( (QUAD_T.^(0:M) * Lmat(i,:)').*norm_rp_quad.*log(vecnorm(r_quad - r_sing,2,2)./abs(QUAD_T-tl)) , QUAD_W);
                if tl > 0
                    %integrate the singular component on [0,tl); variable change so we have integral on (0,1] with singularity on 0
                    %use u= (tl - t)/tl,  du = -1/tl * dt, but we swap the bounds of the integral, so the (-) cancels

                    %note these are equivalent assuming no numerical errors:
                    %I_j = I_j + (1/tl)*integrsolve(@(u) f(tl*(1-u)) *-log(tl*u));
                    %I_j = I_j + (1/tl)*integral_solve(@(u) f(tl*(1-u)) *(-log(tl) -log(u)));
                    %I_j = I_j + tl*(-log(tl)*integral_solve(@(u) f(tl*(1-u))) + log_integral_solve(@(u) f(tl*(1-u))));
                    %I_j = I_j + tl*(-log(tl)*...
                    %    dot( L(i,tl.*(1-QUAD_T)).*vecnorm(Lp(1:(M+1),tl.*(1-QUAD_T))*[x_elem z_elem], 2,2) , QUAD_W) + ...
                    %    dot( L(i,tl.*(1-LOG_QUAD_T)).*vecnorm(Lp(1:(M+1),tl.*(1-LOG_QUAD_T))*[x_elem z_elem], 2,2) , LOG_QUAD_W));
                    I_j = I_j + tl*(-log(tl)*...
                        dot( ((tl.*(1-QUAD_T)).^(0:M) * Lmat(i,:)').*vecnorm(((tl.*(1-QUAD_T)).^(0:(M-1)) * Lpmat')*[x_elem z_elem], 2,2) , QUAD_W) + ...
                        dot( ((tl.*(1-LOG_QUAD_T)).^(0:M) * Lmat(i,:)').*vecnorm(((tl.*(1-LOG_QUAD_T)).^(0:(M-1)) * Lpmat')*[x_elem z_elem], 2,2) , LOG_QUAD_W));
                end
                if tl < 1
                    %integrate the singular component on (tl,1]; variable change so we have integral on (0,1] with singularity on 0
                    %use u= (t-tl)/(1 - tl),  du = 1/(1-tl) * dt
                    %I_j = I_j + (1/(1 - tl))*integral_solve(@(u) f(tl + (1-tl)*u) * -log((1-tl)*u));
                    %I_j = I_j + (1 - tl)*(-log(1-tl)*integral_solve(@(u) f(tl + (1-tl)*u)) + log_integral_solve(@(u) f(tl + (1-tl)*u)) );
                    %I_j = I_j + (1 - tl)*(-log(1-tl)*...
                    %    dot( L(i,tl + (1-tl).*QUAD_T).*vecnorm(Lp(1:(M+1),tl + (1-tl).*QUAD_T)*[x_elem z_elem], 2,2) , QUAD_W) + ...
                    %    dot( L(i,tl + (1-tl).*LOG_QUAD_T).*vecnorm(Lp(1:(M+1),tl + (1-tl).*LOG_QUAD_T)*[x_elem z_elem], 2,2) , LOG_QUAD_W) );
                    I_j = I_j + (1 - tl)*(-log(1-tl)*...
                        dot( ((tl + (1-tl).*QUAD_T).^(0:M) * Lmat(i,:)').*vecnorm(((tl + (1-tl).*QUAD_T).^(0:(M-1)) * Lpmat')*[x_elem z_elem], 2,2) , QUAD_W) + ...
                        dot( ((tl + (1-tl).*LOG_QUAD_T).^(0:M) * Lmat(i,:)').*vecnorm(((tl + (1-tl).*LOG_QUAD_T).^(0:(M-1)) * Lpmat')*[x_elem z_elem], 2,2) , LOG_QUAD_W) );
                end
                %K_d(l,j) = K_d(l,j) + 1/(2*pi) * I_j;
                I_d_iso(iso_elem,l,i) = 1/(2*pi) * I_j;

                %K_n term
                %I_j = integral_solve(@(t) v(t)*Lp(i,t));
                %I_j = dot( v_quad.*Lp(i,QUAD_T) , QUAD_W);
                %K_n(l,j) = K_n(l,j) - 1/(2*pi) * (L(i,1)*v1 - L(i,0)*v0 - I_j);
                
                
                %K_n(l,j) = K_n(l,j) - 1/(2*pi) * (L(i,1)*v1 - L(i,0)*v0 - dot( v_quad.*Lp(i,QUAD_T) , QUAD_W));
                %K_n(l,j) = K_n(l,j) - 1/(2*pi) * (sum(Lmat(i,:),2)*v1 - Lmat(i,1)*v0 - dot( v_quad.*(QUAD_T.^(0:(M-1)) * Lpmat(i,:)') , QUAD_W));
                I_n_iso(iso_elem,l,i) = - 1/(2*pi) * (sum(Lmat(i,:),2)*v1 - Lmat(i,1)*v0 - dot( v_quad.*(QUAD_T.^(0:(M-1)) * Lpmat(i,:)') , QUAD_W));
            end
        elseif sum((r_sing - r0).^2) < quasisingular_thresh || sum((r_sing - r1).^2) < quasisingular_thresh || (0 <= seg_frac && seg_frac <= 1 && (dot([1 -1].*flip(diff), r_sing - r0)*invlength)^2 < quasisingular_thresh)
            %no singularities but we are close to one
            
            [In,Id] = segment_In_Id_integrate_nonFS(x_elem,z_elem,0,1,r_sing,...
                        M,Lmat,Lpmat,cos_alpha_max, subint_min);
            
            for i=1:i_itermax
                %quasisingular log is not that bad
                
                %K_d
                %integrate -1/(2*pi) N_j(t)J(t) log(|r(t) - r_sing|)
                I_d_iso(iso_elem,l,i) = Id(i);
                
                %K_n
                %evaluate -1(2*pi) N_j(t)J(t) <r(t), n(t)>/|r(t)-r_sing|^2
                %TODO better integration scheme
                %I_n_FS_start(fs_elem,l) = -1/(2*pi) * dot( (1-QUAD_T).*(dot(r_quad - r_sing, normal,2)./sum((r_quad - r_sing).^2, 2)).*norm_rp_quad , QUAD_W);
                %I_n_FS_end(fs_elem,l) = -1/(2*pi) * dot( (QUAD_T).*(dot(r_quad - r_sing, normal, 2)./sum((r_quad - r_sing).^2, 2)).*norm_rp_quad , QUAD_W);
                %r = @(t) (t.^(0:M) * Lmat')*[x_elem z_elem];%r
                %rp =@(t) (t.^(0:M-1) * Lpmat')*[x_elem z_elem];%r'
                %norm_func =@(t) flip(rp(t)).*[1 -1];
                %shape_func = @(t) t.^(0:M) * Lmat(i,:)';
                I_n_iso(iso_elem,l,i) = In(i);
            end
        else
            %contains no singularities; we can solve with basic quadrature
            %normal =@(t) flip(rp(t)).*[1 -1]./norm(rp(t));%n
            for i = 1:i_itermax
                %==============node j; shape function L_i(x)========
                %j = k+i-1;
                %if j > N
                %    j = j - N;
                %end
                %K_d(l,j) = K_d(l,j) -1/(2*pi) * integral_solve( @(t) L(i,t)*log(norm(r(t) - r_sing))*norm(rp(t)) );
                %K_d(l,j) = K_d(l,j) -1/(2*pi) * dot( (QUAD_T.^(0:M) * Lmat(i,:)').*log(vecnorm(r_quad - r_sing,2,2)).*norm_rp_quad , QUAD_W);
                I_d_iso(iso_elem,l,i) = -1/(2*pi) * dot( (QUAD_T.^(0:M) * Lmat(i,:)').*log(vecnorm(r_quad - r_sing,2,2)).*norm_rp_quad , QUAD_W);
                
                %K_n(l,j) = K_n(l,j) -1/(2*pi) * integral_solve( @(t) L(i,t)*(dot(r(t) - r_sing, normal(t))/sum((r(t) - r_sing).^2))*norm(rp(t)) );
                %K_n(l,j) = K_n(l,j) -1/(2*pi) * dot( (QUAD_T.^(0:M) * Lmat(i,:)').*(dot(r_quad - r_sing, norm_quad, 2)./sum((r_quad - r_sing).^2 ,2)).*norm_rp_quad, QUAD_W );
                I_n_iso(iso_elem,l,i) = -1/(2*pi) * dot( (QUAD_T.^(0:M) * Lmat(i,:)').*(dot(r_quad - r_sing, norm_quad, 2)./sum((r_quad - r_sing).^2 ,2)).*norm_rp_quad, QUAD_W );
            end
        end
    end
end

K_d = zeros(N); %K_d(l,j) = K_{d}, singularity l, shape of j
K_n = zeros(N); %K_n(l,j) = K_{n}, singularity l, shape of j


for iso_elem = 1:total_k
    %isoparametric elements
    k = valid_K(iso_elem);
    for i = 1:M+1
        j = k + i - 1;
        j = j - N*(j > N);
        K_d(:,j) = K_d(:,j) + I_d_iso(iso_elem,:,i)';
        K_n(:,j) = K_n(:,j) + I_n_iso(iso_elem,:,i)';
    end
end

for fs_elem = 1:FS_size
    %free surface elements
    for i = 1:M_FS+1
        j = sliding_segs(fs_elem) + i - 1;
        j = j - N*(j > N);
        K_d(:,j) = K_d(:,j) + I_d_FS(fs_elem,:,i)';
        K_n(:,j) = K_n(:,j) + I_n_FS(fs_elem,:,i)';
    end
end


COEFS_n = K_n;

for l = 1:N
    COEFS_n(l,l) = -sum(K_n(l,(1:N) ~= l));
end


obj.characteristics_.K_n = K_n;
obj.characteristics_.K_d = K_d;

obj.characteristics_.COEFS_n = COEFS_n;


%===beta values

cosbeta = zeros(1,N);%cosine of the slope at node n
sinbeta = zeros(1,N);%sine of the slope at node n

for k=1:valid_K
    s = (mnode'.^(0:(M-1)) * Lpmat')*[x(k:k+M)' z(k:k+M)'];
    s = s./vecnorm(s,2,2);
    cosbeta(k:k+M) = s(:,1);
    sinbeta(k:k+M) = s(:,2);
end


M_sl = obj.interpolation_sliding.M;
Lpmat_sl = obj.interpolation_sliding.Lagrange_prime;
for k=FS_i:FS_f
    sliding_seg = max(FS_i, min(FS_f - M_sl, k - round(M_sl/2)));
    shift = k - sliding_seg;
    
    deriv_coefs = ((shift/M_sl).^(0:M_sl-1)) * Lpmat_sl';
    
    xp = dot(x(sliding_seg:sliding_seg+M_sl),deriv_coefs);
    zp = dot(z(sliding_seg:sliding_seg+M_sl),deriv_coefs);
    
    %tangent of polynomial at xi=0; corresponds to first node of segment
    len = sqrt(xp^2 + zp^2);
    cosbeta(k) = xp/len;
    sinbeta(k) = zp/len;
end



obj.characteristics_.cosbeta = cosbeta;
obj.characteristics_.sinbeta = sinbeta;

obj.characteristics_.ctime = obj.stepping.t;


%===solve for u,q: generate A and b matrices


FS_size = (FS_f-FS_i + 1);
phi_FS = obj.boundary.phi_FS;
%----------genearate matrix---------------------
%s-indices (dirichlet condition) are given by find(~bt)
%p-indices (neumann condition) are given by find(bt)
bt = ones(1,N);
bt(FS_i:FS_f) = zeros(1,FS_size);
bv = zeros(1,N);
bv(FS_i:FS_f) = phi_FS;


mat_A = bt(1:N).*COEFS_n - (1-bt(1:N)).*K_d;


for i=find(doublenode)
    ip1 = i + 1 - N*(i >= N);

    %if a node has a neumann condition, use phi_i = phi_{i+1}
    if bt(i)
        mat_A(i,:) = zeros(1,N);
        if bt(ip1)                                     %neumann-neumann
            mat_A(i,i) = 1; mat_A(i,ip1) = -1;
        else                                           %neumann-dirichlet
            mat_A(i,i) = 1;
        end
    else                                               %dirichlet-neumann
            mat_A(ip1,:) = zeros(1,N);
            mat_A(ip1,ip1) = 1;
    end
end

mat_b = obj.gen_mat_b(bt,bv);

%TODO allow creation of A,b matrices for different cases



%[mat_L, mat_U, mat_P] = lu(mat_A);
%newvals = (mat_U\(mat_L\(mat_P * mat_b)))';
newvals = (mat_A\mat_b)';

u = bt(1:N).*newvals + (1-bt(1:N)).*bv(1:N);
q = (1-bt(1:N)).*newvals + bt(1:N).*bv(1:N);

obj.characteristics_.phi = u;
obj.characteristics_.phi_n = q;
obj.characteristics_.mat_A = mat_A;
obj.characteristics_.mat_b = mat_b;

%surface derivatives

%solve boundary problem for phi
%dirichlet boundaries on FS, neumann everywhere else (phi_n = 0)


g = obj.meta.g;
%phi = u, phi_n = q

%switch to sliding segments
M = obj.interpolation_sliding.M;
Lmat = obj.interpolation_sliding.Lagrange;
Lpmat = obj.interpolation_sliding.Lagrange_prime;
Lppmat = obj.interpolation_sliding.Lagrange_primeprime;

L = @(i,x) x.^(0:M) * Lmat(i,:)';
Lp = @(i,x) x.^(0:(M-1)) * Lpmat(i,:)';
Lpp = @(i,x) x.^(0:(M-2)) * Lppmat(i,:)';
mnode = linspace(0,1,M+1);

%calculate phi_s, etc. at each collocation point; we do this with first
%derivative
%initialize important information
phi_s_FS = zeros(1,FS_size);
phi_n_FS = zeros(1,FS_size);
phi_ss_FS = zeros(1,FS_size);
phi_ns_FS = zeros(1,FS_size);
phi_t_FS = zeros(1,FS_size);
phi_tn_FS = zeros(1,FS_size);
phi_ts_FS = zeros(1,FS_size);
cosbeta_FS = obj.characteristics_.cosbeta(FS_i:FS_f);
sinbeta_FS = obj.characteristics_.sinbeta(FS_i:FS_f);
beta_s_FS = zeros(1,FS_size);
dxdt_FS = zeros(1,FS_size);
dzdt_FS = zeros(1,FS_size);
d2xdt2_FS = zeros(1,FS_size);
d2zdt2_FS = zeros(1,FS_size);
dphidt_FS = zeros(1,FS_size);
d2phidt2_FS = zeros(1,FS_size);

vnorm2_FS = zeros(1,FS_size);


for j = 1:FS_size
    %place the sliding element, preferably so j is the midpoint
    k = j - floor(M/2);
    %shift so it lies exclusively on the free surface
    if k < 1
        k = 1;
    elseif k+M > FS_size
        k = FS_size-M;
    end
    
    %phi_n is easy, we solved it already
    phi_n_FS(j) = q(j + FS_i - 1);
    
    t_interp = mnode(j - k + 1);
    k_full = k + FS_i - 1;
    jacobian = norm(Lp(1:(M+1),t_interp) * [x(k_full:k_full+M)' z(k_full:k_full+M)']);%r'
    
    %differentiate the interpolation of phi on s to get phi_s
    phi_s_FS(j) = dot(Lp(1:(M+1),t_interp), phi_FS(k:k+M))./jacobian;
    %twice for phi_ss
    shapederiv =  Lp(1:(M+1),t_interp) * [x(k_full:k_full+M)' z(k_full:k_full+M)'];
    shapederiv2= Lpp(1:(M+1),t_interp) * [x(k_full:k_full+M)' z(k_full:k_full+M)'];
    shapeproj = dot(shapederiv,shapederiv2)./(jacobian.^2);
    phi_ss_FS(j) = dot(Lpp(1:(M+1),t_interp) - Lp(1:(M+1),t_interp).*shapeproj, phi_FS(k:k+M))./(jacobian.^2);
    
    
    %beta_s: calculated using chain rule: beta(t) = atan(z'(t)/x'(t))
    % beta'(t) = atan'(z'(t)/x'(t)) * (z'(t)/x'(t))'
    %    = 1/(1 + (z'(t)/x'(t))^2)  * (x'(t)z''(t) - z'(t)x''(t))/(x'(t))^2
    %using the substitution cosbeta = x'(t)/jacobian, sinbeta = z'(t)/jacobian
    %and converting from t space to s-space
    %to:  beta_s = (cosbeta*z_ss - sinbeta * x_ss)/(jacobian^2);
    %derive it first, though; the signs may be wrong
    beta_s_FS(j) = (shapederiv(1)*shapederiv2(2) - shapederiv(2)*shapederiv2(1))./(jacobian.^3);
    
    %use eulerian condition:
    % phi_t   = -1/2 |del phi|^2  - gz - P/rho   (P = 0)
    % Dphi/Dt = +1/2 |del phi|^2  - gz - P/rho   (P = 0)
    v_norm2 = phi_s_FS(j)^2 + phi_n_FS(j)^2;
    vnorm2_FS(j) = v_norm2;
    phi_t_FS(j)  = -v_norm2/2 - g * z(j+FS_i-1);
    dphidt_FS(j) =  v_norm2/2 - g * z(j+FS_i-1);
    
    %space derivatives simply gradient of phi
    dxdt_FS(j) = phi_s_FS(j)*cosbeta_FS(j) + phi_n_FS(j)*sinbeta_FS(j);
    dzdt_FS(j) = phi_s_FS(j)*sinbeta_FS(j) - phi_n_FS(j)*cosbeta_FS(j);
    
    %we will do one more iteration later to get ns,tn,ts and populate 2nd lagrangian derivatives
    
end

%solve BEM again but for phi_t
bt = ones(1,N);
bt(FS_i:FS_f) = zeros(1,FS_size);
bv = zeros(1,N);
bv(FS_i:FS_f) = phi_t_FS;
mat_A = obj.characteristics_.mat_A;
mat_b = obj.gen_mat_b(bt,bv);


%we only had to call BEM once because A is only based on the geometry,
%so when we need to solve for phi_t, we just modify the mat_b values and
%solve.

%newvals = (mat_U\(mat_L\(mat_P * mat_b)))';
newvals = (mat_A \ mat_b)';
phi_tn_FS(:) = newvals(FS_i:FS_f);
%==========================================================================

%last important information items
for j = 1:FS_size
    %place the sliding element, preferably so j is the midpoint
    k = j - floor(M/2);
    %shift so it lies exclusively on the free surface
    if k < 1
        k = 1;
    elseif k+M > FS_size
        k = FS_size-M;
    end
    
    t_interp = mnode(j - k + 1);
    k_full = k + FS_i - 1;
    jacobian = norm(Lp(1:(M+1),t_interp) * [x(k_full:k_full+M)' z(k_full:k_full+M)']);%r'
    
    %differentiate the interpolation of phi_t on s to get phi_ts
    phi_ts_FS(j) = dot(Lp(1:(M+1),t_interp), phi_t_FS(k:k+M))./jacobian;
    
    %differentiate the interpolation of phi_n on s to get phi_ns
    phi_ns_FS(j) = dot(Lp(1:(M+1),t_interp), phi_n_FS(k:k+M))./jacobian;
    
    %D^2phi/Dt^2 = phi_s*phi_st + phi_n*phi_nt + phi_s(phi_s*phi_ss + phi_n*phi_ns)
    %             -phi_n(phi_n*phi_ss + phi_s*phi_ns) - phi_n*beta_s(phi_s^2 + phi_n^2)
    %             - g(phi_s*sinbeta - phi_n*cosbeta) - (1/rho)DP/Dt
    d2phidt2_FS(j) = phi_s_FS(j)*phi_ts_FS(j) + phi_n_FS(j)*phi_tn_FS(j)  ...
        +phi_s_FS(j)*(phi_s_FS(j)*phi_ss_FS(j) + phi_n_FS(j)*phi_ns_FS(j))...
        -phi_n_FS(j)*(phi_n_FS(j)*phi_ss_FS(j) - phi_s_FS(j)*phi_ns_FS(j))...
        -phi_n_FS(j)*beta_s_FS(j)*(phi_s_FS(j)^2 + phi_n_FS(j)^2)...
        - g*(phi_s_FS(j)*sinbeta_FS(j) - phi_n_FS(j)*cosbeta_FS(j));
    
    %D^2x/Dt^2 = cosbeta(phi_ts + phi_s*phi_ss + phi_n*phi_ns)
    %    + sinbeta(phi_tn + phi_s*phi_ns - phi_n*phi_ss - beta_s(phi_n^2 + phi_s^2))
    d2xdt2_FS(j) = cosbeta_FS(j)*(phi_ts_FS(j) + phi_s_FS(j)*phi_ss_FS(j) + phi_n_FS(j)*phi_ns_FS(j))...
        +sinbeta_FS(j)*(phi_tn_FS(j) + phi_s_FS(j)*phi_ns_FS(j) - phi_n_FS(j)*phi_ss_FS(j) - beta_s_FS(j)*(phi_n_FS(j)^2 + phi_s_FS(j)^2));
    %D^2z/Dt^2 = sinbeta(phi_ts + phi_s*phi_ss + phi_n*phi_ns)
    %    + cos(-phi_tn - phi_s*phi_ns + phi_n*phi_ss + beta_s(phi_n^2 + phi_s^2))
    d2zdt2_FS(j) = sinbeta_FS(j)*(phi_ts_FS(j) + phi_s_FS(j)*phi_ss_FS(j) + phi_n_FS(j)*phi_ns_FS(j))...
        +cosbeta_FS(j)*(-phi_tn_FS(j) - phi_s_FS(j)*phi_ns_FS(j) + phi_n_FS(j)*phi_ss_FS(j) + beta_s_FS(j)*(phi_n_FS(j)^2 + phi_s_FS(j)^2));
end

%update either dt or courant number
%first calculate the smallest dx (naiively, only look at adjacent nodes)
xdiff2 = (obj.boundary.x - circshift(obj.boundary.x,-1)).^2 + (obj.boundary.z - circshift(obj.boundary.z,-1)).^2;
dxmin = sqrt(min(xdiff2(~obj.boundary.doublenode)));
% C = u dt/dx
%max free surface velocity = sqrt(max(vnorm2_FS))
umax = sqrt(max(vnorm2_FS));
if obj.stepping.courant_lock
    dt = obj.stepping.courant_number * dxmin / umax;
    obj.stepping.dt = dt;
else
    obj.stepping.courant_number = umax * obj.stepping.dt / dxmin;
end

obj.characteristics_.umax = umax;
obj.characteristics_.dxmin = dxmin;


obj.characteristics_.phi_s_FS = phi_s_FS;
obj.characteristics_.phi_n_FS = phi_n_FS;
obj.characteristics_.phi_ss_FS = phi_ss_FS;
obj.characteristics_.phi_ns_FS = phi_ns_FS;
obj.characteristics_.phi_t_FS = phi_t_FS;
obj.characteristics_.phi_tn_FS = phi_tn_FS;
obj.characteristics_.phi_ts_FS = phi_ts_FS;
obj.characteristics_.cosbeta_FS = cosbeta_FS;
obj.characteristics_.sinbeta_FS = sinbeta_FS;
obj.characteristics_.beta_s_FS = beta_s_FS;
obj.characteristics_.dxdt_FS = dxdt_FS;
obj.characteristics_.dzdt_FS = dzdt_FS;
obj.characteristics_.d2xdt2_FS = d2xdt2_FS;
obj.characteristics_.d2zdt2_FS = d2zdt2_FS;
obj.characteristics_.dphidt_FS = dphidt_FS;
obj.characteristics_.d2phidt2_FS = d2phidt2_FS;
obj.characteristics_.vnorm2_FS = vnorm2_FS;
end

