function did_change = bdry_handle_regrid(bdry, force)
    % Checks if regridding should be done (or forces it if force=1)
    %

    if ~exist("force","var")
        force = 0;
    end
    
    if strcmp(bdry.regridding.mode,'none')
        did_change = 0;
    elseif strcmp(bdry.regridding.mode,'nodeshift')
        did_change = regrid_nodeshift(bdry,force);
    elseif strcmp(bdry.regridding.mode,'linearize_uniform')
        if strcmp(bdry.formulation.type,'parameterized')
            n = bdry.node_count;
            CW = bdry.boundary_nodes(1,:);
            CCW = bdry.boundary_nodes(end,:);
            bdry.boundary_nodes = zeros(n,2);
            bdry.boundary_nodes(1,:) = CW;
            bdry.boundary_nodes(end,:) = CCW;
            param_CW = bdry.formulation.param_CW;
            param_CCW = bdry.formulation.param_CCW;
            for i=1:(n-2)
                t = param_CW + (param_CCW-param_CW)*i/(n-1);
                bdry.boundary_nodes(i+1,:) = ...
                    bdry.formulation.surface_defn(t);
            end
        else
            bdry.boundary_nodes(:,1)=linspace(bdry.boundary_nodes(1,1),...
                bdry.boundary_nodes(end,1),bdry.node_count);
            bdry.boundary_nodes(:,2)=linspace(bdry.boundary_nodes(1,2),...
                bdry.boundary_nodes(end,2),bdry.node_count);
        end
        did_change = 1;
    else
        error("Boundary regridding mode '%s' not valid!", ...
            bdry.regridding.mode);
    end


    bdry.regridding.did_regrid = did_change;
end

function did_change = regrid_nodeshift(bdry, force)
    % Regridding and checks are handled by the parameters in
    % the regridding struct.
    %   regrid_type:
    %     'none' - no regriding is done. force_regrid does not do
    %           anything.
    %     'shift_by_curve3d_normu' - uses curve_renode_3d() to evenly
    %           distribute by curve length in (x,z,phi * a/Umax)
    %           coordinates, where a is specified by
    %           regrid_param{1}, and Umax is calculated as
    %           the maximum sqrt(phi_s^2 + phi_n^2) on the free
    %           surface. The same number of nodes are used, and
    %           regrid_param{2} specifies the interpolation scheme,
    %           which is passed as the interp_scheme argument into
    %           the curve_3d functions.
    %     'shift_by_curve3d_maxphis' - same as
    %           'shift_by_curve3d_normu', but instead of Umax,
    %           max(|phi_s|) is used.
    %     'shift_by_curve2d_normu' -uses curve_renode_by_integral()
    %           using the integrand sqrt(1 + (|u|*a/umax)^2), which
    %           is similar to shift_by_curve3d, except
    %           |nabla(phi)|^2 is used instead of |phi_s|^2. Like
    %           shift_by_curve3d, regrid_param{1} specifies a.
    valid_regrid_types = {'none','shift_by_curve3d_normu',...
        'shift_by_curve3d_maxphis','shift_by_curve2d_normu'};
    % this value is to prevent divide-by-zero on normalization
    min_normalizer = 1e-6;
    %check valid types
    
    if ~any(strcmp(valid_regrid_types,bdry.regridding.regrid_type))
        fprintf(['regrid_type "%s" is invalid. ' ...
            'Change it in the regridding struct!\n'],...
            bdry.regridding.regrid_type);
        return;
    end
    did_change = 0;
    if strcmp(bdry.regridding.regrid_type,'none')
        return;
    end
    
    if isfield(bdry.regridding,"startnode")
        startnode = bdry.regridding.startnode;
    else
        startnode = 1;
    end

    if isfield(bdry.regridding,"endnode")
        endnode = bdry.regridding.endnode;
    elseif isfield(bdry.regridding,"endpadding")
        endnode = bdry.node_count - bdry.regridding.endpadding + 1;
    else
        endnode = bdry.node_count;
    end
    nodes_to_regrid = endnode - startnode + 1;


    %factors
    skip_rg = 1;
    if isfield(bdry.regridding,"global_threshold")
        global_thresh = bdry.regridding.global_threshold;
        skip_rg = 0;
    else
        global_thresh = inf;
    end
    if isfield(bdry.regridding,"local_threshold")
        local_thresh = bdry.regridding.local_threshold;
        skip_rg = 0;
    else
        local_thresh = inf;
    end

    %if no global or local threshold is specified; don't regrid

    %exit conditions: various threshold types fail
    if strcmp(bdry.regridding.regrid_type,'none') || skip_rg
        return;
    elseif strcmp(bdry.regridding.regrid_type, ...
                'shift_by_curve3d_normu')...
                || strcmp(bdry.regridding.regrid_type, ...
                'shift_by_curve3d_maxphis')
        use_normu = strcmp(bdry.regridding.regrid_type, ...
            'shift_by_curve3d_normu');
        if use_normu
            bdry.parent_sim.update_characteristics();
        end
        num_nodes = bdry.node_count;
        if ~isfield(bdry.regridding,"a")
            error("'nodeshift' regridding mode: type '%s' requires" + ...
                " parameter 'a'. Please set regridding.a to a value.", ...
                bdry.regridding.regrid_type);
        end
        a = bdry.regridding.a;
        x = bdry.boundary_nodes(:,1)';
        z = bdry.boundary_nodes(:,2)';
        phi = bdry.characteristics.phi;
        if use_normu
            normalizer = max(bdry.characteristics.umax,min_normalizer);
        else
            %first calculate max |phi_s|
            normalizer = num_nodes;
            M = bdry.regridding.interp.M;
            Lp = bdry.regridding.interp.Lagrange_prime;
            mnode = linspace(0,1,M+1);
            for j = 1:num_nodes
                
                %find the sliding segment
                k = j - floor(M/2);
                if k < 1
                    k = 1;
                elseif k+M > num_nodes
                    k = num_nodes - M;
                end
                %linear combination for lagrange polys
                coefs = mnode(j - k + 1);
                coefs = coefs.^(0:(M-1)) * Lp';
                phis = abs(dot(coefs, phi(k:k+M))) ...
                    /norm(coefs * [x(k:k+M)' z(k:k+M)']);
                if phis > normalizer
                    normalizer = phis;
                end
            end
        end
        V = [x(startnode:endnode); z(startnode:endnode);...
            phi(startnode:endnode)];
        if ~isfield(bdry.regridding,"interp_type")
            error("'nodeshift' regridding mode: type '%s' requires" + ...
                " parameter 'interp_type'. " + ...
                "Please set regridding.interp_type to a value.", ...
                bdry.regridding.regrid_type);
        end
        seglens = curve_lengths_3d(V.*[1;1;a/normalizer], ...
            bdry.regridding.interp_type);
        total_len = sum(seglens);
        %divide by average length
        seglens = (seglens / total_len) * nodes_to_regrid;

        %are the factors within spec?
        neighbor_ratios = seglens(1:end-1) ./ seglens(2:end);
        globalfactor = max(max(seglens), 1/min(seglens));
        localfactor = max(max(neighbor_ratios), 1/min(neighbor_ratios));
        bdry.regridding.local_check_val = localfactor;
        bdry.regridding.global_check_val = globalfactor;

        if (~force) && localfactor < local_thresh...
                   && globalfactor < global_thresh
            return;
        end

    elseif strcmp(bdry.regridding.regrid_type,'shift_by_curve2d_normu')
        if isfield(bdry.regridding,"buffer_length")
            buffer_length = bdry.regridding.buffer_length;
        else
            buffer_length = 0;
        end
        bdry.parent_sim.update_characteristics();
        num_nodes = bdry.node_count;
        if ~isfield(bdry.regridding,"a")
            error("'nodeshift' regridding mode: type '%s' requires" + ...
                " parameter 'a'. Please set regridding.a to a value.", ...
                bdry.regridding.regrid_type);
        end
        a = bdry.regridding.a;
        V = [bdry.boundary_nodes(startnode:endnode,:)';...
            bdry.characteristics.phi(startnode:endnode)];
        vnorm2=max(bdry.characteristics.phi_s(startnode:endnode).^2 ...
            +bdry.characteristics.phi_n(startnode:endnode).^2,...
            min_normalizer);
        seglens = curve_lengths_by_integral(V(1:2,:),...
            sqrt(1 + vnorm2*...
            (a/max(bdry.characteristics.umax,min_normalizer))^2));
        total_len = sum(seglens);
        %divide by average length
        seglens = (seglens / total_len) * nodes_to_regrid;

        %cut out buffer region for testing
        seg_tests = seglens((1+buffer_length):(end-buffer_length));
        seg_tests = (seg_tests / sum(seg_tests)) * length(seg_tests);

        
        
        %are the factors within spec?
        neighbor_ratios = seglens(1:end-1) ./ seglens(2:end);
        globalfactor = max(max(seglens), 1/min(seglens));
        localfactor = max(max(neighbor_ratios), 1/min(neighbor_ratios));
        bdry.regridding.local_check_val = localfactor;
        bdry.regridding.global_check_val = globalfactor;

        if (~force) && localfactor < local_thresh...
                   && globalfactor < global_thresh
            return;
        end
    end

    %do the regrid
    if strcmp(bdry.regridding.regrid_type,'shift_by_curve3d_normu') ||...
            strcmp(bdry.regridding.regrid_type,'shift_by_curve3d_maxphis')
        Vnew = curve_renode_3d(V,nodes_to_regrid, ...
            bdry.regridding.interp_type,[1;1;a/normalizer]);
        bdry.boundary_nodes(startnode:endnode,:) = Vnew(1:2,:)';
        bdry.characteristics.phi(startnode:endnode) = Vnew(3,:);
    elseif strcmp(bdry.regridding.regrid_type,'shift_by_curve2d_normu')
        if isfield(bdry.regridding,"buffer_length")
            buffer_length = bdry.regridding.buffer_length;
        else
            buffer_length = 0;
        end
        
        integrand = sqrt(1 + vnorm2*...
            (a/max(bdry.characteristics.umax,min_normalizer))^2);
        %reweigh: integrand = f(t)  -> f'(t), so that
        %int_0^{start_ref_len} f ds = total_len/nodes_to_regrid
        %int_0^{end_ref_len} f ds = total_len/nodes_to_regrid
        %
        %reweighing changes total_len, though, so estimate that change
        %and build a system of equations to set the LHS and RHS:
        %total_len change estimate is based on a linear interpolation
        %changes by a triangle on each end.
        % len_new - len_old = len1*(f_newstart - f_oldstart)*buffer/2
        %                    +lenN*(f_newend - f_oldend)*buffer/2
        %
        if buffer_length > 0
            sysA = zeros(2); sysb = zeros(2,1);
            %lengths of first and last parts of curve
            % but also multiply by buffer_len / 2 since this is what's used
            len1 = sum(vecnorm(diff(...
                    bdry.boundary_nodes(startnode:(startnode+...
                    buffer_length),:)),2,2))...
                    /2;
            lenN = sum(vecnorm(diff(...
                    bdry.boundary_nodes((endnode-buffer_length):endnode...
                    ,:)),2,2))...
                    /2;
            if startnode == 1
                %force f_new = f_old on start side
                sysA(1,1) = 1; sysA(1,2) = 0;
                sysb(1) = integrand(end);
            else
                %int_0^{end_ref_len} f ds = len_new/nodes_to_regrid
                target_len = vecnorm(diff(...
                    bdry.boundary_nodes((startnode-1):startnode,:)),2);
                sysA(1,1) = target_len*(nodes_to_regrid-1)-len1;
                sysA(1,2) = -lenN;
                sysb(1)=total_len-integrand(1)*len1...
                    -integrand(end)*lenN;
            end
            if endnode == bdry.node_count
                %force f_new = f_old on end side
                sysA(2,1) = 0; sysA(2,2) = 1;
                sysb(2) = integrand(end);
            else
                %int_0^{end_ref_len} f ds = len_new/nodes_to_regrid
                target_len = vecnorm(diff(...
                    bdry.boundary_nodes(endnode:(endnode+1),:)),2);
                sysA(2,2) = target_len*(nodes_to_regrid-1)-lenN;
                sysA(2,1) = -len1;
                sysb(2)=total_len-integrand(1)*len1...
                    -integrand(end)*lenN;
            end
            
            integ_min = 0;
            if isfield(bdry.regridding,'reweight_min')
                integ_min = bdry.regridding.reweight_min;
            end
            integ_new = max(sysA\sysb,integ_min);
            if startnode > 1
                %reweigh by lerp between
                % startnode and startnode + bufer_length
                integrand(1:buffer_length) = integ_new(1) + ...
                    (integrand(1:buffer_length) - integ_new(1))...
                    .*(0:(buffer_length-1))./buffer_length;
            end
            if endnode < bdry.node_count
                %reweigh by lerp between 
                % startnode and startnode + bufer_length
                integrand((end-buffer_length+1):end) = integ_new(2) + ...
                    (integrand((end-buffer_length+1):end) -integ_new(2))...
                    .*((buffer_length-1):-1:0)./buffer_length;
            end
            % Vnew = curve_renode_by_integral(V,integrand, ...
            %     nodes_to_regrid,[1;1;0]);


        end

        Vnew = curve_renode_by_integral(V,integrand, ...
            nodes_to_regrid,[1;1;0]);
        bdry.boundary_nodes(startnode:endnode,:) = Vnew(1:2,:)';
        bdry.characteristics.phi(startnode:endnode) = Vnew(3,:);
    end
    %force a reload of the characteristics
    did_change = 1;
end

