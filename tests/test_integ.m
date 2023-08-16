function test_integ(sim)
%TEST_INTEG performs integral tests on the given sim object

reltol = 1e-2;
abstol = 1e-3;


BIE_integ_calc(sim)

[singularities,sing_inds] = sim.get_nodelist();
N = size(singularities,1);

for Bnum=1:length(sim.boundaries)
shouldskip = 0;
fprintf("---[Boundary %d]---\n",Bnum);
bdry = sim.boundaries{Bnum};
n = bdry.node_count;

fprintf("  shape function rule: '%s'",...
    bdry.formulation.func_interp_rule);
if strcmp(bdry.formulation.func_interp_rule,'sliding') ||...
        strcmp(bdry.formulation.func_interp_rule,'sequential')
    fprintf(" (M=%d)\n",bdry.BIE_shape_interp.M);
else
    fprintf("  (UNSUPPORTED BY TEST)\n")
    shouldskip = 1;
end
fprintf("     bdry interp rule: '%s'",bdry.formulation.bdry_interp_rule);
x_b = bdry.boundary_nodes(:,1);
z_b = bdry.boundary_nodes(:,2);
if strcmp(bdry.formulation.bdry_interp_rule,'sliding') ||...
        strcmp(bdry.formulation.bdry_interp_rule,'sequential')
    fprintf(" (M=%d)\n",bdry.formulation.bdry_interp.M);
elseif strcmp(bdry.formulation.bdry_interp_rule,'spline')
    splx = csape(1:n, x_b).coefs;
    splz = csape(1:n, z_b).coefs;
    fprintf("\n");
else
    fprintf("  (UNSUPPORTED BY TEST)\n")
    shouldskip = 1;
end

if shouldskip
    fprintf("Skipped evaluations!\n");
    continue
end

if isfield(bdry.formulation,'bdry_interp')
    M_bdry = bdry.formulation.bdry_interp.M;
    Lmat_bdry = bdry.formulation.bdry_interp.Lagrange;
    Lpmat_bdry = bdry.formulation.bdry_interp.Lagrange_prime;
end
M = bdry.BIE_shape_interp.M;
Lmat = bdry.BIE_shape_interp.Lagrange;
Lpmat = bdry.BIE_shape_interp.Lagrange_prime;

Kd = zeros(N,n);
Kn = zeros(N,n);
%loop over boundary segments
for j=1:(n-1)
    %get the boundary
    if strcmp(bdry.formulation.bdry_interp_rule,'sliding')
        half_M = floor(M_bdry/2);
        seg_start = max(1,min(n-M_bdry,j-half_M));
        segx = x_b(seg_start:(seg_start + M_bdry))';
        segz = z_b(seg_start:(seg_start + M_bdry))';
        offset = j - seg_start;
        gamma_x = @(t) (segx * Lmat_bdry)...
            * (((offset + t)./M_bdry).^((0:M_bdry)'));
        gamma_z = @(t) (segz * Lmat_bdry)...
            * (((offset + t)./M_bdry).^((0:M_bdry)'));
        gammap_x = @(t) (segx * Lpmat_bdry)...
            * (((offset + t)./M_bdry).^((0:M_bdry-1)'));
        gammap_z = @(t) (segz * Lpmat_bdry)...
            * (((offset + t)./M_bdry).^((0:M_bdry-1)'));

    elseif strcmp(bdry.formulation.bdry_interp_rule,'sequential')
        if mod(n-1,M_bdry) ~= 0
            error("Boundary %d: sequential boundary rule (M=%d) "+ ...
                "is incompatible with %d nodes!",...
                bdry_index,M_bdry,n);
        end
        seg_start = M_bdry * floor((j - 1)/M_bdry) + 1;
        segx = x_b(seg_start:(seg_start + M_bdry))';
        segz = z_b(seg_start:(seg_start + M_bdry))';
        offset = j - seg_start;
        gamma_x = @(t) (segx * Lmat_bdry)...
            * (((offset + t)./M_bdry).^((0:M_bdry)'));
        gamma_z = @(t) (segz * Lmat_bdry)...
            * (((offset + t)./M_bdry).^((0:M_bdry)'));
        gammap_x = @(t) (segx * Lpmat_bdry)...
            * (((offset + t)./M_bdry).^((0:M_bdry-1)'));
        gammap_z = @(t) (segz * Lpmat_bdry)...
            * (((offset + t)./M_bdry).^((0:M_bdry-1)'));
        
    elseif strcmp(bdry.formulation.bdry_interp_rule,'spline')
        ax = splx(j,:); az = splz(j,:);
        gamma_x = @(t) ax(1).*t.^3 + ax(2).*t.^2 + ax(3).*t + ax(4);
        gamma_z = @(t) az(1).*t.^3 + az(2).*t.^2 + az(3).*t + az(4);
        gammap_x = @(t) 3.*ax(1).*t.^2 + 2.*ax(2).*t + ax(3);
        gammap_z = @(t) 3.*az(1).*t.^2 + 2.*az(2).*t + az(3);
    end
    
    %build shape functions
    if strcmp(bdry.formulation.func_interp_rule,'sliding')
        half_M = floor(M/2);
        shape_start = max(1,min(n-M,j-half_M));
    elseif strcmp(bdry.formulation.func_interp_rule,'sequential')
        if mod(n-1,M)~=0
            error("Boundary %d: sequential function rule (M=%d)"...
            +" is incompatible with %d nodes!",...
            bdry_index,M,n);
        end
        shape_start = M * floor((j - 1)/M) + 1;
    end
    shape_off = j - shape_start;

    %loop over shape functions that are nonzero on this interval
    for i = 1:(M+1)
        shape = @(t) Lmat(i,:) * (((shape_off + t)./M).^((0:M)'));

        for l = 1:N
            [Kd_eval,Kn_eval] = test_calc_reference_integral(...
                singularities(l,1),singularities(l,2),...
                gamma_x,gamma_z,gammap_x,gammap_z,shape,0,1);
            Kd(l,shape_start + i - 1)=Kd(l,shape_start + i - 1) + Kd_eval;
            Kn(l,shape_start + i - 1)=Kn(l,shape_start + i - 1) + Kn_eval;
        end
    end
end

%==============Compare answers
errs_Kd = Kd - bdry.characteristics.K_d;
errs_Kn = Kn - bdry.characteristics.K_n;

%pass if |err| < abstol and |err| / max(|true|,1e-5) < reltol

to_spec_Kd = abs(errs_Kd) < abstol & abs(errs_Kd)<reltol*max(abs(Kd),1e-5);
to_spec_Kn = abs(errs_Kn) < abstol & abs(errs_Kn)<reltol*max(abs(Kn),1e-5);

%loop shape functions
for j=1:n
    if all(to_spec_Kd(:,j)) && all(to_spec_Kn(:,j))
        continue
    end
    fprintf("Node %d failures:\n",j);
    for l=1:N
        relative = (l - sing_inds(Bnum) + 1) - j;
        if (~to_spec_Kd(l,j)) || (~to_spec_Kn(l,j))
            fprintf("  sing %d (j%+d)\n",l,relative);
            if ~to_spec_Kd(l,j)
                fprintf("    Kd: %15.7e should be %15.7e (%+e)\n",...
                    bdry.characteristics.K_d(l,j),Kd(l,j),errs_Kd(l,j));
            end
            if ~to_spec_Kn(l,j)
                fprintf("    Kn: %15.7e should be %15.7e (%+e)\n",...
                    bdry.characteristics.K_n(l,j),Kn(l,j),errs_Kn(l,j));
            end
            fprintf("\n");
        end
    end
end

end


end