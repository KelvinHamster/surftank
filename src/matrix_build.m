function [ret1,ret2] = matrix_build(sim, varargin)
%MATRIX_BUILD Generates the matrix / vectors for the BIE linear system

build_A = 1;
build_b = 1;
b_fieldsrc_u = 'phi';
b_fieldsrc_q = 'phi_n';

p = inputParser;
addOptional(p,'build_A',build_A,@(x) islogical(x) || (x==0) || (x==1));
addOptional(p,'build_b',build_b,@(x) islogical(x) || (x==0) || (x==1));
addOptional(p,'b_fieldsrc_u',b_fieldsrc_u,@(x) ischar(x) || isstring(x));
addOptional(p,'b_fieldsrc_q',b_fieldsrc_q,@(x) ischar(x) || isstring(x));
parse(p, varargin{:});

build_A = p.Results.build_A;
build_b = p.Results.build_b;
b_fieldsrc_u = p.Results.b_fieldsrc_u;
b_fieldsrc_q = p.Results.b_fieldsrc_q;




%Kn, Kd are N x S, (N= num singularities, S= N+num_bdrys = num shapefuncs)
S = 1;
NUM_BDRYS = length(sim.boundaries);
shape_inds = zeros(1,NUM_BDRYS);

%alpha(x), boundary shape param, per singularity
alpha = 0;

for i=1:NUM_BDRYS
    shape_inds(i) = S;
    S = S + sim.boundaries{i}.node_count; %start of next bdry
    alpha = alpha - sum(sim.boundaries{i}.characteristics.K_n,2);
end
S = S - 1;%now want end of last bdry
N = length(alpha);

if build_A
    A = zeros(S,S);
end
if build_b
    b = zeros(S,1);  %for phi
end

for i=1:length(sim.boundaries)
    bdry = sim.boundaries{i};
    %build first N rows: based on BIE:
    %       alpha(l)*U_l = sum_j(q_j * Kd(l,j) - u_j * Kn(l,j))
    % for U_l on doublenodes choose one u_j for it to correspond with
    n = bdry.node_count;
    if strcmp(bdry.cond_type,'neumann')
        %alpha(l)*U_l + sum_j(u_j * Kn(l,j)) = sum_j(q_j * Kd(l,j))
        if build_A
            A(1:N,shape_inds(i):(shape_inds(i)+n-1)) = ...
                bdry.characteristics.K_n;
            for l=shape_inds(i):(shape_inds(i)+n-2)
                %shape functions: (n1)  -(n2)  -...-(ni)
                %singularities:   (n1-1)-(n2-1)-...-(ni-1)
                %                                     ^^ offset by i-1
                A(l-i+1,l) = A(l-i+1,l) + alpha(l-i+1);
            end
        end
        %====
        if build_b
            b(1:N) = b(1:N)...
                + bdry.characteristics.K_d * ...
                getfield(bdry.characteristics,b_fieldsrc_q)';
        end

    elseif strcmp(bdry.cond_type,'dirichlet')
        %-sum_j(q_j * Kd(l,j)) = -alpha(l)*U_l - sum_j(u_j * Kn(l,j))
        if build_A
            A(1:N,shape_inds(i):(shape_inds(i)+n-1)) = ...
                -bdry.characteristics.K_d;
        end
        %====
        if build_b
            u = getfield(bdry.characteristics, b_fieldsrc_u)';
            b(1:N) = b(1:N) - bdry.characteristics.K_n * u;
            
            for l=1:(n-1)
                %shape functions: (n1)  -(n2)  -...-(ni)
                %singularities:   (n1-1)-(n2-1)-...-(ni-1)
                %                                     ^^ offset by i-1
                s = shape_inds(i) + l - i;
                b(s) = b(s) - alpha(s) * u(l);
            end
        end

    else
        error("Boundary %d: invalid boundary condition type (%s)",...
            i,bdry.cond_type);
    end

    %build last NUM_BDRYS rows: based on doublenode resolution
    bdCW = bdry.intersect_CW_link;
    if strcmp(bdry.cond_type,'neumann')
        if strcmp(bdCW.cond_type,'neumann')
            %neumann-neumann: equate phi on both sides of DN
            if build_A
                A(N+i,shape_inds(i)) = 1;
                A(N+i,shape_inds(i)-1 + S*(shape_inds(i)==1)) = -1;
            end
        elseif strcmp(bdCW.cond_type,'dirichlet')
            %dirichlet-neumann: equate phi on both sides of DN (CW known)
            if build_A
                A(N+i,shape_inds(i)) = 1;
            end
            if build_b
                u = getfield(bdCW.characteristics,b_fieldsrc_u);
                b(N+i) = u(end);
            end
        end
    elseif strcmp(bdry.cond_type,'dirichlet')
        if strcmp(bdCW.cond_type,'neumann')
            %neumann-dirichlet: equate phi on both sides of DN (CCW known)
            if build_A
                A(N+i,shape_inds(i)-1 + S*(shape_inds(i)==1)) = 1;
            end
            if build_b
                %u was defined already in prev if/elseif block
                b(N+i) = u(1);
            end
        elseif strcmp(bdCW.cond_type,'dirichlet')
            error("Boundary %d-%d intersection: dirichlet-dirichlet" + ...
                " doublenode not implemented yet!",i,i-1+NUM_BDRYS*(i==1))
        end
    end
end

if build_A
    ret1 = A;
    if build_b
        ret2 = b;
    end
else
    if build_b
        ret1 = b;
    end
end

end

