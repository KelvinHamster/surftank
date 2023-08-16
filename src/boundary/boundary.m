classdef (Abstract) boundary < matlab.mixin.Copyable
    %BOUNDARY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        intersection_CCW %struct managing intersection on the CCW side
        intersect_CW_link %ref to boundary on CW side

        formulation %how the surface is formulated, and its formulation

        boundary_nodes %list of nodes for this boundary, including corners
        BIE_shape_interp %interpolation struct for the boundary integral

        characteristics %struct of numerics
        node_count %number of nodes
        parent_sim %ref to bem_sim object
        regridding %struct handling how (if at all) regridding is handled
        meta %other boundary information
    end
    properties (Abstract)
        cond_type
        update_orders
        step_orders
        has_function_handles
    end
    
    methods        
        function init_surface(obj,varargin)
            % Populates fields in the boundary.
            % Call this in the constructor.

            % boundary described by explicit nodes, implicit function,
            % or parameterized eqn
            bd_types = {'nodes','implicit','parameterized'};
            func_interp_rules = {'sliding','sequential'};
            bdry_interp_rules = {'sliding','sequential','spline'};

            bd_type = 'nodes';
            node_count_ = 2;
            surface_nodes = 'none';
            surface_defn = 'none';
            BIE_interp_order = 1;
            func_interp_rule = 'sequential';
            bdry_interp_rule = 'sliding';
            bdry_interp_order = 1;

            p = inputParser;
            addParameter(p,'bd_type',bd_type,...
                @(x) any(validatestring(x,bd_types)));
            addParameter(p,'node_count',node_count_, @(x) isnumeric(x)...
                && isscalar(x) && (x > 0) && (mod(x,1) == 0) );
            addParameter(p,'surface_nodes',surface_nodes, ...
                @(x) isnumeric(x) && length(size(x))== 2 &&...
                dot(size(x),[0,1]) == 2 ||...
                ((ischar(x) || isstring(x)) && strcmp(x,'none')));
            addParameter(p,'surface_defn',surface_defn, ...
                @(x) ((ischar(x) || isstring(x)) && strcmp(x,'none')) ...
                || isa(x,'function_handle'));
            addParameter(p,'BIE_interp_order',BIE_interp_order, ...
                @(x) isnumeric(x)...
                && isscalar(x) && (x > 0) && (mod(x,1) == 0) );
            addParameter(p,'bdry_interp_order',bdry_interp_order, ...
                @(x) isnumeric(x)...
                && isscalar(x) && (x > 0) && (mod(x,1) == 0) );
            addParameter(p,'func_interp_rule',func_interp_rule,...
                @(x) any(validatestring(x,func_interp_rules)));
            addParameter(p,'bdry_interp_rule',bdry_interp_rule,...
                @(x) any(validatestring(x,bdry_interp_rules)));

            parse(p,varargin{:});
            bd_type = p.Results.bd_type;
            node_count_ = p.Results.node_count;
            surface_nodes = p.Results.surface_nodes;
            surface_defn = p.Results.surface_defn;
            BIE_interp_order = p.Results.BIE_interp_order;
            bdry_interp_order = p.Results.bdry_interp_order;
            bdry_interp_rule = p.Results.bdry_interp_rule;
            func_interp_rule = p.Results.func_interp_rule;
            %===

            if strcmp(bd_type,'nodes')
                if (ischar(surface_nodes) || isstring(surface_nodes))...
                        && strcmp(surface_nodes,'none')
                    error("'nodes' boundary type specified, but" +...
                        " no nodes given");
                end
                obj.boundary_nodes = surface_nodes;
                node_count_ = dot(size(surface_nodes),[1,0]);
            end
            obj.BIE_shape_interp = ...
                get_interpolation_struct(BIE_interp_order,0,1);
            
            obj.node_count = node_count_;

            obj.formulation = struct('type',bd_type,...
                'surface_defn',surface_defn,...
                'func_interp_rule',func_interp_rule,...
                'bdry_interp_rule',bdry_interp_rule);

            if strcmp(bd_type,'parameterized')
                obj.formulation.param_CW = 0;
                obj.formulation.param_CCW = 1;
            end

            obj.regridding = struct('mode','none','order',0,...
                'precede_update',0);

            if strcmp(bdry_interp_rule,'sliding')...
                    || strcmp(bdry_interp_rule,'sequential')
                obj.formulation.bdry_interp = ...
                    get_interpolation_struct(bdry_interp_order,0,1);
            end

            obj.meta = struct();

            obj.reassemble_nodes();
        end

        function reassemble_nodes(obj)
            % sets 'boundary_nodes' according to the formulation
            type = obj.formulation.type;
            if strcmp(type,'parameterized')
                n = obj.node_count;
                obj.boundary_nodes = zeros(n,2);
                param_CW = obj.formulation.param_CW;
                param_CCW = obj.formulation.param_CCW;
                for i=1:n
                    t = param_CW + (param_CCW-param_CW)*(i-1)/(n-1);
                    obj.boundary_nodes(i,:) = ...
                        obj.formulation.surface_defn(t);
                end
                
            elseif strcmp(type,'nodes')
                return;
            elseif strcmp(type,'implicit')
                error('Implicit node definition not yet supported.');
            else
                error('Unknown boundary formulation type (%s).',type);
            end
        end
        
    end

    methods (Abstract)
        on_update(obj,order);
            % Called when updating characteristics. order < 0 occurs before
            % the first BEM solve (getting phi,phi_n) and order < 1 occurs
            % before the second BEM solve (getting phi_t,phi_tn).

        on_step(obj,dt,order);
            % Called after post_syssolve, handled alongside
            % the free surface update
    end
end

