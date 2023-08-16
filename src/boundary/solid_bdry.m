classdef solid_bdry < boundary
    %SOLID_BDRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cond_type = "neumann"
        update_orders = [-1]
        step_orders = []
        has_function_handles = 0
    end
    
    methods
        function obj = solid_bdry(varargin)
            %SOLID_BDRY Construct an instance of this class
            %   Detailed explanation goes here
            obj.init_surface(varargin{:});
        end

        function on_update(obj, order)
            obj.characteristics.phi_n = ...
                zeros(1,obj.node_count);
            obj.characteristics.phi_tn = ...
                zeros(1,obj.node_count);
        end

        function on_step(obj,dt,order)
        end
    end
end

