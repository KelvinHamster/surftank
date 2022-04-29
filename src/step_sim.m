function obj = step_sim(obj)
% Does one time step of the simulation with static neumann boundaries; stand-alone and assuming P = 0 (no surface pressure)
%

obj = obj.update_characteristics();
obj = obj.update_history();


x = [obj.boundary.x obj.boundary.x(1)];
z = [obj.boundary.z obj.boundary.z(1)];
doublenode = obj.boundary.doublenode;

FS_i = obj.boundary.FS_start;
FS_f = obj.boundary.FS_end;
phi_FS = obj.boundary.phi_FS;
dt = obj.stepping.dt;

% expected variables to exist:
%
% x: x-coordinates of the collocation points
% z: z-coordinates of the collocation points
% doublenode: array of logic values for if nodes i and i+1 are considered a double-node
%
% FS_i: first point of the free surface
% FS_f: last point of the free surface
%
% phi_FS: phi-values on the free surface phi_FS(1) corresponds with x(FS_i)
%
% dt: size of the timestep


N = length(x) - 1;
FS_size = length(phi_FS);



dxdt_FS = obj.characteristics_.dxdt_FS;
dzdt_FS = obj.characteristics_.dzdt_FS;
d2xdt2_FS = obj.characteristics_.d2xdt2_FS;
d2zdt2_FS = obj.characteristics_.d2zdt2_FS;
dphidt_FS = obj.characteristics_.dphidt_FS;
d2phidt2_FS = obj.characteristics_.d2phidt2_FS;

%do the actual timestep

%replace later; this is to test theory
%dphidt_FS = theo_dphidt(t,x(FS_i:FS_f),z(FS_i:FS_f));
%d2phidt2_FS = theo_d2phidt2(t,x(FS_i:FS_f),z(FS_i:FS_f));
%dxdt_FS = theo_u(t,x(FS_i:FS_f),z(FS_i:FS_f));
%d2xdt2_FS = theo_d2xdt2(t,x(FS_i:FS_f),z(FS_i:FS_f));
%dzdt_FS = theo_w(t,x(FS_i:FS_f),z(FS_i:FS_f));
%d2zdt2_FS = theo_d2zdt2(t,x(FS_i:FS_f),z(FS_i:FS_f));

switch obj.stepping.timestep_method
    
    case 'taylor2'
        x(FS_i:FS_f) = x(FS_i:FS_f) + dt .* dxdt_FS + (dt^2)/2 .* d2xdt2_FS;
        z(FS_i:FS_f) = z(FS_i:FS_f) + dt .* dzdt_FS + (dt^2)/2 .* d2zdt2_FS;


        phi_FS = phi_FS + dt .* dphidt_FS + (dt^2)/2 .* d2phidt2_FS;
    case 'abm4'
        order = 4;
        if ~isfield(obj.characteristics_,'time_derivs')
            obj.characteristics_.time_derivs = zeros(1,3,FS_size);
            obj.characteristics_.time_derivs(1,1,:) = dxdt_FS;
            obj.characteristics_.time_derivs(1,2,:) = dzdt_FS;
            obj.characteristics_.time_derivs(1,3,:) = dphidt_FS;
        else
            d_index = length(obj.characteristics_.time_derivs(:,1,1));
            if d_index >= order
                d_index = order;
            else
                obj.characteristics_.time_derivs(d_index+1,:,:) = 0;
            end
            
            obj.characteristics_.time_derivs(2:min(d_index+1,order),:,:) = obj.characteristics_.time_derivs(1:min(d_index+1,order)-1,:,:);
            obj.characteristics_.time_derivs(1,1,:) = dxdt_FS;
            obj.characteristics_.time_derivs(1,2,:) = dzdt_FS;
            obj.characteristics_.time_derivs(1,3,:) = dphidt_FS;
        end
        
        
        if length(obj.characteristics_.time_derivs(:,1,1)) < order
            %default to euler

            x(FS_i:FS_f) = x(FS_i:FS_f) + dt .* dxdt_FS;% + (dt^2)/2 .* d2xdt2_FS;
            z(FS_i:FS_f) = z(FS_i:FS_f) + dt .* dzdt_FS;% + (dt^2)/2 .* d2zdt2_FS;

            phi_FS = phi_FS + dt .* dphidt_FS + (dt^2)/2 .* d2phidt2_FS;
        else
            step = squeeze(...
                   55 .* obj.characteristics_.time_derivs(1,:,:) ...
                 - 59 .* obj.characteristics_.time_derivs(2,:,:) ...
                 + 37 .* obj.characteristics_.time_derivs(3,:,:) ...
                 -  9 .* obj.characteristics_.time_derivs(4,:,:) );
             
             x(FS_i:FS_f) = x(FS_i:FS_f) + (dt/24) .* step(1,:);
             z(FS_i:FS_f) = z(FS_i:FS_f) + (dt/24) .* step(2,:);
             phi_FS = phi_FS + (dt/24) .* step(3,:);
        end
    case 'abm2'
        order = 2;
        if ~isfield(obj.characteristics_,'time_derivs')
            obj.characteristics_.time_derivs = zeros(1,3,FS_size);
            obj.characteristics_.time_derivs(1,1,:) = dxdt_FS;
            obj.characteristics_.time_derivs(1,2,:) = dzdt_FS;
            obj.characteristics_.time_derivs(1,3,:) = dphidt_FS;
        else
            d_index = length(obj.characteristics_.time_derivs(:,1,1));
            if d_index >= order
                d_index = order;
            else
                obj.characteristics_.time_derivs(d_index+1,:,:) = 0;
            end
            
            obj.characteristics_.time_derivs(2:min(d_index+1,order),:,:) = obj.characteristics_.time_derivs(1:min(d_index+1,order)-1,:,:);
            obj.characteristics_.time_derivs(1,1,:) = dxdt_FS;
            obj.characteristics_.time_derivs(1,2,:) = dzdt_FS;
            obj.characteristics_.time_derivs(1,3,:) = dphidt_FS;
        end
        
        
        if length(obj.characteristics_.time_derivs(:,1,1)) < order
            %default to euler

            x(FS_i:FS_f) = x(FS_i:FS_f) + dt .* dxdt_FS;% + (dt^2)/2 .* d2xdt2_FS;
            z(FS_i:FS_f) = z(FS_i:FS_f) + dt .* dzdt_FS;% + (dt^2)/2 .* d2zdt2_FS;

            phi_FS = phi_FS + dt .* dphidt_FS + (dt^2)/2 .* d2phidt2_FS;
        else
            step = squeeze(...
                   (1.5) .* obj.characteristics_.time_derivs(1,:,:) ...
                 - (0.5) .* obj.characteristics_.time_derivs(2,:,:) );
             
             x(FS_i:FS_f) = x(FS_i:FS_f) + dt .* step(1,:);
             z(FS_i:FS_f) = z(FS_i:FS_f) + dt .* step(2,:);
             phi_FS = phi_FS + dt .* step(3,:);
        end
    case 'euler'
        
        x(FS_i:FS_f) = x(FS_i:FS_f) + dt .* dxdt_FS;
        z(FS_i:FS_f) = z(FS_i:FS_f) + dt .* dzdt_FS;

        phi_FS = phi_FS + dt .* dphidt_FS;
    otherwise
        error('Unknown time-stepping procedure "%s"!',obj.stepping.timestep_method);
        
end

FSim1 = FS_i - 1 + N*(FS_i == 1);
FSfp1 = FS_f + 1 - N*(FS_f == N);

x(FSim1) = x(FS_i); z(FSim1) = z(FS_i);
x(FSfp1) = x(FS_f); z(FSfp1) = z(FS_f);


obj.boundary.x = x(1:N);
obj.boundary.z = z(1:N);
obj.boundary.doublenode = doublenode;
obj.boundary.phi_FS = phi_FS;

obj.stepping.t = obj.stepping.t + dt;

%obj = obj.update_characteristics();
obj = obj.update_history();
end

