function FS_add_absorbing_layer(sim, fs_bdry, x_start, x_end, varargin)
% Adds an absorbing layer to the given free surface
%
% A global event is appended to the simulation that adds a correcting
% pressure to kill high frequency waves. This force starts at x_start,
% ramping up from zero to full force at x_end, in a power curve with
% exponent mu.
%
% sim - a reference to the bem_sim object
% fs_bdry - the index of the free surface to give the absorbing layer to
% x_start - the starting position of the absorbing layer
% x_end - the ending position of the absorbing layer
%
% Optional (keyword) arguments are mu, nu0, and h0.
%
%
% The pressure term is set to nu0 * phi_n * sqrt(gh)


mu = 2;
nu0 = 1;
h = 1;



p = inputParser;
addOptional(p,'mu',mu, @(x) ...
   isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'nu0',nu0, @(x) ...
   isnumeric(x) && isscalar(x) && (x > 0) );
addOptional(p,'h',h, @(x) ...
   isnumeric(x) && isscalar(x) && (x > 0) );

parse(p, varargin{:});

mu = p.Results.mu;
nu0 = p.Results.nu0;
h = p.Results.h;

if ~exist("FS_absorbing_layers","var")
    sim.meta.FS_absorbing_layers = struct("active",0,...
        "x_start",0,"x_end",0,...
        "mu",0,"nu0",0,"h",0);
    sim.append_update_event(@update_pressure, 0.6);
end
sim.meta.FS_absorbing_layers(fs_bdry).active = 1;
sim.meta.FS_absorbing_layers(fs_bdry).x_start = x_start;
sim.meta.FS_absorbing_layers(fs_bdry).x_end = x_end;
sim.meta.FS_absorbing_layers(fs_bdry).mu = mu;
sim.meta.FS_absorbing_layers(fs_bdry).nu0 = nu0;
sim.meta.FS_absorbing_layers(fs_bdry).h = h;


end

function update_pressure(sim)
    for fs_bdry = 1:length(sim.meta.FS_absorbing_layers)
        active = sim.meta.FS_absorbing_layers(fs_bdry).active;
        if isempty(active) || ~active
            continue
        end
        x_start = sim.meta.FS_absorbing_layers(fs_bdry).x_start;
        x_end = sim.meta.FS_absorbing_layers(fs_bdry).x_end;
        mu = sim.meta.FS_absorbing_layers(fs_bdry).mu;
        nu0 = sim.meta.FS_absorbing_layers(fs_bdry).nu0;
        h = sim.meta.FS_absorbing_layers(fs_bdry).h;

        bdry = sim.boundaries{fs_bdry};
        pres = bdry.characteristics.phi_n...
                .* nu0 .* sqrt(sim.meta.g * h) .* ...
                max(0, ...
                (bdry.boundary_nodes(:,1)' - x_start)/(x_end-x_start) ...
                ).^mu;
        
        bdry.characteristics.dphidt =  bdry.characteristics.dphidt - pres;
        bdry.characteristics.phi_t =  bdry.characteristics.phi_t - pres;
    end
end