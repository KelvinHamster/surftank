function intersection_build_default(CW,CCW)
%Initializes the intersection structure between CW and CCW
%   CW is the boundary on the clockwise-side
%   CCW is the boundary on the counter-clockwise-side
%


%resolutions:
% midpoint - place corner node at the average of either side
% match_CW - place corner node at the node of the CW side
% match_CCW - place corner node at the node of the CCW side
% curve_intersect - representing both boundaries as curves, place corner
%                   node at the intersection between the two curves
% ignore - do not do any resolution on this intersection
CW.intersection_CCW = struct(...
    "CCW_link",CCW, ...
    "resolution","ignore"... << default
);
CCW.intersect_CW_link = CW;

%==========================================================================
% Free surface adjacent to node formulation:
% Match to free surface side!
if isa(CW,"free_surface_bdry") && strcmp(CCW.formulation.type,"nodes")
    CW.intersection_CCW.resolution = "match_CW";
    return;
end
if isa(CCW,"free_surface_bdry") && strcmp(CW.formulation.type,"nodes")
    CW.intersection_CCW.resolution = "match_CCW";
    return;
end

%==========================================================================
% Piston-type boundaries adjacent to solid boundaries:
% Match to piston!
if is_piston_type(CW) && isa(CCW,"solid_bdry")
    CW.intersection_CCW.resolution = "match_CW";
    return;
end
if is_piston_type(CCW) && isa(CW,"solid_bdry")
    CW.intersection_CCW.resolution = "match_CCW";
    return;
end

end

function val = is_piston_type(obj)
val = isa(obj,"vertical_wavemaker_bdry") || ...
    isa(obj,"piston_absorber_bdry");
end