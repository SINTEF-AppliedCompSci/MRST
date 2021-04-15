function selection = selectTOFRegion(D, max_tof, min_tof, varargin)
%Select a subset of cells based on TOF criteria
%
% SYNOPSIS:
%   selectTOFRegion(D, max_tof, min_tof, ...);
%
% DESCRIPTION:
%   Selects a subset of a grid based on TOF information and selection
%   criteria.
%
% REQUIRED PARAMETERS:
%
%   D - Time of flight object computed by computeTOFandTracer.
%   max_tof - All TOF-regions higher than this will be discarded.
%   min_tof - All TOF-regions less than this will be discarded.
%
%
% OPTIONAL PARAMETERS:
%       'drain_wells' - List of draining (production) wells to consider.
%                       Other wells will be ignored
%       'flood_wells' - List of flooding (injecting) wells to consider.
%                       Other wells will be ignored
%       'set_op'	  - Operation op to perform, so that the computed
%                       selection is op(flood, drain) where flood and drain
%                       are the flooding and drainage volumes,
%                       respectively.
%       'near_well_max_tof' - Threshold for near well regions to include.
%                             Set to zero to disable including near-well
%                             regions.
%       'psubset' - User specified draining / production region to
%                   consider. When this is enabled, max_tof and min_tof
%                   will not be used to find the selected drainage region.
%       'isubset' - User specified flooding / injection region to consider.
%                   When this is enabled, max_tof and min_tof will not be
%                   used to find the selected flooding region.
%       'tracer_threshold' - Tracer-value threshold for including in
%                   region. For 'intersection', threshold applies to product 
%                   of inj/prod tracers. Default value: 0.05;
%
% RETURNS:
%
%   A list of cell indices
%
% EXAMPLE:
%
% SEE ALSO:
%   

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct('drain_wells', [],...
             'flood_wells', [],...
             'set_op', 'union',...
             'near_well_max_tof', 0.0,...
             'psubset', [],...
             'isubset', [], ...
             'tracer_threshold', 0.05);

opt = merge_options(opt, varargin{:});

% Drainage volumes
if (isempty(opt.psubset))
    data = D.tof(:,2);
    psubset = data >= min_tof & data <= max_tof;
else
    psubset = opt.psubset;
end

% Flooding volumes
if (isempty(opt.isubset))
    data = D.tof(:,1);
    isubset = data >= min_tof & data <= max_tof;
else
    isubset = opt.isubset;
end

%Find the wells we are interested in
%If no injection and no production wells are given, select all
if (isempty(opt.drain_wells) && isempty(opt.flood_wells))
    psubs = D.ppart >= 0;
    isubs = D.ipart >= 0;
else
    psubs  = sum(D.ptracer(:, opt.drain_wells), 2);
    isubs  = sum(D.itracer(:, opt.flood_wells), 2);
end

% Use set operation on the selection
switch(opt.set_op)
    %Union of injection and production volumes
    case 'union' 
        selection = (isubset & isubs > opt.tracer_threshold) | ...
                    (psubset & psubs > opt.tracer_threshold);
    
    %Intersection of injection and production volumes
    case 'intersection'
        selection = (isubset & psubset) & (isubs.*psubs > opt.tracer_threshold);
    
    %Injection / flood volumes only
    case 'flood'
        selection = (isubset & isubs > opt.tracer_threshold);
    
    %Production / drain volumes only
    case 'drain'
        selection = (psubset & psubs > opt.tracer_threshold);
    
    otherwise
        selection = [];
        warning('Invalid set operation.');
end

%Make selection include near-well regions up to near_well_max_tof
if (opt.near_well_max_tof > 0.0)
    isubset = D.tof(:,1) <= opt.near_well_max_tof;
    psubset = D.tof(:,2) <= opt.near_well_max_tof;
    near_well_selection = (isubset & isubs) | (psubset & psubs);
    selection  = near_well_selection | selection;
end

end