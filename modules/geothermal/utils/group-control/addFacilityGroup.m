function groups = addFacilityGroup(groups, members, varargin)
%Insert a facility group into the simulation model.
%
% SYNOPSIS:
%   groups = addFacilityGroup(groups, members)
%   groups = addFacilityGroup(groups, members, 'pn', pv, ...)
%
% REQUIRED PARAMETERS:
%   groups  - Group structure or empty if no other gopus exist.
%             Updated upon return.
%
%   members - Cell array of names of the wells in the group.
%
% OPTIONAL PARAMETERS:
%   name   - Name of facility group.
%            Default value: `sprintf('G%d', numel(groups) + 1)`
%
%   ftype  - Cell array of types for each group memer.
%            Default value: 'well'
%
%   type   - Type of group control. Currently supported types are 'bhp' and
%            'rate'. Support for 
%
%   val    - Group control target value
%
%   compi  - Fluid phase composition for injection group. Vector of phase
%            volume fractions.
%
%   sign   - Group type: Production (sign = -1) or Injection (sign = 1).
%            Default value: 0 (Undetermined sign. Will be derived from
%            rates if possible).
%
%   lims   - Group limits

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    group = struct( ...
        'name'  , sprintf('G%d', numel(groups) + 1), ...
        'ftype' , {{}}                             , ...
        'type'  , 'none'                           , ...
        'val'   , nan                              , ...
        'compi' , 1                                , ...
        'sign'  , 0                                , ...
        'status', true                             , ...
        'lims'  , []                                 ...
    );
    group = merge_options(group, varargin{:});
    group.members = members;
    
    group = checkInput(group);
    
    groups = [groups; group];

end

function group = checkInput(group)

    if ~iscell(group.members), group.members = {group.members}; end
    nm = numel(group.members);
    if isempty(group.ftype), group.ftype = 'well'; end
    if ~iscell(group.ftype), group.ftype = {group.ftype}; end
    if numel(group.ftype) == 1
        group.ftype = repmat(group.ftype, nm, 1);
    end
    assert(numel(group.ftype) == nm);

end