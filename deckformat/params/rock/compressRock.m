function rock = compressRock(rock, active)
%Compress rock properties to active cells only
%
% SYNOPSIS:
%   rock = compressRock(rock, active)
%
% PARAMETERS:
%   rock   - Rock data structure formed, e.g., in function
%            'initEclipseRock'.
%
%   active - Active cell map.  Logical or numeric indices into 'perm' (&c)
%            fields that denotes which cell propertis to extract.  Often
%            equal to the 'cells.indexMap' field of a grid_structure.
%
% SEE ALSO:
%   `grid_structure`, `initEclipseRock`.

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


    assert (isstruct(rock), 'rock must be a structure');

    nc = numel(rock.poro);
    rock = compressFields(rock, nc, active);
end

function str = compressFields(str, nc, subs)
    flds = fieldnames(str);
    
    for i = 1:numel(flds)
        f = flds{i};
        d = str.(f);
        if isstruct(d)
            str.(f) = compressFields(d, nc, subs);
        elseif isnumeric(d) && size(d, 1) == nc
            str.(f) = d(subs, :);
        end
    end
end
