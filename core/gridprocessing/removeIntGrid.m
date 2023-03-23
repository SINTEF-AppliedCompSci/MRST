function G = removeIntGrid(G)
%Cast any grid fields that are presently int32 to double
%
% SYNOPSIS:
%   G = removeIntGrid(G);
%
% PARAMETERS:
%   G - Grid with fields that are possibly int32.
%
% RETURNS:
%   G - Grid where int32 has been removed from cells/faces subfields
%
% NOTE:
%   Previously `int32` was used in the grid structure to conserve memory.
%   This routine converts grids of the old type to the new one.
%
% SEE ALSO:
%   `grid_structure`

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

    intFields = {'cells', 'facePos'; ...
                 'cells', 'faces'; ...
                 'cells', 'indexMap'; ...
                 'faces', 'neighbors'; ...
                 'faces', 'nodePos'; ...
                 'faces', 'nodes'; ...
                 'faces', 'tag' ...
                };
    for i = 1:size(intFields, 1)
        c = intFields{i, 1};
        name = intFields{i, 2};
        if isfield(G.(c), name);
            G.(c).(name) = double(G.(c).(name));
        end
    end
end
