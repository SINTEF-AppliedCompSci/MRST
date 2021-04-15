function v = indirectionSub(i, valuePos, values)
% Look-up in index map of the type G.cells.facePos, G.faces.nodePos, etc
%
% SYNOPSIS:
%   v = indirectionSub(i, valuePos, values)
%
% PARAMETERS:
%   i           - Indices where we want values, for instance a list of
%                 cells if valuePos is G.cells.facePos.
%
%   valuePos    - Mapping from i to values.
%
%   values      - The actual values
%
%
% RETURNS:
%   v           - The set of values corresonding to i.
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


    ind = mcolon(valuePos(i),...
                 valuePos(i+1)-1)';
    v = values(ind);
end
