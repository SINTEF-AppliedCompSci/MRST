function [data, pos] = removeFromPackedData(pos, data, remove)
%Remove all copies of specific elements from packed data structure
%
% SYNOPSIS:
%   [data, p] = removeFromPackedData(p, data, elms)
%
% PARAMETERS:
%   p    - Indirection map into packed data array.  Assumed to satisfy the
%          condition that the data elements pertaining to a numbered
%          entity, e, are stored in rows
%
%              p(e) : p(e + 1) - 1
%
%          of 'data'.  NUMEL(p) is consequently assumed to be one greater
%          than the total number of entities represented in 'data'.
%
%   data - Packed data array elements.
%
%   elms - List of specific elements, all of whose copies should be removed
%          from the array pair (data,p).
%
% RETURNS:
%   data, p - Updated packed data array pair after removing all copies of
%             the elements, 'elms', from the input data array pair.
%
% EXAMPLE:
%   % Remove faces [100, 101, ..., 110] from packed face array pair
%   %   (G.cells.faces, G.cells.facePos)
%   %
%   [G.cells.faces, G.cells.facePos] = ...
%       removeFromPackedData(G.cells.facePos, G.cells.faces, 100:110);
%
%   % Remove nodes [100, 101, ..., 110] from packed node array pair
%   %   (G.faces.nodes, G.faces.nodePos)
%   %
%   [G.faces.nodes, G.faces.nodePos] = ...
%       removeFromPackedData(G.faces.nodePos, G.faces.nodes, 100:110);

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


   n         = numel(pos)-1;
   [r, p]    = reverseLookup(n, pos, data, remove);
   data(p,:) = [];
   pos       = pos - cumsum([0; accumarray(r,1,[n,1])]);
end

%--------------------------------------------------------------------------

function [r, p] = reverseLookup(n, pos, data, v)
   i    = false(max(data(:,1)), 1);
   i(v) = true;
   no   = rldecode(1:n, diff(pos), 2).';
   p    = i(data(:,1));
   r    = no(p);
end
