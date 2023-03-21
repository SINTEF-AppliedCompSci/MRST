function s = cellNodes(g)
%Extract local-to-global vertex numbering for grid cells.
%
% SYNOPSIS:
%   cn = cellNodes(G)
%
% PARAMETERS:
%   G - Grid data structure geometrically discretising a reservoir model.
%
% RETURNS:
%   cn - An m-by-3 array mapping cell numbers to vertex numbers.
%        Specifically, if `cn(i,1)==j` and `cn(i,3)==k`, then global vertex
%        `k` is one of the corners of cell `j`.  The local vertex number of
%        global node `k` within cell `j` may be computed using the
%        following statements::
%
%           n_vert = accumarray(cn(:,1), 1, [G.cells.num, 1]);
%           offset = rldecode(cumsum([0; n_vert(1:end-1)]), n_vert);
%           loc_no = (1 : size(cn,1)).' - offset;
%
%        This calculation is only for illustration purposes.  It is assumed
%        that the local vertex number will be implicitly available in most
%        applications (e.g., finite element methods).
%
%        Alternatively, the local vertex number is found in `cn(i,2). ;)`.
%        In a corner-point grid, local vertex number of zero indicates that
%        the vertex is not part of the original 8 vertices.
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


   dim = numel(g.cartDims);% need to be a corner point grid to get this working

   cellNo = rldecode(1:g.cells.num, double(diff(g.cells.facePos)), 2) .';
   cf     = double(g.cells.faces);
   ne     = double(diff(g.faces.nodePos));

   % expand table of cell number and halfface info
   a   = rldecode([cellNo, cf(:,2)], ne(cf(:,1)), 1);

   % the index i expands g.faces.nodes to 'halffaceNodes'
   i   = mcolon(g.faces.nodePos(cf(:,1)    ), ...
                g.faces.nodePos(cf(:,1) + 1) - 1) .';


   % binary magic applied to halfface-tag to discover nodes shared by
   % three faces.  Remember that sparse adds contributions when
   % (cell, node)-pairs are repeated.
   numbers = 2 .^ (0 : max(a(:,2)) - 1) .';
   mat     = sparse(a(:,1), double(g.faces.nodes(i,1)), numbers(a(:,2)));

   [cellNumber, nodeNumber, tag] = find(mat);

if dim==2
   % For Cartesian grid, these are halfface tags for each of eight
   % vertices
   k = [1, 3; ...  % west-south  = 2^0+2^2+2^4
        2, 3; ...  % east-south
        1, 4; ...  % west-north
        2, 4];     % east-north

elseif dim==3
   % For Cartesian grid, these are halfface tags for each of eight
   % vertices
   k = [1, 3, 5; ...  % west-south-down  = 2^0+2^2+2^4
        2, 3, 5; ...  % east-south-down
        1, 4, 5; ...  % west-north-down
        2, 4, 5; ...  % east-north-down
        1, 3, 6; ...  % west-south-up
        2, 3, 6; ...  % east-south-up
        1, 4, 6; ...  % west-north-up
        2, 4, 6];     % east-north-up
end
   localnode = sparse(sum(2 .^ (k-1), 2), 1, 1 : size(k,1));

   s = sortrows([cellNumber(:), full(localnode(tag(:))), nodeNumber(:)]);
   assert (all(s(:,3) > 0));
end
