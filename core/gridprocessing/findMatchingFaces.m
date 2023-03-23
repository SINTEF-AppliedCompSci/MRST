function N = findMatchingFaces(G)
%Find indices of faces sharing the same set of nodes in a grid
%
% SYNOPSIS:
%   N = findMatchingFaces(G)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
% RETURNS:
%   N       - An n x 2 array of face numbers. Each row represents a pair of
%             faces in G that share the same set of nodes.
%
% SEE ALSO:
%  `makeInternalBoundary`, `removeInternalBoundary`
%
% EXAMPLES:
%   % Make internal boundary, find matching faces, remove boundary again
%   G = cartGrid([3,2]);
%   G = makeInternalBoundary(G, 13);
%   N = findMatchingFaces(G);
%   G = removeInternalBoundary(G, N);

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

    % Reshape face nodes to matrix with zero-padding
    n     = reshape(diff(G.faces.nodePos), [], 1);
    nodes = zeros([G.faces.num, max(n)]);
    ix    = sub2ind(size(nodes)                 , ...
                 rldecode(1:G.faces.num, n, 2).', ...
                 reshape(mcolon(1, n), [], 1)   );
    nodes(ix) = G.faces.nodes(:,1);
    % Ensure equal node ordering for all faces
    nodes = sort(nodes, 2);
    % Find matching rows
    [tmp, faces] = sortrows(nodes);
    [tmp, num  ] = rlencode(tmp,1); %#ok
    % Check that at most two faces share the same set of nodes
    assert(max(num) <= 2, ['This function currently does not support ', ...
         'the case when more than two faces share the same set of nodes']);
    fix = rldecode(num==2, num, 1);
    % Reshape to pairs of matching faces
    N = reshape(faces(fix), 2, [])';
    
end