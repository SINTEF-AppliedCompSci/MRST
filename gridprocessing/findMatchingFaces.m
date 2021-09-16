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

    % Reshape face nodes to matrix with nan-padding
    nodes = G.faces.nodes;
    n     = diff(G.faces.nodePos);
    stop  = cumsum(n);
    start = stop - n + 1;
    nmax  = max(n);
    % Get index and map for values to be replaced by nans
    ix     = bsxfun(@plus, start, 0:nmax-1);
    nanmap = bsxfun(@gt, ix, stop);
    % Pad imssing dimension at the end with nans
    nodes(end+1:end+1+nmax-1) = nan;
    % Expand to matrix
    nodes = nodes(ix);
    % Replace extra values with nan
    nodes(nanmap) = nan;
    % Find matching rows
    [tmp, i  ] = sortrows(nodes);
    [tmp, num] = rlencode(tmp,1); %#ok
    % Check that at most two faces share the same set of nodes
    assert(max(num) <= 2, ['This function currently does not support ', ...
         'the case when more than two faces share the same set of nodes']);
    fix = rldecode(num==2, num, 1);
    % Find correspoing face numbers
    faces = (1:G.faces.num)'; faces = faces(i); faces = faces(fix);
    % Reshape to adjacency list
    N = reshape(faces, 2, [])';
    
end