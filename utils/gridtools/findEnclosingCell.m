function c = findEnclosingCell(G, pt)
% Find cell(s) containing points. The cells must be convex.
%
% SYNOPSIS:
%   c = findEnclosingCell(G, pt)
%
% PARAMETERS:
%   G  - Valid grid structure with normals, face and cell centroids. The 
%        cells must be convex.
%
%   pt - Set of point coordinates, represented as an n-by-dim array.
%
% RETURNS:
%   c  - Set of grid cells. Specifically, c(k), k=1:n, contains pt(k,:). If
%        a point lies on the boundary between two cells, the function 
%        returns the smallest index. c(k) is zero if the point is not
%        found.
%
% SEE ALSO:
%   `pebi`.

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

% Get the normals for each face in each cell
cells = (1:G.cells.num)';

% Map for expanding faces to cell faces
cfmap = G.cells.faces(:, 1);

% Map normals and face centroids to cell face numbering
n = G.faces.normals(cfmap, :);
cfc = G.faces.centroids(cfmap, :);

% Expand cell centroids
j = rldecode(cells, diff(G.cells.facePos));
cc = G.cells.centroids(j, :);

% Compute sign of dot product
sgn = sign(dot(cfc-cc, n, 2));

% Adjust sign to make n point in the same direction as cfc-cc
n = bsxfun(@times, sgn, n);

% Expand cell numbers 
cellno = rldecode(cells, diff(G.cells.facePos));

% Initialize
c = zeros(size(pt, 1), 1);

for k = 1:size(pt, 1)
    % Calculate vector between face centroids and pt
    b = bsxfun(@minus, cfc, pt(k, :));
    
    % Check if n and b are pointing in the same direction
    v = dot(n, b, 2);
    V = v >= 0.0;
    
    % Find cell indices with all non-negative dot products
    i = accumarray(cellno, V, [G.cells.num, 1], @(x) all(x));
    ind = find(i);
    if ~isempty(ind)
        c(k) = min(ind);
    end
end
