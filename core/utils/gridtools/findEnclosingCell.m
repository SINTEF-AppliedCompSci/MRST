function [cc, ii, jj] = findEnclosingCell(G, pt, cells)
% Find cell(s) containing points. The cells must be convex.
%
% SYNOPSIS:
%   c = findEnclosingCell(G, pt)
%   c = findEnclosingCell(G, pt, cells)
%
% PARAMETERS:
%   G  - Valid grid structure with normals, face and cell centroids. The
%        cells must be convex.
%
%   pt - Set of point coordinates, represented as an n-by-dim array.
%
% OPTIONAL PARAMETERS
%   cells - Subset of grid cells, either as cell indices or a logical mask.
%           Default is cells = 1:G.cells.num
%
% RETURNS:
%   cc - Set of grid cells. Specifically, c(k), k=1:n, contains pt(k,:). If
%        a point lies on the boundary between two cells, the function
%        returns the smallest index. c(k) is zero if the point is not
%        found.
%   ii - Indices of the points in pt that are found in the grid cells.
%   jj - Indices of the cells in the input cells that contain the points in pt.
%
% EXAMPLE:
%
%   G = cartGrid([2,3,4]);
%   G = computeGeometry(G);
%   pt = [[0.5,0.5,0.5];    % inside cell 1
%         [1.5, 2.5, 3.5]]; % inside cell 24
%   cells = [2:24];         % Don't include cell 1
%   [cc,ii,jj] = findEnclosingCell2(G, pt, cells);
%
%   Then
%   cc = [0, because point 1 is in cell 1, which is not in the set of cells
%          24], because point 2 is in cell 24
%
%   ii = [0, because point 1 is not in the set of cells
%         2], because point 2 is in cell 24
%
%   jj = [0, because point 1 is not in the set of cells
%         23], because point 2 is in cells(23)
%
% NOTE:
%   The function can be very slow for large grids and many points. It is
%   recommended to reduce the grid, for example by using
%   Gsub = extractSubgrid(G, cells)
%   and use the field Gsub.cells.global to map the cells to the original grid.
%
% SEE ALSO:
%   `pebi`.

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

% Check cell input
if nargin < 3,       cells = (1:G.cells.num)';  end
if islogical(cells), cells = find(cells);       end
assert(min(cells) >= 1 && max(cells) <= G.cells.num, ...
       'Cell subset is ouside [1,G.cells.num]');
cells = cells(:);

% Get the normals for each face in each cell
% Map for expanding faces to cell faces
fno = mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1);
cfmap = G.cells.faces(fno, 1);

% Map normals and face centroids to cell face numbering
if ~isfield(G.faces, 'normals')
    error('Grid structure does not contain normals. Please run eg. ''computeGeometry''.');
end
n = G.faces.normals(cfmap, :);
cfc = G.faces.centroids(cfmap, :);

% Expand cell centroids
ncf = diff(G.cells.facePos);
ncf = ncf(cells);
j = rldecode(cells, ncf);
cc = G.cells.centroids(j, :);

% Compute sign of dot product
sgn = sign(dot(cfc-cc, n, 2));

% Adjust sign to make n point in the same direction as cfc-cc
n = bsxfun(@times, sgn, n);

% Expand cell numbers
cellno = rldecode(cells, ncf);

% Initialize
cc = zeros(size(pt, 1), 1);
ii = zeros(size(pt, 1), 1);
jj = zeros(size(pt, 1), 1);

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
        cc(k) = min(ind);
        ii(k) = k;
        jj(k) = find(cells == cc(k));
    end
end

end
