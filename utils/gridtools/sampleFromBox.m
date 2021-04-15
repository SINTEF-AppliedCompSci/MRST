function q = sampleFromBox(G, p, c)
%Sample from data on a uniform Cartesian grid that covers the bounding box
%
% SYNOPSIS:
%    q = sampleFromBox(G, p)
%    q = sampleFromBox(G, p, c)
%
%
% PARAMETERS:
%    G  - grid structure
%    p  - input array to sample from. The array is assumed to be values
%         on a uniform Cartesian grid covers the bounding box of grid `G`.
%         The values should be given as a Nx by Ny matrix where Nx and Ny
%         is the number of values in the x and y direction, respectively.
%    c  - array specifying a subset of cells in `G` in which values will be
%         sampled
%
% RETURNS:
%    q  - output vector sampled at the centroids of grid cells in `G`

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


assert(any(G.griddim==[2,3]), 'Number of space dimensions must be 2 or 3');

% Check number of dimensions, using a hack to avoid problems if the last
% the size of the array should have been (n,m,1)
nd = size(p);
if numel(nd)==G.griddim-1
   nd = [nd 1];
end
assert(G.griddim == numel(nd), 'Number of space dimensions do not match');

% Process input parameters
if nargin<3
   c = 1:G.cells.num;
   X = G.cells.centroids;
else
   X = G.cells.centroids(c,:);
end

% -- Scale coordinates
% First we find the cells with largest/smallest centroid values in each
% spatial direction.
[buf,i] = min(X,[],1); %#ok<*ASGLU>
[buf,j] = max(X,[],1);

% The bounding box is then found among the nodes of these cells and mapped
% onto the box that contains the sample data.
n  = unique(gridCellNodes(G, c([i j]), 'unique', false));
Y  = G.nodes.coords(n,:);
mC = min(Y);
coords = bsxfun(@minus, X, mC);
coords = bsxfun(@rdivide, coords, max(Y)-mC);
ijk = num2cell(ceil(bsxfun(@times, coords, nd)),1);

% Sample
q = p(sub2ind(nd,ijk{:}));
end
