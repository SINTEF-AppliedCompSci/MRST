function dual = partitionUIdual(CG, blockSizes)
% Create coarse dual partition/grid
%
% SYNOPSIS:
%    Create coarse dual partition/grid using logical partitioning
% DESCRIPTION:
%
%
% REQUIRED PARAMETERS:
%    CG -- Coarse grid as defined by partitionUI
%    blockSizes -- block sizes for partitionUI
% RETURNS:
%    dual -- Categorization of the nodes
% NOTE:
%    Use dualPartition for non-Cartesian grids or grids with erosion

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


if numel(blockSizes) == 2
    blockSizes(3) = 1;
end
G = CG.parent;
n = G.cells.num;

if G.griddim == 3 && G.cartDims(3) > 1
    dimension = 3;
else
    dimension = 2;
end

% Decrement all positions because working with zero indexing is easier
% to work with when partitioning...
cells = G.cells.indexMap-1;
% Find positions of all the nodes in ijk space
spaces = cell(dimension,1);
[spaces{:}] = ind2sub(G.cartDims, G.cells.indexMap);
uniques = cell(dimension,1);


for d = 1:dimension
    ms = min(spaces{d});
    Ms = max(spaces{d});
    % Find the number of positions in the current dimension
    M = Ms - ms + 1;
    % Do a load balanced distribution of the positions in the same manner
    % as in the primal grid. Cast to double to avoid integer division
    % problems.
    balanced = lbLinDist(double(0:(Ms-ms)), double(M), double(blockSizes(d)));
    u = unique(balanced);
    uniques{d} = zeros(numel(u),1);
    index = 1;
    for i = 1:numel(u)
        % Keep a running index and find the midpoints of the intervals
        % to use as centerpoints for the dual grid
        Ni = sum(balanced == u(i));
        uniques{d}(i) = index + round(Ni/2) - 1;
        index = index + Ni;
    end
end

if ~any(strcmp(G.type,'cartGrid'))
    warning('warn:logicalcart','Logical partitioning on non-Cartesian grid!');
end
% Gather all the different positions in each D where the grid has dual edges
tmp = zeros(n,dimension);
for d = 1:dimension
    ii = uniques{d};
    for i = 1:numel(ii)
        tmp(:,d) = tmp(:,d) | spaces{d} == ii(i);
    end
end
% The central nodes are those which are in two sets of edges in 2D, and 3
% sets of edges in 3D. Increment to get one indexing.
dual.nn =  cells(sum(tmp,2)==dimension)+1;
% Any node which exists in more than one edge must be filtered to
% decouple the edge system. Increment to get original indexing.
dual.lineedge = cells(sum(tmp,2)>1)+1';
% Do an OR on each dimension to find all nodes corresponding to the
% midpoints of each primal coarse block
if dimension == 3
    tmp = tmp(:,1) | tmp(:,2) | tmp(:,3);
else
    tmp = tmp(:,1) | tmp(:,2);
end

% Find the indices of nodes on the edge, increment to get back to one indexing
dual.ee = cells(tmp == 1)+1;

% Use the new behavior of setdiff to ensure forward compatability
% Ensure that the different edges are distinct

dual.lineedge = setdiff(dual.lineedge, dual.nn);
dual.ee = setdiff(setdiff(dual.ee, dual.lineedge), dual.nn);
dual.ii = setdiff(G.cells.indexMap, [dual.ee; dual.nn; dual.lineedge]);


function f = lbLinDist(f, M, B)
% lbLinDist -- Load-balanced linear distribution
% See Eric F. Van de Velde, Concurrent Scientific Computing,
% 1994, Springer Verlag, p. 54 (Sect. 2.3) for details.
%
% Maps index set (0..M-1) to blocks (0..B-1).

L = floor(M ./ B);  % Tentative number of cells per coarse block.
R = mod(M, B);      % Additional cells not previously accounted for.
f = max(floor(f ./ (L + 1)), floor((f - R) ./ L));
