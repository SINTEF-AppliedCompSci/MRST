function [cell, largestVol, allFaces, point, traps] = findOptimalInjectionPoint(G, res)
%Find the optimal point to inject CO2
%
% SYNOPSIS:
%   [cell, largestVol, allFaces, point] = findOptimalInjectionPoint(G, res)
%
% DESCRIPTION:
%   Finds the best possible point to inject CO2. It is done by assuming
%   that the best injection point is at the interface between two trap
%   regions which together correspond to the largest trap trees.
%
% REQUIRED PARAMETERS:
%   G   - Top surface grid.
%   res - trapAnalysis output.
%
% OPTIONAL PARAMETERS:
%   None.
%
% RETURNS:
%   cell - Global cell index of best possible injection point.
%
%   largestVol - Structural trapping volume reachable from this point.
%
%   allFaces   - The fine faces corresponding to the interface between the
%                two trapping regions of the optimal point.
%
%   point      - The centroid of cell.

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

    [trees, trapVols] = maximizeTrapping(G, 'removeOverlap', false,...
                                            'res', res, ...
                                            'calculateAll', true);
    % Sort the trees on startpoint so that trees(i) corresponds to trap i
    [v, sortInd] = sort([trees.root]); %#ok
    trees = trees(sortInd);
    
    N = G.faces.neighbors;
    
    % Find the region at each side of each fine face
    tmp = [0; res.trap_regions];
    region = sort(tmp(N+1), 2);
    % Find all interface indices that boundry to two different regions
    candidates = find(region(:, 1) ~= region(:, 2));
    % Map out the UNIQUE pairs, and index into candidates so that
    % candidates remain a map into g.faces.neighbors.
    [v, i, j] = unique(region(candidates, :), 'rows');%#ok
    candidates = candidates(i);
    
    % Find the trap-pair corresponding to the largest pair of injectors
    faceVols = zeros(size(candidates, 1), 1);
    for i = 1:size(candidates, 1)
        traps = region(candidates(i), :);
        % Ignore zero outside
        traps = traps(traps~=0);
        % Strip overlap between trees
        downstream = unique(vertcat(trees(traps).traps));
        faceVols(i) = sum(trapVols(downstream));
    end
    [largestVol, largestInterface] = max(faceVols);
    
    % We now know the neighborship of the optimal face - evaluate all these
    % faces to find the lowest possible cell in the region. This is the
    % optimal place of injection
    n = region(candidates(largestInterface), :);
    allFaces = find(region(:, 1) == n(1) & region(:, 2) == n(2));
    cells = N(allFaces, :);
    cells = unique(cells(:));
    cells = cells(cells~=0);
    
    [v, ind] = min(G.cells.z(cells));%#ok
    cell = cells(ind);
    point = G.cells.centroids(cell, :);
    traps = n;
end
