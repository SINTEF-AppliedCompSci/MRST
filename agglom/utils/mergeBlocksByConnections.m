function p = mergeBlocksByConnections(G, p, T, minBlockSize)
%Merge blocks based on connection strength
%
% SYNOPSIS:
%   p = mergeBlocksByConnections(G, partition, T, minBlockSize)
%
% DESCRIPTION:
%   An alternate merging algorithm based on face connections.
%
% REQUIRED PARAMETERS:
%   G    - Grid structure
%
%   p    - Partition vector with one entry per cell in G.
%
%   T    - Connection strength. One entry per face in G, including boundary
%          faces. Negative weights are interpreted as zero connections.
%
%   mz   - Minimum block size. The algorithm terminates when all coarse
%          blocks have at least mz fine cells.
%
% RETURNS:
%   p    - Modified partition vector.
%
% SEE ALSO:
%   `mergeBlocks`, `mergeBlocks2`

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

% Ignore boundary connections
T(any(G.faces.neighbors == 0, 2)) = -1;

ok = false;
bad = getBlockCount(p) < minBlockSize;
dispif(mrstVerbose(), 'Initially: %d blocks are under limit...\n', sum(bad))
while ~ok
    % Set up coarse grid
    CG = generateCoarseGrid(G, p);
    p_next = p;
    % Get coarse connections and find coarse connection strength as the
    % sum of all fine scale connections in that face.
    coarseFaceNo = rldecode(1:CG.faces.num, diff(CG.faces.connPos), 2)';
    Tc = accumarray(coarseFaceNo, T(CG.faces.fconn));

    % Targets are all blocks with less than the desired number of
    % cells.
    targets = find(bad);
    disconnected = 0;
    for i = 1:numel(targets)
        C = targets(i);
        fa = gridCellFaces(CG, C);
        fa = fa(~(Tc(fa)<0));
        if isempty(fa)
            disconnected = disconnected+1; continue;
        end
        % Find strongest connection for this local block, and merge
        % into that neighbor.
        [v, strongest] = max(Tc(fa));                                      %#ok<ASGLU>

        % Set partition to that of the strongest neighbor
        other = CG.faces.neighbors(fa(strongest), :);
        other = other(other ~= C);
        p_next(p_next == C) = other;
    end
    % Ensure ordering goes 1, ..., n
    p = compressPartition(p_next);
    % Update bad entries
    bad = getBlockCount(p) < minBlockSize;
    if disconnected == sum(bad)
        dispif(mrstVerbose(), ...
            '%d disconnected blocks remain under limit...\n', sum(bad))
        ok = 1;
    else
        ok = ~any(bad);
        dispif(mrstVerbose(), '%d blocks remain under limit...\n', sum(bad))
    end
end
end

function n = getBlockCount(p)
    n = accumarray(p, 1);
end
