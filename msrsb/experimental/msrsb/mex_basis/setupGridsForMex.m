function g = setupGridsForMex(A, CG, varargin)
%Undocumented Utility Function

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

    A = A - diag(sum(A, 2));
    d = diag(A);
    
    
    bad = find(d == 0);
    assert(numel(bad) == 0, 'Zero diagonal entries present!');
    % A = A + sparse(bad, bad, max(mean(d), 1e-8), size(A, 1), size(A, 2));
    
    N = CG.cells.num;
    G = CG.parent;
    
    [subsys, cells, isBnd] = deal(cell(N, 1));
    
    globalBnd = false(G.cells.num, 1);
    
    globalBndToLocal = cell(G.cells.num, 1);
    
    
    outerInd = 0;
    for i = 1:N
        interaction = CG.cells.interaction{i};

        
        isInt = false(G.cells.num, 1);
        locPart = CG.partition(interaction);
        tmp = bndNeighbors(G, find(CG.partition == i));
        
        % Local interaction region partition: The set of all coarse blocks
        % that are either direct two-point neighbors to the current coarse
        % block, or is covered by the interaction region somehow.
        locPart = unique(vertcat(locPart, CG.partition(tmp)));

        for j = 1:numel(locPart)
            isInt(CG.partition == locPart(j)) = true;
        end
        % Pad with zero for boundary
        isInt = [false; isInt];
        
        [c_self, c_other] = boundaryNeighbors(G, interaction);
        
        [isCell, isBndLoc] = deal(false(G.cells.num, 1));
        % If the neighbor over a face is outside our local collection of
        % coarse blocks, we instead convert the current interaction cell to
        % a boundary cell. This is done to avoid strange renormalizations.
        isBndLoc(c_other( isInt(c_other+1))) = true;
        isBndLoc(c_self( ~isInt(c_other+1))) = true;
        isCell(interaction) = true;

        c = find(isCell | isBndLoc);
        c = sort(c);
        % Local boundary indicator
        isBnd{i}  = ismember(c, find(isBndLoc));
        % Matrix subsystem in local indices
        subsys{i} = A(c, c);
        
        % Local to global mapping of these cells (zero indexing!)
        cells{i} = c - 1;
        
       
        globalBnd(c(isBnd{i})) = true;
        outerInd = outerInd + numel(c);
    end
    
    tmp = vertcat(cells{:});
    
    for i = 1:numel(tmp);
        if globalBnd(tmp(i)+1)
            globalBndToLocal{tmp(i)+1} = [globalBndToLocal{tmp(i)+1}; i-1];
        end
    end
    
    bndIndicator = cell(N, 1); 
    for i = 1:N
        subs = ~isBnd{i} & globalBnd(cells{i}+1);
        
        bndIndicator{i} = zeros(numel(cells{i}), 1);
        bndIndicator{i}(subs) = 2;
        bndIndicator{i}(isBnd{i}) = 1;
    end
    isBnd = bndIndicator;
        
    g.subsys = subsys;
    g.cells = cells;
    g.isBnd = isBnd;
    
    g.nc = CG.cells.num;
    g.nf =  G.cells.num;
    
    g.partition = CG.partition;
    
    g.globalBoundary = find(globalBnd) - 1;
    g.globalBoundaryToLocal = globalBndToLocal(g.globalBoundary+1);
    
    g.centers = CG.cells.centers - 1;
end

function c = bndNeighbors(g, subs)
    fa = gridCellFaces(g, subs);
    
    c = g.faces.neighbors(fa, :);
    c = unique(c(:));
    c = c(c ~= 0);
end


function [c_self, c_other] = boundaryNeighbors(g, subs)
    fa = boundaryFaces(g, subs);
    
    N = g.faces.neighbors(fa, :);
    
    bad = any(N == 0, 2);
    fa = fa(~bad);
    N = N(~bad, :);
    
    isInt = ismember(N(:, 1), subs);
    
    [c_self, c_other] = deal(zeros(numel(fa), 1));
    
    c_self( isInt) = N( isInt, 1);
    c_self(~isInt) = N(~isInt, 2);
    
    c_other( isInt) = N( isInt, 2);
    c_other(~isInt) = N(~isInt, 1);

end
