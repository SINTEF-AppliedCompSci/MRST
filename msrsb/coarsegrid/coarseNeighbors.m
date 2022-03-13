function [cellNeighbors, coarsefaceNeighbors] = coarseNeighbors(CG, cell, isMPFA, faceNodePos)
% Determine coarse neighbors for a given coarse cell (MPFA or TPFA-like).

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
if nargin < 3
    isMPFA = false;
end
cellNeighbors = getCellNeighbors(CG, cell);
if isMPFA
    G = CG.parent;
    if nargin < 4
        faceNodePos = rldecode(1 : G.faces.num, diff(G.faces.nodePos), 2) .';
    end
    % Fine the boundary of the coarse tpfa-like region
    fineInside = find(ismember(CG.partition, cellNeighbors));
    bf = boundaryFaces(G, fineInside);

    % Find those nodes
    bfNodes = G.faces.nodes(mcolon(G.faces.nodePos(bf), G.faces.nodePos(bf+1)-1));


    % Find coarse face nodes at the fine level
    coarsefaces = gridCellFaces(CG, cell);
    fa = CG.faces.fconn(mcolon(CG.faces.connPos(coarsefaces), CG.faces.connPos(coarsefaces+1)-1));
    nodes = gridFaceNodes(G, fa);

    % Find interfaces and do the neighborhood song and dance because
    % this is MATLAB
    commonNodes = intersect(bfNodes, nodes);
    tmp = false(G.nodes.num, 1);
    tmp(commonNodes) = true;
    commonNodePos = tmp(G.faces.nodes);

    faceNeighbors = unique(faceNodePos(commonNodePos));

    fa = G.faces.neighbors(faceNeighbors, :);
    fa = fa(:);
    fa = fa(fa ~= 0);
    partition = CG.partition(fa);

    cellNeighbors = unique([partition; cellNeighbors]);
    cellNeighbors = cellNeighbors(cellNeighbors~=0 & cellNeighbors~=cell);  

    cfno = rldecode((1:CG.faces.num)', diff(CG.faces.connPos), 1);
    coarsefaceNeighbors = unique(cfno(ismember(CG.faces.fconn, faceNeighbors)));
else
    N = CG.faces.neighbors;
    coarsefaceNeighbors = find((ismember(N(:,1), cellNeighbors) & N(:,2) == cell) | ...
                               (ismember(N(:,2), cellNeighbors) & N(:,1) == cell));
end
end
% find(ismember(CG.faces.fconn, faceNeighbors))
