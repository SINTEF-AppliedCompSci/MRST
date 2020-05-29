function out = getMexDiscreteDivergenceJacPrecomputes(model)
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

    G = model.G;
    N = model.operators.N;
    if 0
        % Use grid
        interior = false(G.faces.num, 1);
        interior(model.operators.internalConn) = true;

        facePos = G.cells.facePos;
        cellNo = rldecode(1 : G.cells.num, diff(facePos), 2) .';
        faces = G.cells.faces(:, 1);

        active = interior(faces);

        cellNo = cellNo(active);
        faces = faces(active);

        internal = model.operators.internalConn;
        remap = zeros(G.faces.num, 1);
        remap(internal) = (1:sum(internal))';
        faces = remap(faces);
    else
        % Use operators.N
        f = (1:size(N, 1))';
        pos = sortrows([[N(:, 1); N(:, 2)], [f; f]]); % Sort cell -> face maps
        faces = pos(:, 2);
        cellNo = pos(:, 1);
    end
    
    lc = N(faces, 1);
    rc = N(faces, 2);
    left = (lc ~= cellNo);
    cells = left.*lc + ~left.*rc;
    cells = cells - 1;
    cells(~left) = -cells(~left);



    [c, ps] = rlencode(cellNo);
    if numel(ps) ~= G.cells.num
        % we have some neighborless cells, insert zeros in ps accordingly
        cix = false(G.cells.num, 1);
        cix(cellNo) = true;
        tmp = ps;
        ps = zeros(G.cells.num,1);
        ps(cix) = tmp;
    end
        
    facePos = cumsum([1; ps]);

    localCellIndex = zeros(G.cells.num, 1);
    for i = 1:G.cells.num
        loc = facePos(i):facePos(i+1)-1;
        ca = cells(loc);
        [~, ix] = sort(abs(ca));

        tmp = faces(loc);
        faces(loc) = tmp(ix);
        cells(loc) = ca(ix);

        tmp = find(abs(ca(ix))+1 > i, 1, 'first')-1;
        if isempty(tmp)
            tmp = numel(ca);
        end
        localCellIndex(i) = tmp;

    end
    % Need to sort cells as well
    fpos0 = facePos-1;
    faces = faces - 1;
    out = struct('facePos', fpos0, 'faces', faces, 'cells', cells, 'cellIndex', localCellIndex);
end
