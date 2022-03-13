function [cellcenters, facecenters] = mapCenters(CG, blockCenters, faceCenters)
% Map fine cells as centers of coarse cells and faces.

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
        faceCenters = CG.faces.centroids;
        if nargin < 2
            blockCenters = CG.cells.centroids;
        end
    end
    
    
    G = CG.parent;
    cellcenters = zeros(CG.cells.num, 1);
    for i = 1:CG.cells.num
        local = find(CG.partition == i);
        cellcenters(i) = doMap(i, local, blockCenters, G.cells.centroids);
    end
    
    if nargout > 1
        facecenters = zeros(CG.faces.num, 1);
        for i = 1:CG.faces.num
            local = CG.faces.fconn(CG.faces.connPos(i):CG.faces.connPos(i+1)-1);
            facecenters(i) = doMap(i, local, faceCenters, G.faces.centroids);
        end
    end
end


function l = doMap(coarseind, localfine, coarsepts, finepts)
    % Add a tiny bit of bias to the selection so that we avoid
    % consistent mappings on cartesian grids where the center is exactly
    % between two cells
    dist = bsxfun(@(x, y) x-y, finepts(localfine,:), ...
                               coarsepts(coarseind,:) - 10*sqrt(eps));
    [v, ind] = min(sqrt(sum(dist.^2, 2))); %#ok
    l = localfine(ind);
end
