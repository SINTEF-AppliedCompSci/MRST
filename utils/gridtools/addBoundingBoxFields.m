function G = addBoundingBoxFields(G, varargin)
% Add minimal bounding box for each cell/face
% Useful for quickly picking subset of candidate cells/faces before more 
% costly geometry is invoked
%
% SEE ALSO:
%  computeVerticalGridIntersection, computeTraversedCells,
%  addBoundingBoxFields

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

opt = struct('cells', false, ...
             'faces', true);
opt = merge_options(opt, varargin{:});

gdim = size(G.nodes.coords, 2);
% cells (cellNodes fails for e.g., Olympus)
if opt.cells
    cno  = cellNodes(G);
    nn   = accumarray(cno(:,1), ones(size(cno,1), 1));
    npos = cumsum([1;nn]);
    % max and min ix after sorting
    minIx = npos(1:end-1);
    maxIx = npos(2:end)-1;
    bbox = nan(G.cells.num, G.griddim);
    for d = 1:gdim
        tmp = sortrows([cno(:,1), G.nodes.coords(cno(:,3),d)]);
        bbox(:,d) = tmp(maxIx, 2) - tmp(minIx, 2);
    end
    G.cells.bbox = bbox;
end

% faces
if opt.faces
    npos = G.faces.nodePos;
    fno  = rldecode((1:G.faces.num)', diff(npos));
    % max and min ix after sorting
    minIx = npos(1:end-1);
    maxIx = npos(2:end)-1;
    bbox = nan(G.faces.num, G.griddim);
    for d = 1:gdim
        tmp = sortrows([fno(:,1), G.nodes.coords(G.faces.nodes, d)]);
        bbox(:,d) = tmp(maxIx, 2) - tmp(minIx, 2);
    end
    G.faces.bbox = bbox;
end
end
