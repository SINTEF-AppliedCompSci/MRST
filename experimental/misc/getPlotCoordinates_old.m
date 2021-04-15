function coords = getPlotCoordinates_old(G, varargin)
%Undocumented Utility Function

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

    opt = struct('n'      , 100  , ...
                 'phaseNo', 1    , ...
                 'plot'   , true , ...
                 'plot1d' , false);
    
    [opt, ~] = merge_options(opt, varargin{:});

    xmax = max(G.nodes.coords, [], 1);
    xmin = min(G.nodes.coords, [], 1);

    dx   = xmax - xmin;
    n = round(opt.n.*dx./max(dx));
    
    x = cell(1,G.griddim);
    if opt.plot1d
        x{1} = linspace(xmin(1), xmax(1), n(1))';
        for d = 2:G.griddim
            x{d} = (xmin(d) + xmax(d))/2;
            n(d) = 1;
        end
    else
        for d = 1:G.griddim
            x{d} = linspace(xmin(d), xmax(d), n(d))';
        end
    end

    x0 = cell(size(x));
    [x0{:}] = ndgrid(x{:});

    nPts = numel(x0{1});
    points = zeros(nPts, G.griddim);
    for d = 1:G.griddim
        points(:,d) = x0{d}(:);
    end
    
    cells = ones(nPts,1);
    keep = false(nPts,1);
    for cNo = 1:G.cells.num
        if G.griddim == 2
            nodes = getCellNodes(cNo, G);
        else
            ix = G.cells.nodePos(cNo):G.cells.nodePos(cNo+1)-1;
            nodes = G.cells.nodes(ix);
        end
        xn       = G.nodes.coords(nodes,:);
        [in, on] = deal(false(numel(keep),1));
        [in(~keep), on(~keep)] = inCell(G, points(~keep,:), xn);
%         [in, on] = inpolygon(points(:,1), points(:,2), xn(:,1), xn(:,2));
        cells(in | on) = cNo;
        keep = keep | in | on;
    end
    coords = struct('points', points, 'cells', cells, 'keep', keep, 'n', n);
    
end

function nodes = getCellNodes(cell, G)
    
    faces = G.cells.faces(G.cells.facePos(cell):G.cells.facePos(cell+1)-1);
    nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
    nodes = reshape(nodes, 2, [])';
    swap = G.faces.neighbors(faces,1) ~= cell;
    nodes(swap,:) = nodes(swap, [2,1]); nodes = nodes(:,1);
%     nodes = unique(nodes);
    
end

function [in, on] = inCell(G, x, xv)
    on = false(size(x,1),1);
    if G.griddim == 2
        [in, on] = inpolygon(x(:,1), x(:,2), xv(:,1), xv(:,2));
    else
        tri = delaunay(xv);
        tn  = tsearchn(xv, tri, x);
        in  = ~isnan(tn);
    end
end
