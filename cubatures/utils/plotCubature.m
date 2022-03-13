function plotCubature(G, cubature, elements, varargin)
%Undocumented Utility Function

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

    opt = struct('plotBoundingBox', true, 'plotTri', true, 'smin', 5, 'smax', 20, 'markerStyle', 'ok');
    [opt, extra] = merge_options(opt, varargin{:});
    
    if G.griddim - cubature.dim == 1
        faces = elements;
        cells = reshape(G.faces.neighbors(faces,:), [], 1);
        [~, x, w, cellNo, faceNo] = cubature.getCubature(faces, 'face');
    else
        cells = elements;
        [~, x, w, cellNo] = cubature.getCubature(cells, 'cell');
    end

    hold on
    plotGrid(G, cells, extra{:}); axis off
    if opt.plotBoundingBox
        for i = 1:numel(cells)
            rectangle('Position', [G.cells.centroids(cells(i),:) - G.cells.dx(cells(i),:)/2, G.cells.dx(cells(i),:);]);
        end
    end
    
    if any(strcmpi(class(cubature), {'TriangleCubature', 'TetrahedronCubature'})) && opt.plotTri
        tri = cubature.triangulation.triPos(cellNo):cubature.triangulation.triPos(cellNo+1)-1;
        if G.griddim == 2
            trimesh(cubature.triangulation.ConnectivityList(tri,:), ...
                    cubature.triangulation.Points(:,1), ...
                    cubature.triangulation.Points(:,2), ...
                    'linestyle', '-', 'color', 'k')
        else
            trimesh(cubature.triangulation.ConnectivityList(tri,:), ...
                    cubature.triangulation.Points(:,1), ...
                    cubature.triangulation.Points(:,2), ...
                    cubature.triangulation.Points(:,3), ...
                    'linestyle', '-', 'edgecolor', 'k')
        end
    end

    alpha = (w - min(w))/(max(w) - min(w));
    alpha(isnan(alpha)) = 1;
    if isempty(opt.smax)
        opt.smax = opt.smin*max(w)/min(w);
    end
    scaling = opt.smin*(1-alpha) + opt.smax*alpha;
    if G.griddim == 2
        plotPoint = @(x, varargin) plot(x(:,1), x(:,2), varargin{:});
    else
        plotPoint = @(x, varargin) plot3(x(:,1), x(:,2), x(:,3), varargin{:});
    end
    for i = 1:numel(w)
        plotPoint(x(i,:), opt.markerStyle, 'markerSize', scaling(i), 'markerFaceColor', 'k')
%         plot(x(i,1), x(i,2), opt.markerStyle, 'markerSize', scaling(i), 'markerFaceColor', 'k')
    end

end
