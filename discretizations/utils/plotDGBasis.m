function h = plotDGBasis(G, basis, cellNo, varargin)
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

    opt          = struct('n', 15);
    [opt, extra] = merge_options(opt, varargin{:});
    
    nodes = G.cells.nodes(G.cells.nodePos(cellNo):G.cells.nodePos(cellNo+1)-1);
    xn    = G.nodes.coords(nodes,:);
    xc = G.cells.centroids(cellNo,:);
    
    xmin = min(xn,[],1);
    xmax = max(xn,[],1);
    
    dx   = max(xmax - xc, xc - xmin)*2;
    
    g = computeGeometry(pebiGrid2D(max(dx)/opt.n, dx, 'polyBdr', xn));
    x = g.nodes.coords;
    x = (x - G.cells.centroids(cellNo,:))./(dx/2);
    
    z = cell(1, basis.nDof);
    
    xlim = round([xc(1) - dx(1)/2, xc(1) + dx(1)/2, xc(2) - dx(2)/2, xc(2) + dx(2)/2, -1.1, 1.1]*100)/100;
    
    
    faces = boundaryFaces(g);
    for p = 1:basis.nDof
        h(p) = figure();
        plotFaces(g, faces);
        z{p} = basis.psi{p}(x);
        unstructuredCarpetPlot(g, z{p}, 'nodeval', true, extra{:});
        view(3);
        axis tight
        zlim([-1.1, 1.1]);
        pbaspect([1,1,.5]);
        ax = gca;
        ax.XTick = round(linspace(xlim(1), xlim(2), 5)*100)/100;
        ax.YTick = round(linspace(xlim(3), xlim(4), 5)*100)/100;
        caxis([-1,1]);
        colormap('pink')
        box on
        ax.ZDir = 'normal';
    end
    
end
