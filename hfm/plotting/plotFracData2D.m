function varargout = plotFracData2D(G, data, varargin)
% plotFracData plots data inside fracture cells only, given a grid with
% matrix and fractures. The function is designed for 2D grids.
%
% SYNOPSIS:
%       plotFracData2D(G, data)
%       plotFracData2D(G, data, 'pn1', pv1, ...)
%   h = plotFracData2D(...)
%
% REQUIRED PARAMETERS:
%
%   G    - Grid data structure with fractures as defined by
%          assembleGlobalGrid.
%
%   data - data to plot with values inside fracture cells only. Matrix data
%          is ignored.
%
% OPTIONAL PARAMETERS:
%
%   wide        - This option can be used to redefine a fracture grid (by
%                 recalling FracTensorGrid2D) with larger aperture for
%                 plotting purposes. Can be thought of as an
%                 exaggerated/magnified fracture grid.
%
%   width       - fracture width/aperture in m if key 'wide' is used.
%
%   F           - Structure containing information about partitioned
%                 fracture lines as returned by assembleFracNodes2D. See
%                 also gridFracture.
%
%   CG          - Fracture coarse grid as returned by getRsbGrids_HFM.
%
%   cmap        - Colormap desired. default = jet.
%
%   outline     - If true, matrix coarse grid will be outlined.
%
%   outlineCoarseNodes - If true, will outline fracture coarse nodes.
%
%   plotMatrix  - If true, plots matrix grid as well.
%
% RETURNS:
%   h - Handle to polygonal patch structure as defined by function
%       plotFaces.  OPTIONAL.
%
% SEE ALSO:
%   outlineCoarseGrid, plotToolbar, patch

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


opt = struct(...
    'wide'               ,false  ,...
    'width'              , 2     ,...
    'F'                  ,[]     ,...
    'CG'                 ,[]     ,...
    'cmap'               ,jet    ,...
    'outline'            ,false  ,...
    'outlineCoarseNodes' ,true   ,...
    'plotMatrix'         ,false  );

[opt,patchopt] = merge_options(opt, varargin{:});
if opt.wide
    assert(~isempty(opt.F) && ~isempty(opt.CG),'Provide ''F'' and ''CG'' struct');
end
if numel(data) == G.cells.num, data = data(G.Matrix.cells.num+1:G.cells.num); end
if opt.wide
    Gfwide = FracTensorGrid2D(G.Matrix,opt.F,opt.width); 
    Gfwide.nnc = G.nnc; Gfwide.Matrix = G.Matrix; Gfwide = assembleFracGrid(Gfwide);
    plotToolbar(Gfwide, data, 'EdgeAlpha', 0.1, patchopt{:}); hold on
    if opt.outline, outlineCoarseGrid(Gfwide,opt.CG.partition); end, axis tight off; c = colormap(opt.cmap);
    if opt.plotMatrix, plotGrid(G.Matrix,'FaceColor','none','EdgeAlpha',0.03); hold on; end
    caxis([0 1]); c(1,:) = [1 1 1]; colormap(c); colorbar
    if opt.outlineCoarseNodes
    for i = 1:numel(opt.CG.cells.centers)
        ci = opt.CG.cells.centers(i);
        x = Gfwide.cells.centroids(ci,1);
        y = Gfwide.cells.centroids(ci,2);
        faces = Gfwide.cells.faces(Gfwide.cells.facePos(ci):Gfwide.cells.facePos(ci+1)-1,1);
        cnodes = [];
        for k = 1:numel(faces)
            fnodes = Gfwide.faces.nodes(Gfwide.faces.nodePos(faces(k)):Gfwide.faces.nodePos(faces(k)+1)-1);
            cnodes = [cnodes;fnodes]; %#ok
        end
        x = [x;Gfwide.faces.centroids(faces,1);Gfwide.nodes.coords(cnodes,1)]; %#ok
        y = [y;Gfwide.faces.centroids(faces,2);Gfwide.nodes.coords(cnodes,2)]; %#ok
        K = convhull(x,y);
        plot(x(K),y(K),'m','LineWidth',2.5);
%         text(x+0.1,y,num2str(i));
    end
    end
else
    Gf = assembleFracGrid(G);
    plotToolbar(Gf, data, 'EdgeAlpha', 0.05, patchopt{:}); hold on
    if opt.plotMatrix, plotGrid(G.Matrix,'FaceColor','none','EdgeAlpha',0.05); end 
    axis tight off; c = colormap(opt.cmap);
    caxis([0 1]); c(1,:) = [1 1 1]; colormap(c);
end

if nargout > 0, varargout{1} = h; end

return