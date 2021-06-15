function plotOutlinedGrid(G, W, bc, tag, varargin)
%Utility to plot a grid outline that highlights certain faces
%
% SYNOPSIS:
%   plotOutlinedGrid(G, W, bc, tag)
%
% REQUIRED PARAMETERS:
%   G   - Grid structure
%   W   - Wells in model
%   bc  - BC in model
%   tag - Faces that should be plotted as highlighted
%
% RETURNS:
%   Nothing.
%
% NOTE:
%   Mostly useful for 2D grids.

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

    opt = struct('plotFaces', {{'xmin', 'xmax', 'zmin', 'zmax'}}, ...
                 'markersize', 12, ...
                 'woffset', [0, 10, 0]);
    opt = merge_options(opt, varargin{:});
    
    faces = false(G.faces.num, 1);
    for i = 1:numel(opt.plotFaces)
        faces(boundaryFaceIndices(G, opt.plotFaces{i})) = true;
    end
    
    if isempty(bc)
        bc.face = [];
    end
    faces(bc.face) = false;

    plotFaces(G, faces, 'linewidth', 2);
    
    if isnumeric(tag)
        faces = find(tag == 0);
    else
        faces = find(tag);
    end
    
    if any(faces)
        plotFaces(G, faces, 'facec', 'k', 'linewidth', 2)
    end
    plotFaces(G, bc.face, 'edgecolor', [1, 0, 0], 'linewidth', 2)

    hold on
    for i = 1:numel(W)
        x = G.cells.centroids(W(i).cells, 1);
        y = G.cells.centroids(W(i).cells, 2);
        z = G.cells.centroids(W(i).cells, 3);
        plot3(x - opt.woffset(1), y - opt.woffset(2), z - opt.woffset(3), 'kO', 'markerfacecolor', 'r', 'markersize', opt.markersize);
    end
end