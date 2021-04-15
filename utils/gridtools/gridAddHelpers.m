function G2 = gridAddHelpers(G)
% Add helpers to existing grid structure for cleaner code structure.
%
% SYNOPSIS:
%   G = gridAddHelpers(G);
%
% PARAMETERS:
%   G    - Grid structure
%
% RETURNS:
%   G    - Grid structure with added fields:
%
%             - plot contains handles to plotting functions
%             - helpers contains various helper functions commonly used
%               when plotting/debugging grid structures.
%

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


    assert(isfield(G.cells, 'centroids'), ['Grids passed to gridAddHelpers', ...
    'must have geometry information. Please call computeGeometry first...'])
    G2 = G;

    % Plotting
    G2.plot.grid      = @(varargin) plotGrid(G, varargin{:});
    G2.plot.cellData  = @(data, varargin) plotCellData(G, data, varargin{:});
    G2.plot.faces     = @(varargin) plotFaces(G, varargin{:});
    G2.plot.faceData  = @(data, varargin) plotFaceData(G, data, varargin{:});
    G2.plot.points    = @(pts) plotPts(pts, numel(G.cartDims));


    % Helpers for grid operations
    G2.helpers.cellNo            = @(varargin) gridCellNo(G, varargin{:});
    G2.helpers.getLogicalIndices = @(varargin) gridLogicalIndices(G, varargin{:});
    G2.helpers.getFaceNodes      = @(faces) gridFaceNodes(G, faces);
    G2.helpers.getCellNodes      = @(cells) gridCellNodes(G, cells);
    G2.helpers.getCellFaces      = @(cells) gridCellFaces(G, cells);

    G2.type = [G2.type, {'gridAddHelpers'}];
end

function plotPts(pts, d)
    hold on;
    if d==3
        plot3(pts(:, 1), pts(:, 2), pts(:, 3), '.k');
    else
        plot(pts(:, 1), pts(:, 2), '.k');
    end
    hold off;
end
