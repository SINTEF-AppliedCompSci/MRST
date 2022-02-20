function handles = plotCellVectorData(G, data, varargin)
%Plot a vector field defined per cell
%
% SYNOPSIS:
%   h = plotCellVectorData(G, vectordata);
%
%   % Increase line width, 10 different colors.
%   h = plotCellVectorData(G, vectordata, 'linewidth', 3, 'n', 10);
%
% DESCRIPTION:
%   plotCellVectorData plots cell-wise vector fields of the same dimension
%   as the grid. The values are normalized to the size of the grid, giving
%   relatively useful output even on fully unstructured datasets. Current
%   colormap is used to colorize the vectors, but as they do not use CData
%   the vector plot can be combined with e.g. a pressure field plotted by
%   plotCellData.
%
% REQUIRED PARAMETERS:
%       G    - Valid grid structure
%
%       data - Dataset of size(G.cells.num, G.griddim) interpreted as a
%              vector field.
%
%   cells   - Vector of cell indices defining sub grid.  The graphical
%             output of function 'plotCellVectorData' will be restricted to
%             the subset of cells from 'G' represented by 'cells'.
%
%             If unspecified, function 'plotCellVectorData' will behave as
%             if the user defined
%
%                 cells = 1 : G.cells.num
%
%             meaning graphical output will be produced for all cells in
%             the grid model 'G'.  If 'cells' is empty (i.e., if
%             ISEMPTY(cells)), then no graphical output will be produced.
%
% OPTIONAL PARAMETERS:
%
%   'linewidth' - Width of the underlying quiver plot vectors. Default 2.
%
%   'Scale'     - Adjustment parameter for scaling of the vector arrows.
%                 If the default is too small/large, this can be
%                 increased/decreased from the default value of 1.
%
%   'Colormap'  - Function handle @(N) which produces the colormap for N
%                 entries used to colorize the vector field subset.
%
%   'N'         - The vector data is binned using the number of unique
%                 elements, up to a maximum of 64. Increasing this will
%                 increase color resolution at the cost of execution speed.
%
% RETURNS:
%   h           - Handles to up to N different quiverplots.
%
% EXAMPLE:
%    cellFlux = faceFlux2cellVelocity(G, state.flux);
%    plotCellVectorData(G, cellFlux);
%
% SEE ALSO:
%   `plotCellData`, `plotToolbar`, `quiver3`, `quiver`

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


    % Set renderer to OpenGL to avoid painters panicking completely
    set(gcf, 'renderer', 'opengl')
    % We are in fact plotting a reservoir
    set(gca, 'ZDir', 'reverse')

    assert(isfield(G, 'cartDims'), ['This function requires cartDims'...
                                    ' field to be present in grid.']);
    subset = true(G.cells.num, 1);
    if mod(numel(varargin), 2) == 1
        subset = varargin{1};
        if ~islogical(subset)
            tmp = subset;
            subset = false(G.cells.num, 1);
            subset(tmp) = true;
        end
        varargin = varargin(2:end);
    end
    dfcmap = colormap();
    opt = struct('LineWidth', 2,...
                 'Scale', 1,...
                 'Colormap', @(n) interp1(1:size(dfcmap, 1), dfcmap, 1:size(dfcmap, 1)),...
                 'N', 64);

    opt = merge_options(opt, varargin{:});


    gc = G.cells.centroids(subset, :);

    % Try to work out the scale of the grid using cartDims
    scalefac = opt.Scale*2*(max(gc) - min(gc))...
                ./G.cartDims;

    ndata = data(subset, :);
    is3d = G.griddim == 3;

    N = min(numel(unique(ndata)), opt.N);

    % Scale by maximal vector norm
    vnorm = sqrt(sum(ndata.^2, 2));
    ndata = ndata./max(vnorm);
    pos = vnorm./max(vnorm);

    X = gc(:, 1);
    Y = gc(:, 2);
    if is3d
        Z = gc(:, 3);
    end

    H = 1/N;

    colors = opt.Colormap(N);

    handles = [];
    % Iterate over the bins, colorize and plot
    wasHold = ishold();
    hold on
    for i = 1:N
        subz = abs(pos - (i-1)*H) <= H;
        if ~any(subz); continue; end

        U = scalefac(1)*ndata(subz, 1);
        V = scalefac(2)*ndata(subz, 2);
        if is3d
            W = scalefac(3)*ndata(subz, 3);
        end

        if is3d
            h = quiver3(X(subz), Y(subz), Z(subz), ...
                       U, V, W, 0, 'linewidth', opt.LineWidth, 'color', colors(i, :));
        else
            h = quiver(X(subz), Y(subz), ...
                       U, V, 0, 'linewidth', opt.LineWidth, 'color', colors(i, :));
        end
        handles = [handles; h]; %#ok
    end
    if ~wasHold
        hold off
    end
end