function [h, saturation] = plotSaturationDG(disc, state, varargin)
    % Visualize dG saturation with higher-order dG basis functions

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

    opt = struct('n'      , 100  , ...
                 'phaseNo', 1    , ...
                 'plot1d' , false, ...
                 'coords' , []   , ...
                 'contour', false, ...
                 'makePlot', true, ...
                 'pbaspect', [1,1,1]);
    [opt, extra] = merge_options(opt, varargin{:});
    coords = opt.coords;
    if isempty(coords)
        coords = getPlotCoordinates2(disc.G, varargin{:});
    end
    x     = coords.points;
    cells = coords.cells;
    xs    = cell(1, disc.G.griddim);
    saturation = @(state) disc.evaluateDGVariable(x, cells, state, state.sdof(:,opt.phaseNo));
    h = [];
    if opt.makePlot && ~isempty(state)
        s = saturation(state);
        if opt.plot1d
            h = plot(x(:,1), s, extra{:});
        elseif disc.G.griddim == 2
            if opt.contour
                h = contourf(xs{1}, xs{2}, s, extra{:});
            else
                x = [coords.points, s];
                h = patch('faces', coords.faces, 'vertices', x, 'FaceVertexCData', s, 'FaceColor','interp', extra{:});
                zlim([0,1]);
                pbaspect(opt.pbaspect);
            end
        else
            error('Only supported for 1D and 2D')
        end
    end

end
