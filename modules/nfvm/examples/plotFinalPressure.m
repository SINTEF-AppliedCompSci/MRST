function plotFinalPressure(G, states, name)
%Undocumented Helper Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    figure
    plotToolbar(G, states)
    axis tight equal
    colorbar
    title(name)

    if isfield('G', 'cartDims')
        figure, hold on
        plotCellData(G, states{end}.pressure, 'edgealpha', 0);
        contour(reshape(G.cells.centroids(:, 1), G.cartDims), ...
                reshape(G.cells.centroids(:, 2), G.cartDims), ...
                reshape(states{end}.pressure, G.cartDims), ...
                'linewidth', 1, 'color', 'k');
        axis tight equal
        colorbar
        title([name, ' at endtime'])
    else
        unstructuredContour(G, states{end}.pressure, 100);
        title(name)
    end

end
