function plotGmsh(G)
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

    rng(1);

    if G.griddim == 2

        % Cells
        figure, hold on
        ut = unique(G.cells.tag);
        colors = [lines(7); rand(numel(ut)-7, 3)];

        for k = 1:numel(ut)
            idx = G.cells.tag == ut(k);
            plotGrid(G, idx, 'facecolor', colors(k, :));
        end

        % Faces
        figure, hold on
        plotGrid(G)
        ut = unique(G.faces.tag);
        colors = [lines(7); rand(numel(ut)-7, 3)];

        for k = 1:numel(ut)
            idx = G.faces.tag == ut(k);
            plotFaces(G, idx, 'linewidth', min(2, k), 'edgecolor', colors(k, :));
        end

    elseif G.griddim == 3

        % Cells
        figure; hold on
        ut = unique(G.cells.tag);
        colors = [lines(7); rand(numel(ut)-7, 3)];
        plotGrid(G, 'facecolor', 'none', 'facealpha', 0.0)
        view(3)

        for k = 1:numel(ut)
            if ut(k) > 0
                idx = G.cells.tag == ut(k);
                plotGrid(G, idx, 'facecolor', colors(k, :));
            end
        end

        % Faces
        fig = figure; hold on;
        ut = unique(G.faces.tag);
        colors = [lines(7); rand(numel(ut)-7, 3)];
        plotGrid(G, 'facecolor', 'none', 'facealpha', 0.0)

        for k = 1:numel(ut)
            if ut(k) > 0
                figure
                idx = G.faces.tag == ut(k);
                plotFaces(G, idx, 'linewidth', 2, 'facecolor', colors(k, :));
                view(3)
                title(sprintf('tag %d: %d entities', ut(k), sum(idx)));

                figure(fig)
                plotFaces(G, idx, 'linewidth', 2, 'facecolor', colors(k, :));
            end
        end

    end

end
