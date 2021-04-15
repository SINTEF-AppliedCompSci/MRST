%% Create a grid and add helpers
G = computeGeometry(cartGrid([10, 20, 3]));
G = gridAddHelpers(G);

%% Demonstrate plotting
% We can now use the helper structures to plot data.
clf;
subplot(2,1,1);
G.plot.grid();
% Which is equivialent to
% plotGrid(G);

subplot(2,1,2);
data = randn(G.cells.num,1);
G.plot.cellData(data, 'FaceAlpha', .5);
% Which is equivialent to
% plotCellData(G, data, 'FaceAlpha', .5);

%% Demonstrate helpers
% Pick out a subset of cells, and plot them along with their corresponding
% points.
cells = [1, 2, 10];
[nodes, nodePos] = G.helpers.getCellNodes(cells);

clf;
G.plot.grid(cells);
G.plot.points(G.nodes.coords(nodes, :))

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
