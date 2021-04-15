%% Streamlines for a 3D Quater Five-Spot Well Pattern
% The setup consist of a vertical injector and a vertical producer located
% in the southwest and northeast corners, respectively.  This gives a
% symmetric areal flow pattern with high flow along the diagonal and
% stagnant points at the southeast and northwest corners. The setup is
% almost the same as in the <matlab:edit('testPollock2D.m') 2D test>,
% except that we now use a consistent discretization from the |mimetic|
% module.

mrstModule add mimetic incomp streamlines

%% Set up model
G = cartGrid([25,25,2]);
G = computeGeometry(G);
rock = makeRock(G, 10*milli*darcy, 0.3);

fluid = initSimpleFluid('mu',  [1, 1]*centi*poise', ...
                        'rho', [1000, 1000]*kilogram/meter^3, ...
                        'n',   [2, 2]);

IP = computeMimeticIP(G, rock);
src = addSource([], [1, G.cells.num], [1, -1]);
x = initResSol(G, 0, 0);
x = incompMimetic(x, G, IP, fluid, 'src', src);

%% Trace streamline
clf
h = plotGrid(G, 'facea', 0.3, 'edgea',0.1);
hold on;
cells = (G.cartDims(1):G.cartDims(1)-1:prod(G.cartDims)/2-1)';
cells = [cells; cells+prod(G.cartDims)/2];
streamline(pollock(G, x, cells));
x.flux = -x.flux;
streamline(pollock(G, x, cells, 'substeps', 1));
plot3(G.cells.centroids(cells,1), G.cells.centroids(cells,2), ...
     G.cells.centroids(cells,3), ...
    'or','MarkerSize',8,'MarkerFaceColor',[.6 .6 .6]);
hold off
axis tight, view(30,40)

%% Copyright notice

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
