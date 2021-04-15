%% Streamlines for a 2D Quater Five-Spot Well Pattern
% The quarter five-spot problem is one of the cases that are most widely
% used to test numerical methods. The setup consist of an injector and a
% producer located in diagonally oposite corners, with no-flow conditions
% along the model perimeter. This gives a symmetric flow pattern with high
% flow along the diagonal and stagnant points at the corners with no wells.

mrstModule add mimetic incomp streamlines

%% Construct model and compute flow field
G = cartGrid([25,25]);
G = computeGeometry(G);

rock = makeRock(G, 10*milli*darcy, 0.3);

fluid = initSimpleFluid('mu',  [1, 1]*centi*poise', ...
                        'rho', [1000, 1000]*kilogram/meter^3, ...
                        'n',   [2, 2]);

T = computeTrans(G, rock);
src = addSource([], [1, G.cells.num], [1.1, -1.1]);

x = initResSol(G, 0, 0);
x = incompTPFA(x, G, T, fluid, 'src', src);

%% Trace streamlines
% Pick start points on the diagonal orthogonal to the well-pair direction
% and trace streamlines forward and backward from these points. This will
% ensure maximal accuracy in the near-well regions where the streamlines
% are converging.
clf
h = plotGrid(G, 'facea', 0.3, 'edgea',0.1);
hold on;
cells = (G.cartDims(1):G.cartDims(1)-1:prod(G.cartDims)-1)';
streamline(pollock(G, x, cells));
x.flux = -x.flux;
streamline(pollock(G, x, cells, 'substeps', 1));
plot(G.cells.centroids(cells,1), G.cells.centroids(cells,2), ...
    'or','MarkerSize',8,'MarkerFaceColor',[.6 .6 .6]);
axis equal tight
hold off

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
