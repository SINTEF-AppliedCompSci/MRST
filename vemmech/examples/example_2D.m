%% Linear Elasticity on Complex Grid
% We consider a case with a complex block-structured grid, in which each
% individual block is gridded using different grid types: curvilinear
% quadrilateral blocks, unstructured Voronoi grids, and triangular grids.
% The purpose of the example is to highlight the flexibility of the virtual
% element method.

mrstModule add vemmech

%% Set up the rock parameters and the grid
%

E  = 1;    % Young's modulus
nu = 0.3;  % Poisson's ratio

L = [1 1]; % Grid domain
dims = [16 16]; % The grid size

% Parameter used to twist the grid so that the Cartesian symmetries are
% broken and we can test the method on more irregular grid.
disturb = 0.03;

test_name = 'original_2D_test'; % list of test names to be found in file
                                % squareTest.
grid_type = 'mixed4'; % Run the script exploreSquareGrid to see the grids
                      % that have been set up to test this VEM code with
                      % respect to different grid features.

% We use the utility function squareTest to set up the problem, that is
% defining the load and boundary conditions.
[G, bc, test_cases] = squareTest('E', E, 'nu', nu, 'cartDims', dims, ...
                                 'L', L, 'grid_type', grid_type, 'disturb', disturb, ...
                                 'make_testcases', false, 'test_name', test_name);
G = computeGeometry(G);
clf, plotGrid(G);

% We recover the problem parameters using the structure test_cases
testcase = test_cases{1};
el_bc    = testcase.el_bc; % The boundary conditions
load     = testcase.load;  % The load
C        = testcase.C;     % The material parameters in Voigts notations

%% We assemble and solve the equations
%
[uVEM, extra] = VEM_linElast(G, C, el_bc, load);

%% We plot the results
%
plotopts = {'EdgeAlpha', 0.0, 'EdgeColor', 'none'};
plotopts1 = {'EdgeAlpha', 0.1};
figure('Position',get(0, 'DefaultFigurePosition').*[1 1 1 2]);
subplot(3, 1, 1);
plotNodeData(G, uVEM(:, 1), plotopts{:}); colorbar();
title('x displacement');
subplot(3, 1, 2);
plotNodeData(G, uVEM(:, 2), plotopts{:}); colorbar();
title('y displacement');
subplot(3, 1, 3);
% We compute the divergence
vdiv = VEM2D_div(G);
mdiv = vdiv*reshape(uVEM', [], 1)./G.cells.volumes;
title('Divergence');
plotCellDataDeformed(G, mdiv, uVEM, plotopts1{:}); colorbar();

%%
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
