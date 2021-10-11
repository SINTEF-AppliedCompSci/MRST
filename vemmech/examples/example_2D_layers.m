%% Example using three layers with different parameters in each layer
% We set up a 

mrstModule add vemmech

%% Set up the rock parameters and the grid
%
% Use the utility function squareLayersTest to set up the problem. You may
% look at the documentation of this function and try the other settings
% presented there by modifying this example.

E  = [1; 10; 1];      % Define the Youngs moduli in the layers
nu = [0.3; 0.3; 0.3]; % Define the Poisson's ratio in the layers

zlayers = [0.3; 0.3; 0.3]; % Define the width of the layers.
zlayers = 1/sum(zlayers)*zlayers;
L    = [1 1];   % Define the domain size

% Define the grid type, for example 'mixed4'
grid_type = 'mixed4';
dims= [4 4];

% Define test case type. It determines the boundary conditions and the forces
test_name = 'hanging_rod_2D';
% Parameter used to twist the grid so that the Cartesian symmetries are
% broken and we can test the method on more irregular grid.
disturb = 0.042;

% The utility function squareLayersTest sets up the problem.
[G, bc, test_cases] = squareLayersTest('E', E, 'nu', nu, 'zlayers', ...
    zlayers, 'cartDims', dims, 'L', L, 'grid_type', grid_type, 'disturb', ...
    disturb, 'make_testcases', false, 'test_name', test_name);
G = computeGeometry(G);
figure,
plotGrid(G); 

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
plotCellDataDeformed(G, mdiv, uVEM, plotopts1{:}); colorbar(); axis tight

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
