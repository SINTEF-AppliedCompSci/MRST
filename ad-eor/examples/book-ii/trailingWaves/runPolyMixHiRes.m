%% Resolution of trailing waves with high-resolution schemes
% The example runs a simple displacement for which the c-wave should be a
% shock (Todd-Longstaff mixing parameter w=0.85). The purpose of the
% example is to demonstrate that using a high-resolution spatial
% discretization significantly improves the resolution of the trailing
% chemical wave.
mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat

%% Run fine-grid simulation
gravity reset on;
bookdir = getDatasetPath('eor_book_ii');
fn = fullfile(bookdir,'trailingWaves','MixPar400.DATA');
[state0, model, schedule, nlsolver] = initEclipseProblemAD(fn);
model.fluid.mixPar = 0.85;
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);

clf
subplot(1,2,1), hold all
plot(model.G.cells.centroids(:,1),states{end}.s(:,1),'LineWidth',2);
subplot(1,2,2), hold all
plot(model.G.cells.centroids(:,1),states{end}.cp(:,1),'LineWidth',2);

%% Run simulation with 50 cells
fn = fullfile(bookdir,'trailingWaves','MixPar50.DATA');
[state0, model, schedule, nlsolver] = initEclipseProblemAD(fn);
model.fluid.mixPar = 0.85;
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);

pargs = {'-','Marker','.','MarkerSize', 14, 'LineWidth',1};
subplot(1,2,1)
plot(model.G.cells.centroids(:,1),states{end}.s(:,1),pargs{:});
subplot(1,2,2)
plot(model.G.cells.centroids(:,1),states{end}.cp(:,1),pargs{:});

%% Run the same with WENO and with WENO + AIM
model = setWENODiscretization(model);
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);
subplot(1,2,1)
plot(model.G.cells.centroids(:,1),states{end}.s(:,1),pargs{:});
subplot(1,2,2)
plot(model.G.cells.centroids(:,1),states{end}.cp(:,1),pargs{:});

model = setTimeDiscretization(model, 'aim', 'saturationCFL', 0.75);
[~, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nlsolver);
subplot(1,2,1)
plot(model.G.cells.centroids(:,1),states{end}.s(:,1),pargs{:});
subplot(1,2,2)
plot(model.G.cells.centroids(:,1),states{end}.cp(:,1));
legend('SPU: \Delta x = 1 m','SPU: \Delta x = 8 m', ...
    'WENO: \Delta x = 8 m', 'WENO+AIM: \Delta x = 8 m');

%% Copyright Notice
%
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
