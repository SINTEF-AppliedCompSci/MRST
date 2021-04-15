%% Set up a compositional test model
mrstModule add ad-props compositional
G = cartGrid([200, 200, 1]);
G = computeGeometry(G);
cmodel = GenericOverallCompositionModel(G, makeRock(G, 1, 1), initSimpleADIFluid(), getBenchmarkMixture('spe5'));
cmodel = cmodel.validateModel(); % Set up groups
%% Set up a compositional state
useAD = false;
statec = initCompositionalState(G, 10*barsa,  273.15 + 30, [0.3, 0.4, 0.3], rand(1, 6), cmodel.EOSModel);
stateAD = cmodel.validateState(statec);
stateAD = cmodel.getStateAD(stateAD, useAD); % No AD initialization
% Show the initialized container
disp(stateAD.PVTProps)
%% Evaluate density twice - second call is cached
fprintf('First evaluation: '); tic(); rho = cmodel.getProp(stateAD, 'Density'); toc();
fprintf('Second evaluation: '); tic(); rho = cmodel.getProp(stateAD, 'Density'); toc();
%% Show the PVT property cache after evaluation
disp(stateAD.PVTProps)
%% Plot density relations for a compositional model
figure;
[~, gcomp] = plotStateFunctionGroupings(cmodel,'Stop', 'Density', 'includeState', false);
clc
printStateFunctionGroupingTikz(gcomp);

%% Set up a black-oil model
pth = getDatasetPath('spe1');
fn  = fullfile(pth, 'BENCH_SPE1.DATA');
% deck = readEclipseDeck(fn);
[state0, bomodel, schedule, nonlinear] = initEclipseProblemAD(fn);
bomodel = bomodel.validateModel();

%% Plot density relations for a black-oil model
figure;
[h, gbo] = plotStateFunctionGroupings(bomodel, 'Stop', 'Density', 'includeState', false);
clc
printStateFunctionGroupingTikz(gbo);

%% Plot all functions
figure;
[~, gbos] = plotStateFunctionGroupings(bomodel, 'includeState', false);
clc
printStateFunctionGroupingTikz(gbos);
%% Plot all functions (as symbols)
figure;
[h2, gbos2, g] = plotStateFunctionGroupings(bomodel, 'includeState', false, 'label', 'label');
clc
printStateFunctionGroupingTikz(gbos2, h2);
%% Plot all functions for flux discretization (as symbols)
groups = bomodel.getStateFunctionGroupings();
% groups = groups(~cellfun(@(x) isa(x, 'FacilityFlowDiscretization'), groups));
figure;
[h2, gbos2, g] = plotStateFunctionGroupings(groups, 'includeState', false, 'label', 'label');
%%
clc
printStateFunctionGroupingTikz(gbos2, h2, 'outer', 'tree layout, grow=right', 'inner', 'tree layout, components go down, grow = down');

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
