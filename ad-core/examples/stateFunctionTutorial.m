%% Introduction to StateFunctions for the AD-OO framework
% The StateFunction class, together with the StateFunctionGrouping class,
% is the main mechanism for evaluating properties during a AD-OO
% simulation. These functions act directly on the "state" object and fully
% support automatic differentiation, lazy evaluation and different regions.
mrstModule add ad-core ad-blackoil deckformat ad-props
%% Create a black-oil test model
% We set up the SPE1 model, a black-oil model with disgas (but no vapoil).
pth = getDatasetPath('spe1');
fn  = fullfile(pth, 'BENCH_SPE1.DATA');
[state0, model, schedule, nonlinear] = initEclipseProblemAD(fn);
%% Validate model to set up defaults
% State functions are set up at the start of a simulation, with reasonable
% defaults given. Here, we manually call validateModel to invoke these
% defaults.
model = model.validateModel();
% We can now get the groupings belonging to this model
groups = model.getStateFunctionGroupings();
%% Test the evaluation
% Each grouping is a collection of one or more StateFunctions, which can
% depend on each other, or the state of the physical system itself.
% Dependencies between these are automatically handled during simulations.
% For instance, let us get the mobility of the initial state. We then plot
% the oil mobility on the grid.

% We can get this with getProp, just as we would with the pressure or
% saturations which are stored directly in the state. Normally, we do not
% have to think too hard about where exactly a property is defined, since
% the framework can figure it out for us.
relperm = model.getProp(state0, 'RelativePermeability');
% Alternatively, we could invoke the FlowPropertyFunctions grouping's get
% function to the same effect.
kr = model.FlowPropertyFunctions.get(model, state0, 'RelativePermeability');

figure;
plotCellData(model.G, relperm{2});
title('Initial oil relperm')
view(30, 30);
colorbar
%% We can also initialize AD-variables in state to get derivatives
% Initialize pressure as a AD-variable
stateAD = state0;
stateAD.pressure = initVariablesADI(state0.pressure);
% Outputs are now ADI, with differentiation with respect to pressure
mobAD = model.getProp(stateAD, 'Mobility');
disp(mobAD)
%% Show a state function grouping
% We can inspect the FlowPropertyFunctions directly. The
% FlowPropertyFunctions consist of the basic properties required for flow
% in MRST. Additional properties can be added dynamically to the class
% instance itself.
disp(model.FlowPropertyFunctions)
%% Plot dependencies
% Generally, only a few values are needed to assemble the final linearized
% equation. Many other functions may be required to get intermediate
% results, however. For instance, capillary pressure does not occur
% directly in the flow equations, but standard functions for density, phase
% pressures and viscosity can depend strongly on the capillary pressure
% between phases.
%
% We can plot the dependency graph of all the state function groupings in
% order to understand the relationships between the different functions.
for i = 1:numel(groups)
    figure;
    plotStateFunctionGroupings(groups{i},  'label', 'label')
    title(class(groups{i}));
end
%% Plot a specific dependency
% Plot all dependencies for the mobility we just evaluated
figure(1); clf
plotStateFunctionGroupings(model.FlowPropertyFunctions, 'Stop', 'Mobility')
title('All dependencies required to evaluate mobility')
%% Plot all dependencies of a property
% Let's see all properties which directly or indirectly depend on on the
% pressure in the state. We get the PVT properties together with the flow
% properties for this purpose.
groups = {model.FlowPropertyFunctions, model.PVTPropertyFunctions};
figure(1); clf
plotStateFunctionGroupings(groups, 'Start', 'pressure')
title('Flow properties that depend on pressure') 
%% Plot everything which either depends upon a property, or uses that property
% Let's see all properties which directly or indirectly depend on on the
% pressure in the state
figure(1); clf
plotStateFunctionGroupings(groups, 'Center', 'Mobility')
title('Upstream and downstream mobility dependencies');
%% We can combine graphs
% Here, we plot all the flow properties together with the discretization of
% the flux, which depends on many of these properties.
df = get(0, 'DefaultFigurePosition');
figure('Position', df.*[1, 1, 2, 2]);
tmp = {model.FlowPropertyFunctions, model.FlowDiscretization, model.PVTPropertyFunctions};
plotStateFunctionGroupings(tmp);
title('Flow properties + flow discretization')
%% Create a compositional model and visualize the property graphs
% We create a compositional model, which has a few major differences from
% the black-oil model: Note the addition of equation-of-state properties
% such as fugacity or phase compressibility factors. In addition, the
% relationship between the ShrinkageFactors and Density
mrstModule add compositional
eos = getBenchmarkMixture('simple');
cmodel = GenericNaturalVariablesModel(model.G, model.rock, model.fluid, eos, 'water', true);
cmodel = cmodel.validateModel();

clf;
[h, g] = plotStateFunctionGroupings(cmodel.FlowPropertyFunctions);
%% Print as TiKZ figure
% If you have a LaTeX+PGF compiler, you can also make publication-quality plots
% using the experimental printStateFunctionGroupingTikz routine. The code
% is printed to the command window.
clc
printStateFunctionGroupingTikz(g)
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
