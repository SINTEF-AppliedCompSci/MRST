% Script for simulating core flooding experiment with CO2 foam
% Results presented at the GHGT-15 conference in 2021:
%
% Grimstad, Alv-Arne and Bergmo, Per Eirik Strand and Barrabino, Albert and
% Holt, Torleif
% Modelling of CO2-foam Core Flooding Experiments
% Proceedings of the 15th Greenhouse Gas Control Technologies
% Conference 15-18 March 2021, http://dx.doi.org/10.2139/ssrn.3816402

%% Check script options
mrstModule add co2-foam
mrstModule add deckformat ad-core ad-blackoil ad-props

%% Initialize from deck
dataPath = getDatasetPath('co2_foam');
deck = readEclipseDeck(fullfile(dataPath, 'coreflooding', '1DFOAM_Z250_PC.DATA'));
deck = convertDeckUnits(deck);

%% Set up grid, rock, and fluid
G     = initEclipseGrid(deck);
G     = computeGeometry(G);
rock  = initEclipseRock(deck);
fluid = initDeckADIFluid(deck);
% For flooding from top of vertical core specify 'gravity x on'
% Vertical core makes interpretation of pore pressure difficult,
% so simulate as if horizontal core with 'gravity reset on'
%gravity reset on
gravity on z

%% Set up CO2 foam model
mrstModule add ad-eor

% Mobility reduction function parameters
foam = getFoam(15,'noShear',true); % Brij L23 foam properties

% Partitioning
Cpart = 0.02; % Data for Brij L23 at 200 bar, 40 °C

% Adsorption
adsMax = 1e-4; % 0.5 mg surfactant per g reservoir rock
adsSat = 5e-4; % Adsorption saturates at 0.05 wt% concentration in water

fluidF = addSimpleFoamProperties(fluid, ...
    'foam', foam, ...
    'cCrit', 0.001, ...         % Critical concentration ~0.1% by weight
    'cDecl', 50,    ...         % Foam effect zero at 0.002% by weight 
    'adsMax',adsMax, ...        % Maximum adsorption 
    'adsSat',adsSat, ...        % Threshold concentration for adsorption
    'Cpart', Cpart, ...         % Partitioning constant
    'surfingas', false, ...     % Gas-soluble surfactant
    'mobMultIsLinear',true, ... % Is the mobility multiplier a linear function of gas velocity?
    'noShear', true);           % No velocity dependence
% Turn off adsorption
%fluidF = rmfield(fluidF, 'adsMax');

model = GasWaterSurfactantFoamModel(G, rock, fluidF);
modelOld.OutputStateFunctions = [model.OutputStateFunctions, 'PhaseFlux'];

% No partitioning if commented
model.fluid.Cpart = Cpart;
% Initial state
state0 = initStateDeck(model, deck);

%% Add foam-specific pproperties to initial state
% Initial surfactant:
%   0.5 wt% in water throughout core
% Equilibrated with surfactant adsorbed to rock in case of adsorption,
% and with surfactant in gas in case of partitioning.
% Total concentration calculated based on the above.
% I.e., calculate the amount adsorbed when concentration in water is 0.5
% wt%, and add this to the total amount per m3 of space.

% Surfactant in water:
initialcW = 0.5/100; % 0.5 wt%
% Mass of surfactant
mSf = model.operators.pv.*model.fluid.rhoWS.*model.fluid.bW(state0.pressure).* ...
    state0.s(:,1)*initialcW/(1-initialcW);
% Surfactant in gas:
initialcG = Cpart*initialcW;
% No initial gas saturation, so no additional surfactant in gas phase
% Initial surfactant concentration, in kg per m3 pore volume
% Adsorbed surfactant
%state0.cA = zeros(G.cells.num,1);
state0.cA = ones(G.cells.num,1)*model.fluid.surfads(initialcW);
% Mass of surfactant
mSf = mSf + state0.cA.*G.cells.volumes.*(1-model.rock.poro).*model.fluid.rhoRSft;
% Initial volumetric concentration
c0 = (mSf)./model.operators.pv;
[state0.cs, state0.csmax] = deal(ones(G.cells.num,1).*c0);

% Stored surfactant concentrations in water, gas and adsorbed
state0.cW = ones(G.cells.num,1)*initialcW;
state0.cG = zeros(G.cells.num,1);

%% Set up simulation schedule
schedule = convertDeckScheduleToMRST(model, deck);
% Initialize surfactant concentration in wells
[schedule.control.W.cs] = deal(0);
% No surfactant injection
schedule.control(1).W(1).cs = 0;

% Simulate 30 hours with maximum time step 30 minutes, and 5 time step
% ramp-up of time step length.
tvec = rampupTimesteps(10*hour, 10*minute, 5);
schedule.step.control = ones(size(tvec));
schedule.step.val = tvec;

reservoirTime = cumsum(schedule.step.val);

%% Set up problem and solvers

mrstModule add linearsolvers

lsol = selectLinearSolverAD(model);
nls = NonLinearSolver();
nls.useRelaxation = true;
nls.LinearSolver = lsol;

model.dpMaxRel = 0.1;
model.dsMaxAbs = 0.1;

caseName = 'coreflooding-example';
problem = packSimulationProblem(state0, model, schedule, ...
    'Core flooding', 'nonlinearsolver', nls, 'Name', caseName);

%% Solve the schedule
simulatePackedProblem(problem, 'restartStep', 1);
%% Read results form disk
[wellSols, states, reports] = getPackedSimulatorOutput(problem);

%% Plot results

% Load results
[ws, states, reports] = getPackedSimulatorOutput(problem);
% Generate plots
mrstModule add mrst-gui
figure;
plotToolbar(G, states,'plot1d',true,'field','s:2');
title(caseName,'Interpreter','none');
plotWellSols(ws, reservoirTime)
title(caseName,'Interpreter','none');

%% Copyright Notice 
%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF
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