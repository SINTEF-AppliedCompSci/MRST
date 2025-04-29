% Script for simulating field-scale application of CO2 foam
%
% An early version of this was presented at the GHGT-14 conference in 2018
% (Grimstad et al. (2018), DOI: 10.2139/ssrn.3365966)

%% Check script options and add modules
mrstModule add co2-foam ad-eor deckformat ad-core ad-blackoil ad-props

%% Initialize from deck
dataPath = getDatasetPath('co2_foam');
fn = fullfile(dataPath, 'field', 'RUN', 'QBOX2_GAS.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

%% Set up grid, rock, and fluid
G     = initEclipseGrid(deck);
G     = computeGeometry(G);
rock  = initEclipseRock(deck);
fluid = initDeckADIFluid(deck);

gravity reset on

%% Set up CO2 foam model
surfactantSystem = 15; % Brij L23
case_name = 'BOX_K1_C1';

% Mobility reduction function parameters
foam = getFoam(surfactantSystem,'noShear',true);  % Foam properties

% Surfactant properties
surf = getSurfactant(surfactantSystem); % Dissolution properties

% Override partition coefficient from getSurfactant().
surf.Cpart = 0.1;

fluidF = addSimpleFoamProperties(fluid, ...
    'foam', foam, ...
    'cCrit', 2e-3, ...          % Critical concentration by weight
    'cDecl', 2.5e3, ...         % Foam effect decline interval
    'adsMax',surf.adsMax, ...   % Maximum adsorption 
    'adsSat',surf.adsSat, ...   % Threshold concentration for adsorption
    'surfingas', false, ...     % Gas-soluble surfactant
    'mobMultIsLinear',true, ... % Is the mobility multiplier a linear function of gas velocity?
    'usePermDep',false, ...     % Is foam strength dependent on absolute permeability?
    'noShear', true);           % No velocity dependence

model = GasWaterSurfactantFoamModel(G, rock, fluidF);
model.OutputStateFunctions = [model.OutputStateFunctions,'PhaseFlux'];
% No partitioning if commented
model.fluid.Cpart = surf.Cpart;
% Set up inital state
state0 = initStateDeck(model, deck);

%% Add foam-specific pproperties to initial state
% This can be used e.g. to initialize a region around the injection well
% with some surfactant dissolved in the water, as if a slug of surfactant
% solution has already been injected.
%
% Initial surfactant:
% X wt% in water
% This should be equilibrated with surfactant adsorbed to rock in case of
% adsorption, and with surfactant in gas in case of partitioning.
%
% Total concentration calculated based on concentration in water and
% adsorbed. I.e., calculate the amount adsorbed when concentration in water
% is X wt%, and add this to the total amount per m3 of space.
%
% Insert calculations below
% Surfactant in water:
%initialcW = . . .
%mSf = model.operators.pv.*model.fluid.rhoWS.*model.fluid.bW(state0.pressure).* ...
%    state0.s(:,1)*initialcW/(1-initialcW);
%
% Surfactant in gas:
% initialcG = Cpart*initialcW;
%
% No initial gas saturation, so no additional surfactant in gas phase
% Initial surfactant concentration, in kg per m3 pore volume

% Adsorbed surfactant
state0.cA = zeros(G.cells.num,1);
%state0.cA = ones(G.cells.num,1)*model.fluid.surfads(initialcW);
%mSf = mSf + state0.cA.*G.cells.volumes.*(1-model.rock.poro).*model.fluid.rhoRSft;

% Initial volumetric concentration
[state0.cs, state0.csmax] = deal(zeros(G.cells.num,1));

% Stored surfactant concentrations in water, gas and adsorbed
state0.cW = zeros(G.cells.num,1);
state0.cG = zeros(G.cells.num,1);


%% Set up simulation schedule
schedule = convertDeckScheduleToMRST(model, deck);
% Initialize surfactant concentration in wells
% Check units for this parameter!
[schedule.control.W.cs] = deal(0); % Default to no surfactant injection
schedule.control(1).W(1).cs = 0.01; % 1.0 wt% surfactant.
schedule.control(2) = schedule.control(1);
schedule.control(2).W(1).cs = 0; % Stop surfactant injection.

% Time step vector
% Initial period with surfactant injection (control 1) followed by period
% without surfactant (control 2).
tvec = rampupTimesteps(4*year, 90*day);
n1 = length(tvec);
tvec = [tvec; rampupTimesteps(30*year, 180*day, 3)];
schedule.step.control = ones(size(tvec))*2;
schedule.step.control(1:n1) = 1;
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

problem = packSimulationProblem(state0, model, schedule, ...
    'fieldscale-co2-foam', 'nonlinearsolver', nls, 'Name', case_name);

%% Solve the schedule
simulatePackedProblem(problem, 'restartStep', 1);

%% Read results form disk
[wellSols, states, reports] = getPackedSimulatorOutput(problem);

%% Plot results
% Generate plots
mrstModule add mrst-gui
figure;
plotToolbar(G, states,'plot1d',false,'field','s:2');
title(case_name,'Interpreter','none');
plotWellSols(wellSols, reservoirTime)
title(case_name,'Interpreter','none');

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