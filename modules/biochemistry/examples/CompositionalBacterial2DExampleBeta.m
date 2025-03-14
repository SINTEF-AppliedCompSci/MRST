%% Illustrative Case with Compositional Model: Efficiency of Structural Trap for Hydrogen Storage
%
% This example adapts the illustrative case from the H2store black-oil simulator 
% into a compositional model with bio-chemistry (disabled here). The fluid and rock properties remain consistent 
% with the original model.
%--------------------------------------------------------------------------
clearvars; 
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr

%% Define the case name and read the Eclipse deck file
name = 'H2_STORAGE_COMPOSITIONAL_2D_TRAP_BACT_50_50';
%% Use H2STORAGE_RS_SALT.DATA for brine
deck = readEclipseDeck('\\wsl.localhost\Ubuntu\home\elyesa\Projects\MRST\modules\H2store\data\Illustrative_example\H2STORAGE_RS.DATA');

%% Set up the simulation parameters and model components from the black-oil model
[~, ~, state0Bo, modelBo, scheduleBo, ~] = H2_illustration_storage_example(deck);
bacteriamodel = true;
% Convert black-oil model to compositional model
model = convertBlackOilModelToCompositionalModel(modelBo);
%model.fluid.rhoGS = 1.4.*kilogram*meter^3;
state0=convertBlackOilStateToCompositional(modelBo,state0Bo);
state0.components = ensureMinimumFraction(state0.components,6.0e-4);

%% Set up the compositional model with Water and Hydrogen
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, {'H2O', 'H2', 'CO2', 'C1'});
EOS = SoreideWhitsonEquationOfStateModel([], compFluid, 'sw');
T0 = 40 * Kelvin; % Initial temperature in Kelvin
T0 = 273.15 * Kelvin + T0; % Convert to absolute temperature
if (compFluid.getNumberOfComponents>2)
    state0.components(:,3) = state0.components(:,2)/3;
    state0.components(:,4) = state0.components(:,3);
    state0.components(:,2) = state0.components(:,3);
    state0.components(:,1) = 1- sum(state0.components(:,2:end),2);
    state0.components =state0.components.*0+ [0.8480,  1.0e-5,  1.0e-5,   0.1760-0.0240];

    EOS = SoreideWhitsonEquationOfStateModel([], compFluid, 'sw');
    T0 = 44.35 * Kelvin; % Initial temperature in Kelvin
    T0 = 273.15 * Kelvin + T0; % Convert to absolute temperature
end

model.EOSModel = EOS;
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);

% Set up additional model properties for biochemistry and flow
arg = {model.G, model.rock, model.fluid, compFluid,...
       false, diagonal_backend, 'oil', true, 'gas', true, ... % Define phases for water-oil system
    'bacteriamodel', bacteriamodel, 'diffusioneffect', false, 'liquidPhase', 'O', ...
    'vaporPhase', 'G', 'eos', model.EOSModel}; % Set phases and EOS model
model = BiochemistryModel(arg{:});
model.gravity = modelBo.gravity;

%% Manually set up well compositions, temperature, and bacteria
schedule = scheduleBo; % Load schedule from the black-oil model
nbact0 = 1.0e2; % Initial number of bacteria
state0.T = state0.T.*0+T0;
% Initialize compositional state with bacteria model
if bacteriamodel
    state0 = initCompositionalStateBacteria(model, state0.pressure, T0, state0.s, state0.components, nbact0, EOS);
else
     state0 = initCompositionalState(model, state0.pressure, T0, state0.s, state0.components, EOS);
end
bc = schedule.control(1).bc; % Get boundary condition from the schedule
cells_bc = sum(model.G.faces.neighbors(bc.face, :), 2);
%% Update schedule controls and boundary conditions
if (compFluid.getNumberOfComponents>2)
    for i = 1:length(schedule.control)
        % We inject 50% H2 and 50% CO2 in the build-up phase
        schedule.control(i).W.compi = [0, 1]; % Well components
        if (strcmp(schedule.control(i).W.name, 'cushion')&&i<11)
            schedule.control(i).W.components = [0.0, 0.1, 0.9 0.0];
        else
            schedule.control(i).W.components = [0.0, 0.95, 0.05, 0.0];
        end
        schedule.control(i).W.T = T0; % Set temperature (not used, but necessary)

        % Update boundary conditions for each control
        schedule.control(i).bc.components = repmat(state0.components(1,:), numel(cells_bc), 1);
        schedule.control(i).bc.sat = repmat(state0.s(1,:), numel(cells_bc), 1);
        schedule.control(i).bc = [];
    end
else
    for i = 1:length(schedule.control)
        % Set component mix for each control, adjusting for specific conditions
        schedule.control(i).W.compi = [0, 1]; % Well components
        if (strcmp(schedule.control(i).W.name, 'cushion'))
            schedule.control(i).W.components = [0.0, 1.0];
        else
            schedule.control(i).W.components = [0.0, 1.0];
        end
        schedule.control(i).W.T = T0; % Set temperature (not used, but necessary)

        % Update boundary conditions for each control
        schedule.control(i).bc.components = repmat(state0.components(1,:), numel(cells_bc), 1); % Set boundary component concentrations
        schedule.control(i).bc.sat = repmat(state0.s(1,:), numel(cells_bc), 1); % Set boundary saturation
    end
end

model.outputFluxes = false;
%% Plot Grid with Wells, Permeability, and Porosity
figure;
subplot(1, 2, 1); 
plotCellData(model.G, model.rock.poro);
hold on;
plotGrid(model.G, schedule.control(1).W.cells, 'LineStyle', 'none', 'FaceColor', 'red');
title('Porosity');
axis off tight; 
subplot(1, 2, 2); 
plotCellData(model.G, log10(model.rock.perm(:, 1)));
hold on;
plotGrid(model.G, schedule.control(1).W.cells, 'LineStyle', 'none', 'FaceColor', 'red');
title('Permeability (log10)');
axis off tight; 
sgtitle('2D Dome-Shaped Aquifer'); 

%% Add output functions to the model for various properties
% model.OutputStateFunctions{end + 1} = 'CapillaryPressure';
% model.OutputStateFunctions{end + 1} = 'SurfaceDensity';
% model.OutputStateFunctions{end + 1} = 'ShrinkageFactors';
model.outputFluxes = false;

%% Initialize the nonlinear solver and select the linear solver
nls = NonLinearSolver(); 
lsolve = selectLinearSolverAD(model); 
nls.LinearSolver = lsolve;
%[ws, states, rep] = simulateScheduleAD(state0, model, schedule);
%% Pack the simulation problem with the defined components
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);

%% Execute the simulation of the packed problem
simulatePackedProblem(problem);
%% Compare with and without bectrial effects
modelNoBact = model;
modelNoBact.bacteriamodel = false;
state0NoBact = state0;
state0NoBact.nbact = 0;
nameNoBact = 'H2_STORAGE_COMPOSITIONAL_2D_TRAP_NOBACT_50_50';
problemNoBact = packSimulationProblem(state0, modelNoBact, schedule, nameNoBact, 'NonLinearSolver', nls);
%% Run the simulation
simulatePackedProblem(problemNoBact, 'restartStep',1);
%% Get reservoir and well states
[ws,states] = getPackedSimulatorOutput(problem);
[wsNoBact,statesNoBact] = getPackedSimulatorOutput(problemNoBact);
% Plot states
figure;
plotToolbar(model.G, states);
plotToolbar(model.G, statesNoBact);

%% Plot well output
figure;
plotWellSols({ws,wsNoBact});
%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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
