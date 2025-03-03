%==========================================================================
%% 1D Compositional Simulation with Bio-Clogging
%==========================================================================
% This script simulates a 1D two-phase compositional problem using MRST.
% Key features:
%   - Injection of 95% H2 and 5% CO2 into a reservoir with CO2, Methane, H2.
%   - Incorporates bio-clogging effects (porosity and permeability reduction).
%   - Uses Soreide-Whitson EOS asimnd LBC viscosity correlation.
%   - Relative permeability modifications for residual saturations.
%   - Bottom-hole pressure controls for wells.
%==========================================================================
clear all;
%% Initialize MRST and Load Modules
%--------------------------------------------------------------------------
% Add required MRST modules for compositional simulation.
%--------------------------------------------------------------------------
mrstModule add compositional deckformat ad-core ad-props;

%% Set Up Simulation Model
%--------------------------------------------------------------------------
% Initialize the model, schedule, and initial state from MRST validation example.
%--------------------------------------------------------------------------
% Use overall composition mode
[state0, model, schedule, ref] = setupSimpleCompositionalExample(false);
schedule.step.val = schedule.step.val/5;
schedule.step.val = schedule.step.val(1:400);
schedule.step.control = schedule.step.control(1:400);

% Set simulation name based on mode
name = 'comp-bio-clogging';

%% Define Fluid Properties and EOS
%--------------------------------------------------------------------------
% Use the Soreide-Whitson EOS for fluid modeling.
%--------------------------------------------------------------------------
eosname = 'sw';
model.EOSModel = SoreideWhitsonEquationOfStateModel(model.G, model.EOSModel.CompositionalMixture, eosname);

%% We modify the Relative Permeability
%--------------------------------------------------------------------------
% Adjust relative permeability curves to account for residual saturations.
%--------------------------------------------------------------------------
% Define synthetic saturation ranges
SW = linspace(0, 1, 100);  % Water saturation
SO = linspace(0, 1, 100);  % Oil saturation
SG = linspace(0, 1, 100);  % Gas saturation

% Evaluate original relative permeability functions
krW_original = arrayfun(model.fluid.krW, SW);  % Water rel perm
krOW_original = arrayfun(model.fluid.krOW, SO);  % Oil rel perm (water-oil)
krG_original = arrayfun(model.fluid.krG, SG);  % Gas rel perm
krOG_original = arrayfun(model.fluid.krOG, SO);  % Oil rel perm (gas-oil)

% Define residual saturations
swc = 0.2;  % Residual water saturation
sgr = 0.1;  % Residual gas saturation
sor = 0.2;  % Residual oil saturation

% Update water relative permeability
SW_shifted = max(SW - swc, 0);  % Shift water saturation
krW_new = interp1(SW, krW_original, SW_shifted, 'linear', 0);  % Interpolate
model.fluid.krW = @(sw) interp1(SW, krW_new, value(max(sw - swc, 0)));

% Update gas relative permeability
SG_shifted = max(SG - sgr, 0);  % Shift gas saturation
krG_new = interp1(SG, krG_original, SG_shifted, 'linear', 0);  % Interpolate
model.fluid.krG = @(sg) interp1(SG, krG_new, value(max(sg - sgr, 0)));

% Update oil relative permeability
SO_shifted = max(SO - sor, 0);  % Shift oil saturation
krOW_new = interp1(SO, krOW_original, SO_shifted, 'linear', 0);  % Interpolate
krOG_new = interp1(SO, krOG_original, SO_shifted, 'linear', 0);  % Interpolate
model.fluid.krOW = @(so) interp1(SO, krOW_new, value(max(so - sor, 0)));
model.fluid.krOG = @(so) interp1(SO, krOG_new, value(max(so - sor, 0)));

% Define relative permeability endpoints
model.fluid.krPts.w = [swc, 0.9];  % Water endpoints
model.fluid.krPts.g = [sgr, 0.8];  % Gas endpoints
model.fluid.krPts.ow = [0, 1 - swc];  % Oil endpoints (water-oil)
model.fluid.krPts.og = [0, 1 - sgr];  % Oil endpoints (gas-oil)

%% Add Bio-Clogging Effects
%--------------------------------------------------------------------------
% Define porosity and permeability reduction due to bacterial growth.
%--------------------------------------------------------------------------
% Initial rock properties
poro0 = model.rock.poro;  % Initial porosity
perm0 = model.rock.perm(:, 1);  % Initial permeability

% Bacterial concentration parameters
nc = 5e6;  % Critical bacterial concentration
cp = 1.5;  % Clogging coefficient

% Define porosity multiplier as a function of bacterial concentration
pvMult_nbact = @(nbact) 1 ./ (1 + cp .* log(1 + (nbact / nc).^2));
model.fluid.pvMultR = @(p, nbact) pvMult_nbact(nbact);
poro = @(p, nbact) poro0 .* pvMult_nbact(nbact);
model.rock.poro = poro;

% Define permeability as a function of porosity
tau = @(p, nbact) ((1 - poro0) ./ (1 - poro(p, nbact))).^2 .* (poro(p, nbact) ./ poro0).^3;
perm = @(p, nbact) perm0 .* tau(p, nbact);
model.rock.perm = perm;

% Define densities
model.fluid.rhoGS = 0.1;
model.fluid.rhoOS = 999.0140;


% Initialize bacterial concentration in the state
state0.nbact = zeros(model.G.cells.num, 1);  % Initial bacterial concentration

%% Set Up Bio-Compositional Fluid Model
%--------------------------------------------------------------------------
% Define the compositional fluid model and initialize the simulation.
%--------------------------------------------------------------------------
compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
    {'H2O', 'H2', 'CO2', 'C1'});

% Define backend for ADI (Automatic Differentiation)
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);

% Set up the biochemistry model
arg = {model.G, model.rock, model.fluid, compFluid, true, diagonal_backend, ...
    'water', false, 'oil', true, 'gas', true, 'bacteriamodel', true, ...
    'bDiffusionEffect', false, 'moleculardiffusion', false, ...
    'liquidPhase', 'O', 'vaporPhase', 'G'};
model = BiochemistryModel(arg{:});
minComp = 1.0e-8;
% Set initial conditions
z0 = [0.84, 0.005, 0.005,0.15];  % Initial composition
schedule.control.W(1).components = [0.001, 0.798, 0.2, 0.001];  % Inject 95% H2, 5% CO2
schedule.control.W(2).components = [0.001, 0.798, 0.2, 0.001];  % Same for second well
schedule.control.W(1).compi = [0,1];
schedule.control.W(2).compi = [0,1];
% schedule.control.W(1).val = 70;
% schedule.control.W(2).val = 40;
T = 40 + 273.15;  % Temperature (K)
p = 82 * barsa;  % Pressure (Pa)
nbact0 = 1.0e6;  % Initial bacterial concentration
state0 = initCompositionalStateBacteria(model, p, T, [0, 1], z0, nbact0, model.EOSModel);

%% Simulate the Schedule
%--------------------------------------------------------------------------
% Pack and simulate the problem.
%--------------------------------------------------------------------------
problem = packSimulationProblem(state0, model, schedule, 'simple_comp_SW_bact_clogging_test', 'name', name);
simulatePackedProblem(problem, 'restartStep', 1);

% Retrieve simulation results
[ws, states, rep] = getPackedSimulatorOutput(problem);


% %% Simulate without clogging effects
BaseName = 'simple_comp_SW_bact_noclogging';
% Set up the biochemistry model
cp = 0;
model.rock.perm = perm0;
model.rock.poro = poro0;
model.fluid.pvMultR =  @(p, nbact) 1;
arg = {model.G, model.rock, model.fluid, compFluid, true, diagonal_backend, ...
    'water', false, 'oil', true, 'gas', true, 'bacteriamodel', true, ...
    'bDiffusionEffect', false, 'moleculardiffusion', false, ...
    'liquidPhase', 'O', 'vaporPhase', 'G'};
model = BiochemistryModel(arg{:});
problemNoClog = packSimulationProblem(state0, model, schedule, BaseName, 'name', name);
simulatePackedProblem(problemNoClog, 'restartStep', 1);

[wsNoClog,statesNoClog] = getPackedSimulatorOutput(problemNoClog);

%% Simulate without bacterial effects
BaseName = 'simple_comp_SW_nobact';
arg = {model.G, model.rock, model.fluid, compFluid, true, diagonal_backend, ...
    'water', false, 'oil', true, 'gas', true, 'bacteriamodel', false, ...
    'bDiffusionEffect', false, 'moleculardiffusion', false, ...
    'liquidPhase', 'O', 'vaporPhase', 'G'};
model = BiochemistryModel(arg{:});
state0.nbact = 0;
problemNoBact = packSimulationProblem(state0, model, schedule, BaseName, 'name', name);
problemNoBact.SimulatorSetup.model.bacteriamodel = false;
simulatePackedProblem(problemNoBact, 'restartStep', 1);
[wsNoBact,statesNoBact] = getPackedSimulatorOutput(problemNoBact);

%% Compare Three Simulations: H2 Loss, CO2 Consumption, and CH4 Production
%--------------------------------------------------------------------------
% This section compares the results of three simulations:
%   1. With bacterial effects (bio-clogging).
%   2. Without clogging effects.
%   3. Without bacterial effects.
% The comparison includes:
%   - H2 loss percentage
%   - CO2 consumption percentage
%   - CH4 production percentage
% Results are displayed in a table for easy comparison.
%--------------------------------------------------------------------------

% Get component indices
componentNames = model.EOSModel.getComponentNames();
idxH2  = find(strcmp(componentNames, 'H2'));   % Index of H2
idxCO2 = find(strcmp(componentNames, 'CO2'));  % Index of CO2
idxCH4 = find(strcmp(componentNames, 'C1'));   % Index of CH4

% Number of time steps
nTimeSteps = numel(states);

% Initialize arrays to store component masses
totalH2_bact = zeros(nTimeSteps, 1);    % H2 mass with bacterial effects
totalH2_noClog = zeros(nTimeSteps, 1);  % H2 mass without clogging effects
totalH2_noBact = zeros(nTimeSteps, 1);  % H2 mass without bacterial effects

totalCO2_bact = zeros(nTimeSteps, 1);   % CO2 mass with bacterial effects
totalCO2_noClog = zeros(nTimeSteps, 1); % CO2 mass without clogging effects
totalCO2_noBact = zeros(nTimeSteps, 1); % CO2 mass without bacterial effects

totalCH4_bact = zeros(nTimeSteps, 1);   % CH4 mass with bacterial effects
totalCH4_noClog = zeros(nTimeSteps, 1); % CH4 mass without clogging effects
totalCH4_noBact = zeros(nTimeSteps, 1); % CH4 mass without bacterial effects

% Loop through time steps to calculate component masses
for i = 1:nTimeSteps
    % With bacterial effects
    totalH2_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{idxH2});
    totalCO2_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{idxCO2});
    totalCH4_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{idxCH4});

    % Without clogging effects
    totalH2_noClog(i) = sum(statesNoClog{i}.FlowProps.ComponentTotalMass{idxH2});
    totalCO2_noClog(i) = sum(statesNoClog{i}.FlowProps.ComponentTotalMass{idxCO2});
    totalCH4_noClog(i) = sum(statesNoClog{i}.FlowProps.ComponentTotalMass{idxCH4});

    % Without bacterial effects
    totalH2_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{idxH2});
    totalCO2_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{idxCO2});
    totalCH4_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{idxCH4});
end

%% Calculate Percentage Changes
% H2 Loss Percentage
H2_loss_bact = ((totalH2_noBact - totalH2_bact) ./ totalH2_noBact) * 100;
H2_loss_noClog = ((totalH2_noBact - totalH2_noClog) ./ totalH2_noBact) * 100;

% CO2 Consumption Percentage
CO2_consumption_bact = ((totalCO2_noBact - totalCO2_bact) ./ totalCO2_noBact) * 100;
CO2_consumption_noClog = ((totalCO2_noBact - totalCO2_noClog) ./ totalCO2_noBact) * 100;

% CH4 Production Percentage
CH4_production_bact = ((totalCH4_bact - totalCH4_noBact) ./ totalCH4_noBact) * 100;
CH4_production_noClog = ((totalCH4_noClog - totalCH4_noBact) ./ totalCH4_noBact) * 100;

%% Create a Comparison Table
% Extract final values for the table
final_H2_loss_bact = H2_loss_bact(end);
final_H2_loss_noClog = H2_loss_noClog(end);

final_CO2_consumption_bact = CO2_consumption_bact(end);
final_CO2_consumption_noClog = CO2_consumption_noClog(end);

final_CH4_production_bact = CH4_production_bact(end);
final_CH4_production_noClog = CH4_production_noClog(end);

% Create a table for comparison
comparisonTable = table(...
    [final_H2_loss_bact; final_H2_loss_noClog], ...
    [final_CO2_consumption_bact; final_CO2_consumption_noClog], ...
    [final_CH4_production_bact; final_CH4_production_noClog], ...
    'VariableNames', {'H2_Loss_Percentage', 'CO2_Consumption_Percentage', 'CH4_Production_Percentage'}, ...
    'RowNames', {'With_Bacterial_Effects', 'Without_Clogging_Effects'});

% Display the table
disp('Comparison of Simulation Results:');
disp(comparisonTable);
%% Comparison Plots for Abiotic, Bio, and Bio-Clog Scenarios
%--------------------------------------------------------------------------
% This script compares simulation results from three scenarios:
%   1. Abiotic: No bacterial effects or clogging.
%   2. Bio: Includes bacterial effects but no clogging.
%   3. Bio-Clog: Includes both bacterial effects and clogging.
% The script generates plots for component mole fractions, bacterial
% concentration, pressure, and saturation profiles.
%--------------------------------------------------------------------------

% Get the number of components and their names
ncomp = model.EOSModel.getNumberOfComponents();
cnames = model.EOSModel.getComponentNames();

% Load precomputed states for the three scenarios
data = {statesNoBact, states, statesNoClog}; % Abiotic, Bio, Bio-Clog
n = min(cellfun(@numel, data)); % Number of time steps to plot
names = {'ABIOTIC', 'BIO', 'BIO-CLOG'}; % Scenario names
markers = {'-', '--', ':'}; % Line styles for each scenario
colors = lines(ncomp + 1); % Color scheme for components and bacteria

%% Plot Component Mole Fractions and Bacterial Concentration
%--------------------------------------------------------------------------
% Plot the mole fractions of each component and normalized bacterial
% concentration for all three scenarios.
%--------------------------------------------------------------------------
figure('Position', get(0, 'DefaultFigurePosition') + [0, 0, 350, 0]);
hold on;

% Generate legend labels
legendLabels = cell(numel(data) * ncomp, 1);
for i = 1:numel(data)
    for j = 1:ncomp
        legendLabels{(i-1)*ncomp + j} = [names{i}, ' ', cnames{j}];
    end
end

% Plot data for each time step
for step = 1:n
    cla; hold on;
    for i = 1:numel(data)
        s = data{i}{step}; % State at current time step
        comp = s.components; % Component mole fractions
        nbact = s.nbact; % Bacterial concentration
        
        % Ensure components are in matrix form
        if iscell(comp)
            comp = [comp{:}];
        end
        
        % Plot component mole fractions
        for j = 1:ncomp
            plot(comp(:, j), markers{i}, 'LineWidth', 2, 'Color', colors(j, :));
        end
        
        % Plot normalized bacterial concentration
        plot(nbact / max(nbact), markers{i}, 'LineWidth', 2, 'Color', colors(ncomp + 1, :));
    end
    
    % Add legend and labels
    legend(legendLabels, 'Location', 'north', 'NumColumns', 3);
    ylim([0, 1]);
    ylabel('Mole Fraction / Normalized Concentration');
    xlabel('Cell Index');
    title('Component Mole Fractions and Bacterial Concentration');
    drawnow;
end

%% Compare Pressure and Saturation Profiles
%--------------------------------------------------------------------------
% Plot normalized pressure and gas saturation profiles for the Abiotic and
% Bio scenarios.
%--------------------------------------------------------------------------
figure;
hold on;

% Plot data for each time step
for step = 1:n
    clf; hold on;
    handles = []; % Handles for legend
    
    for i = 1:numel(data) % Loop through all three scenarios
        s = data{i}{step}; % State at current time step
        marker = markers{i}; % Line style
        linewidth = 2; % Line width
        
        % Plot gas saturation
        hs = plot(s.s(:, 2), marker, 'LineWidth', linewidth, 'Color', colors(i, :));
        
        % Plot normalized pressure
        p = s.pressure / max(s.pressure);
        hp = plot(p, marker, 'LineWidth', linewidth, 'Color', colors(i, :));
        
        % Store handles for legend
        if i == 1
            handles = [hs; hp]; % Include saturation and pressure
        end
    end
    
    % Add legend and labels
    legend(handles, 'Gas Saturation (ABIOTIC)', 'Normalized Pressure (ABIOTIC)', ...
        'Gas Saturation (BIO)', 'Normalized Pressure (BIO)', ...
        'Gas Saturation (BIO-CLOG)', 'Normalized Pressure (BIO-CLOG)', ...
        'Location', 'northoutside', 'Orientation', 'horizontal');
    ylim([0, 1]);
    ylabel('Normalized Value');
    xlabel('Cell Index');
    title('Gas Saturation and Pressure Profiles: Abiotic, Bio, and Bio-Clog');
    drawnow;
end
%% Set up interactive plotting
% Finally we set up interactive plots to make it easy to look at the
% results from the different simulators.

mrstModule add mrst-gui
for i = 1:numel(data)
    figure;
    plotToolbar(model.G, data{i}, 'plot1d', true);
    title(names{i});
end
%% Plot Ternary Diagrams for Three Simulations
%--------------------------------------------------------------------------
% This script plots ternary diagrams for three simulations:
%   1. Abiotic: No bacterial effects or clogging.
%   2. Bio: Includes bacterial effects but no clogging.
%   3. Bio-Clog: Includes both bacterial effects and clogging.
% The ternary diagrams show the displacement lines for the components in
% each simulation, plotted horizontally for easy comparison.
%--------------------------------------------------------------------------

% Define mapping functions for ternary plot
mapx = @(x, y, z) (1/2) * (2*y + z) ./ (x + y + z); % Map x-coordinate
mapy = @(x, y, z) (sqrt(3)/2) * z ./ (x + y + z);   % Map y-coordinate

colors = parula(numel(states)); % Color scheme for time steps

% Create a figure with three horizontal subplots
figure;
set(gcf, 'Position', [100, 100, 1200, 400]); % Adjust figure size

% Loop through each scenario and plot ternary diagram
for scenario = 1:numel(data)
    subplot(1, 3, scenario); % Create subplot for current scenario
    hold on;
    
    % Plot ternary diagram outline
    plot([0, 0.5, 1, 0], [0, sqrt(3)/2, 0, 0], 'k', 'LineWidth', 1.5);
    
    % Plot displacement lines for current scenario
    for i = 1:20:numel(data{scenario})
        C = data{scenario}{i}.components; % Component mole fractions
        plot(mapx(C(:, 1), C(:, 2), C(:, 3)), mapy(C(:, 1), C(:, 2), C(:, 3)), ...
            '-', 'color', colors(i, :), 'LineWidth', 1.5);
    end
    
    % Add component labels
    text(0, 0, cnames{1}, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 12);
    text(1, 0, cnames{2}, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 12);
    text(0.5, sqrt(3)/2, cnames{3}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % Add 0.5 labels
    text(mapx(0.5, 0.5, 0), mapy(0.5, 0.5, 0), '0.5', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', 10);
    text(mapx(0, 0.5, 0.5), mapy(0, 0.5, 0.5), '0.5', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
    text(mapx(0.5, 0.0, 0.5), mapy(0.5, 0.0, 0.5), '0.5', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
    
    % Add title for subplot
    title(names{scenario}, 'FontSize', 14);
    axis off;
end
%%
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
