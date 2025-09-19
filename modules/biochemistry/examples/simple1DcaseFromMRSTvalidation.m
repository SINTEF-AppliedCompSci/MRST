%% 1D Compositional Simulation with Bio-Clogging - Enhanced Version
%==========================================================================
% This script simulates a 1D two-phase compositional problem using MRST with:
% - Injection of 95% H2 and 5% CO2 into a reservoir with CO2, Methane, H2
% - Bio-clogging effects (porosity and permeability reduction)
% - Soreide-Whitson EOS and LBC viscosity correlation
% - Relative permeability modifications for residual saturations
% - Comprehensive comparison of three scenarios:
%   1. With bacterial effects and clogging
%   2. With bacterial effects but no clogging
%   3. Without bacterial effects (abiotic)
%
% Enhanced features:
% - Better code organization with clear sections
% - Improved documentation
% - Additional validation checks
% - More efficient plotting
% - Better variable naming
%==========================================================================

clear all;
close all;

%% 1. Initialize MRST and Load Modules
%--------------------------------------------------------------------------
mrstModule add  biochemistry compositional deckformat ad-core ad-props mrst-gui;
mrstVerbose off;  % Reduce command window output for cleaner execution

%% 2. Set Up Base Simulation Model
%--------------------------------------------------------------------------
% Initialize from MRST validation example and adjust simulation parameters
[state0, model, schedule, ref] = setupSimpleCompositionalExample(false);

% Adjust time steps for faster simulation while maintaining resolution
schedule.step.val = sort(schedule.step.val/3);
schedule.step.val = schedule.step.val(1:400);
schedule.step.control = schedule.step.control(1:400);

simulationName = 'comp-bio-clogging';

%% 3. Define Fluid Properties and EOS
%--------------------------------------------------------------------------
% Use Soreide-Whitson EOS for fluid modeling
eosname = 'sw';
model.EOSModel = SoreideWhitsonEquationOfStateModel(...
    model.G, model.EOSModel.CompositionalMixture, eosname);
model.EOSModel.msalt = 5;  % Salt parameter for EOS

%% 4. Modify Relative Permeability
%--------------------------------------------------------------------------
% Adjust relperm curves to account for residual saturations
swc = 0.2;  % Connate water saturation
sgr = 0.1;  % Residual gas saturation
sor = 0.2;  % Residual oil saturation

% Create shifted saturation functions with residual saturations
model.fluid = modifyRelPermForResidualSaturations(model.fluid, swc, sgr, sor);

% Initial conditions
%initialComposition = [0.4, 0.445, 0.055, 0.1];  % High H2 content case
initialComposition = [0.9, 0.045, 0.005,0.05]; % Low H2 content
initialTemperature = 40 + 273.15;  % K
initialPressure = 82 * barsa;      % Pa
initialBacteria = 9;              % Initial bacterial concentration
%% 5. Define Bio-Clogging Effects
%--------------------------------------------------------------------------
% Porosity and permeability reduction due to bacterial growth
[model, poro0, perm0] = setupBioCloggingModel(model, initialBacteria);

%% 6. Set Up Bio-Compositional Fluid Model
%--------------------------------------------------------------------------
% Define components and initialize the biochemical model
compFluid = TableCompositionalMixture(...
    {'Water', 'Hydrogen', 'Methane', 'CarbonDioxide'}, ...
    {'H2O', 'H2', 'C1', 'CO2'});

% Set up ADI backend
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);

% Biochemical model arguments
modelArgs = {model.G, model.rock, model.fluid, compFluid, true, diagonal_backend, ...
    'water', false, 'oil', true, 'gas', true, 'bacteriamodel', true, ...
    'liquidPhase', 'O', 'vaporPhase', 'G'};

model = BiochemistryModel(modelArgs{:});
model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', false);

% Set well compositions
schedule.control.W(1).components = [0.001, 0.958, 0.001, 0.05];
schedule.control.W(2).components = [0.001, 0.998, 0.001, 0.05];
schedule.control.W(1).compi = [0, 1];
schedule.control.W(2).compi = [0, 1];
schedule.control.W(2).type = 'rate';
schedule.control.W(2).val = 0;
schedule.control.W(1).type = 'rate';
schedule.control.W(1).val = 0;

% Initialize state
state0 = initCompositionalStateBacteria(model, initialPressure, ...
    initialTemperature, [0, 1], initialComposition, initialBacteria, model.EOSModel);

%% 7. Run Simulations for Three Scenarios
%--------------------------------------------------------------------------
% Scenario 1: With bacterial effects and clogging
problemBioClog = packSimulationProblem(state0, model, schedule, ...
    'bio_clogging', 'name', simulationName);
problemBioClog.SimulatorSetup.model.OutputStateFunctions{end} = 'ComponentPhaseMass';
nls = NonLinearSolver();
problemBioClog.SimulatorSetup.NonLinearSolver = nls;
simulatePackedProblem(problemBioClog, 'restartStep', 1);
[wsBioClog, statesBioClog, repBioClog] = getPackedSimulatorOutput(problemBioClog);

% Scenario 2: With bacterial effects but no clogging
modelNoClog = model;
modelNoClog.rock.perm = model.rock.perm(1,state0.nbact);
modelNoClog.rock.poro = model.rock.poro(1,state0.nbact);
modelNoClog.fluid.pvMultR = @(p, nbact) 1;  % Remove clogging effect
problemBioNoClog = packSimulationProblem(state0, modelNoClog, schedule, ...
    'bio_no_clogging', 'name', simulationName);
problemBioNoClog.SimulatorSetup.model.OutputStateFunctions{end} = 'ComponentPhaseMass';
simulatePackedProblem(problemBioNoClog,'restartStep', 1);
[wsBioNoClog, statesBioNoClog] = getPackedSimulatorOutput(problemBioNoClog);

% Scenario 3: Without bacterial effects (abiotic)
% Biochemical model arguments
modelNoBact = modelNoClog;
modelNoBact.bacteriamodel = false;   
state0NoBact = state0;
problemNoBact = packSimulationProblem(state0NoBact, modelNoBact, schedule, ...
    'no_bacteria', 'name', simulationName);
problemNoBact.SimulatorSetup.model.OutputStateFunctions{end} = 'ComponentPhaseMass';
simulatePackedProblem(problemNoBact,'restartStep', 1);
[wsNoBact, statesNoBact] = getPackedSimulatorOutput(problemNoBact);

%% 8. Analyze and Compare Results
%--------------------------------------------------------------------------
% Get component indices
componentNames = model.EOSModel.getComponentNames();
idxH2 = find(strcmp(componentNames, 'H2'));
idxCO2 = find(strcmp(componentNames, 'CO2'));
idxCH4 = find(strcmp(componentNames, 'C1'));

% Prepare data for analysis
scenarios = {
    struct('name', 'Bio-Clog', 'model', model,'states', {statesBioClog}, 'color', [0 0.5 0], 'line', '-'), 
    struct('name', 'Bio-NoClog','model',modelNoClog, 'states', {statesBioNoClog}, 'color', [0 0.7 0.7], 'line', '--'), 
    struct('name', 'Abiotic', 'model', modelNoBact,'states', {statesNoBact}, 'color', [0.7 0 0], 'line', ':')
};

% Time information
timeDays = cumsum(schedule.step.val)/day;
timeYears = timeDays/365;

%% 9. Plot Results
%--------------------------------------------------------------------------
% 9.1 Plot component mole fractions over time
plotComponentProfiles(scenarios, componentNames, timeYears);

% 9.2 Plot bacterial concentration and clogging effects
plotBacterialEffects(scenarios, timeYears, poro0, perm0);

% 9.3 Plot pressure and saturation profiles
plotPressureAndSaturation(scenarios, timeYears);

% 9.4 Plot comparison of H2 loss, CO2 consumption, and CH4 production
plotComponentComparison(scenarios, timeYears, idxH2, idxCO2, idxCH4);

% 9.5 Interactive plotting
plotToolbar(model.G, statesBioClog, 'plot1d', true, 'field', 'pressure');
title('Bio-Clogging Scenario - Pressure');

%% Helper Functions
%--------------------------------------------------------------------------

function fluid = modifyRelPermForResidualSaturations(fluid, swc, sgr, sor)
    % Modify relative permeability functions for residual saturations
    
    % Define synthetic saturation ranges
    SW = linspace(0, 1, 100);
    SO = linspace(0, 1, 100);
    SG = linspace(0, 1, 100);
    
    % Evaluate original relative permeability functions
    krW_original = arrayfun(fluid.krW, SW);
    krOW_original = arrayfun(fluid.krOW, SO);
    krG_original = arrayfun(fluid.krG, SG);
    krOG_original = arrayfun(fluid.krOG, SO);
    
    % Update water relative permeability
    SW_shifted = max(SW - swc, 0);
    krW_new = interp1(SW, krW_original, SW_shifted, 'linear', 0);
    fluid.krW = @(sw) interp1(SW, krW_new, value(max(sw - swc, 0)));
    
    % Update gas relative permeability
    SG_shifted = max(SG - sgr, 0);
    krG_new = interp1(SG, krG_original, SG_shifted, 'linear', 0);
    fluid.krG = @(sg) interp1(SG, krG_new, value(max(sg - sgr, 0)));
    
    % Update oil relative permeability
    SO_shifted = max(SO - sor, 0);
    krOW_new = interp1(SO, krOW_original, SO_shifted, 'linear', 0);
    krOG_new = interp1(SO, krOG_original, SO_shifted, 'linear', 0);
    fluid.krOW = @(so) interp1(SO, krOW_new, value(max(so - sor, 0)));
    fluid.krOG = @(so) interp1(SO, krOG_new, value(max(so - sor, 0)));
    
    % Define relative permeability endpoints
    fluid.krPts.w = [swc, 0.9];
    fluid.krPts.g = [sgr, 0.8];
    fluid.krPts.ow = [0, 1 - swc];
    fluid.krPts.og = [0, 1 - sgr];
end

function [model, poro0, perm0] = setupBioCloggingModel(model, nbact0)
    % Set up bio-clogging model with porosity and permeability reduction
    
    % Initial rock properties
    poro0 = model.rock.poro;
    perm0 = model.rock.perm(:, 1);
    poro0(1) = 0.9;
    poro0(end) = 0.9;
    
    % Bacterial concentration parameters
    nc = 180;  % Critical bacterial concentration
    cp = 1.0;   % Clogging coefficient
    
    % Define porosity multiplier
        
    scale = 1 + 0.*(nbact0 / nc).^2;
    pvMult_nbact = @(nbact) 1 ./ (1 + cp.*(nbact/nc).^2);
    model.fluid.pvMultR = @(p, nbact) pvMult_nbact(nbact);
    poro = @(p, nbact) scale.*poro0 .* pvMult_nbact(nbact);
    model.rock.poro = poro;
    
    % Define permeability as a function of porosity
    tau = @(p, nbact) ((1 - poro0) ./ (1 - poro(p, nbact))).^2 .* (poro(p, nbact) ./ poro0).^3;
    perm = @(p, nbact) perm0 .* tau(p, nbact);
    model.rock.perm = perm;
    
    % Define densities
    model.fluid.rhoGS = 0.08988;
    model.fluid.rhoOS = 999.0140;
end

function plotComponentProfiles(scenarios, componentNames, timeYears)
    % Plot component mole fractions for all scenarios
    
    figure('Position', [100, 100, 1200, 800]);
    nComp = numel(componentNames);
    
    for compIdx = 1:nComp
        subplot(2, 2, compIdx);
        hold on;
        
        for scenIdx = 1:numel(scenarios)
            scen = scenarios{scenIdx};
            nSteps = length(timeYears);
            compData = zeros(nSteps, 1);
            
            for step = 1:nSteps
                if iscell(scen.states{step}.components)
                    compData(step) = mean([scen.states{step}.components{:, compIdx}]);
                else
                    compData(step) = mean(scen.states{step}.components(:, compIdx));
                end
            end
            
            plot(timeYears, compData, scen.line, ...
                'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
        end
        
        title(componentNames{compIdx});
        xlabel('Time (years)');
        ylabel('Average Mole Fraction');
        grid on;
        legend('Location', 'best');
    end
    
    sgtitle('Component Mole Fraction Evolution');
end

function plotBacterialEffects(scenarios, timeYears, poro0, perm0)
    % Plot bacterial concentration and clogging effects
    
    figure('Position', [100, 100, 1000, 800]);
    
    % Find which scenario has bacteria
    hasBacteria = cellfun(@(x) isfield(x.states{1}, 'nbact'), scenarios);
    bioScenarios = scenarios(hasBacteria);
    
    if isempty(bioScenarios)
        return;  % No bacterial scenarios to plot
    end
    
    % Plot bacterial concentration
    subplot(2, 2, 1);
    hold on;
    for scenIdx = 1:numel(bioScenarios)
        scen = bioScenarios{scenIdx};
        nSteps = numel(scen.states);
        bactData = zeros(nSteps, 1);
        
        for step = 1:nSteps
            bactData(step) = mean(scen.states{step}.nbact);
        end
        
        plot(timeYears, bactData, scen.line, ...
            'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
    end
    title('Average Bacterial Concentration');
    xlabel('Time (years)');
    ylabel('Concentration');
    grid on;
    legend;
    
    % Plot porosity reduction
    subplot(2, 2, 2);
    hold on;
    for scenIdx = 1:numel(bioScenarios)
        scen = bioScenarios{scenIdx};
        nSteps = numel(scen.states);
        poroData = zeros(nSteps, 1);
        
        for step = 1:nSteps
            currentPoro = scen.model.rock.poro;
            if isa(currentPoro, 'function_handle')
                nbact = scen.states{step}.nbact;
                poroData(step) = mean(currentPoro(1, nbact));
            else
                poroData(step) = mean(currentPoro);
            end
        end
        
        plot(timeYears, poroData./mean(poro0), scen.line, ...
            'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
    end
    title('Porosity Reduction (Normalized)');
    xlabel('Time (years)');
    ylabel('Porosity / Initial Porosity');
    grid on;
    
    % Plot permeability reduction
    subplot(2, 2, 3);
    hold on;
    for scenIdx = 1:numel(bioScenarios)
        scen = bioScenarios{scenIdx};
        nSteps = numel(scen.states);
        permData = zeros(nSteps, 1);
        
        for step = 1:nSteps
            currentPerm = scen.model.rock.perm;
            if isa(currentPerm, 'function_handle')
                nbact = scen.states{step}.nbact;
                permData(step) = mean(currentPerm(1, nbact));
            else
                permData(step) = mean(currentPerm(:,1));
            end
        end
        
        plot(timeYears, permData./mean(perm0), scen.line, ...
            'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
    end
    title('Permeability Reduction (Normalized)');
    xlabel('Time (years)');
    ylabel('Permeability / Initial Permeability');
    grid on;
    
    sgtitle('Bio-Clogging Effects');
end

function plotPressureAndSaturation(scenarios, timeYears)
    % Plot pressure and saturation profiles
    
    figure('Position', [100, 100, 1200, 500]);
    
    % Plot pressure
    subplot(1, 2, 1);
    hold on;
    for scenIdx = 1:numel(scenarios)
        scen = scenarios{scenIdx};
        nSteps = numel(scen.states);
        pressureData = zeros(nSteps, 1);
        
        for step = 1:nSteps
            pressureData(step) = mean(scen.states{step}.pressure)/barsa;
        end
        
        plot(timeYears, pressureData, scen.line, ...
            'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
    end
    title('Average Pressure');
    xlabel('Time (years)');
    ylabel('Pressure (bar)');
    grid on;
    legend;
    
    % Plot gas saturation
    subplot(1, 2, 2);
    hold on;
    for scenIdx = 1:numel(scenarios)
        scen = scenarios{scenIdx};
        nSteps = numel(scen.states);
        satData = zeros(nSteps, 1);
        
        for step = 1:nSteps
            satData(step) = mean(scen.states{step}.s(:,2)); % Gas saturation
        end
        
        plot(timeYears, satData, scen.line, ...
            'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
    end
    title('Average Gas Saturation');
    xlabel('Time (years)');
    ylabel('Gas Saturation');
    grid on;
    
    sgtitle('Pressure and Saturation Profiles');
end

function plotComponentComparison(scenarios, timeYears, idxH2, idxCO2, idxCH4)
    % Plot comparison of H2 loss, CO2 consumption, and CH4 production
    
    % Get abiotic scenario as reference
    abioticIdx = find(strcmp(cellfun(@(x) x.name, scenarios, 'UniformOutput', false), 'Abiotic'));
    if isempty(abioticIdx)
        return;  % Need abiotic scenario for comparison
    end
    
    abioticStates = scenarios{abioticIdx}.states;
    nSteps = numel(abioticStates);
    
    % Initialize arrays
    totalH2_abiotic = zeros(nSteps, 1);
    totalCO2_abiotic = zeros(nSteps, 1);
    totalCH4_abiotic = zeros(nSteps, 1);
    
    % Calculate total masses for abiotic scenario
    for step = 1:nSteps
        totalH2_abiotic(step) = sum(abioticStates{step}.FlowProps.ComponentTotalMass{idxH2});
        totalCO2_abiotic(step) = sum(abioticStates{step}.FlowProps.ComponentTotalMass{idxCO2});
        totalCH4_abiotic(step) = sum(abioticStates{step}.FlowProps.ComponentTotalMass{idxCH4});
    end
    
    % Create comparison plots
    figure('Position', [100, 100, 1200, 400]);
    
    % H2 Loss
    subplot(1, 3, 1);
    hold on;
    for scenIdx = 1:numel(scenarios)
        if scenIdx == abioticIdx, continue; end  % Skip abiotic scenario
        
        scen = scenarios{scenIdx};
        totalH2 = zeros(nSteps, 1);
        
        for step = 1:nSteps
            totalH2(step) = sum(scen.states{step}.FlowProps.ComponentTotalMass{idxH2});
        end
        
        H2_loss = ((totalH2_abiotic - totalH2) ./ totalH2_abiotic) * 100;
        plot(timeYears, H2_loss, scen.line, ...
            'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
    end
    title('H2 Loss Percentage');
    xlabel('Time (years)');
    ylabel('H2 Loss (%)');
    grid on;
    legend;
    
    % CO2 Consumption
    subplot(1, 3, 2);
    hold on;
    for scenIdx = 1:numel(scenarios)
        if scenIdx == abioticIdx, continue; end  % Skip abiotic scenario
        
        scen = scenarios{scenIdx};
        totalCO2 = zeros(nSteps, 1);
        
        for step = 1:nSteps
            totalCO2(step) = sum(scen.states{step}.FlowProps.ComponentTotalMass{idxCO2});
        end
        
        CO2_consumption = ((totalCO2_abiotic - totalCO2) ./ totalCO2_abiotic) * 100;
        plot(timeYears, CO2_consumption, scen.line, ...
            'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
    end
    title('CO2 Consumption Percentage');
    xlabel('Time (years)');
    ylabel('CO2 Consumption (%)');
    grid on;
    
    % CH4 Production
    subplot(1, 3, 3);
    hold on;
    for scenIdx = 1:numel(scenarios)
        if scenIdx == abioticIdx, continue; end  % Skip abiotic scenario
        
        scen = scenarios{scenIdx};
        totalCH4 = zeros(nSteps, 1);
        
        for step = 1:nSteps
            totalCH4(step) = sum(scen.states{step}.FlowProps.ComponentTotalMass{idxCH4});
        end
        
        CH4_production = ((totalCH4 - totalCH4_abiotic) ./ totalCH4_abiotic) * 100;
        plot(timeYears, CH4_production, scen.line, ...
            'Color', scen.color, 'LineWidth', 2, 'DisplayName', scen.name);
    end
    title('CH4 Production Percentage');
    xlabel('Time (years)');
    ylabel('CH4 Production (%)');
    grid on;
    
    sgtitle('Component Mass Comparison Relative to Abiotic Scenario');
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST. If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>