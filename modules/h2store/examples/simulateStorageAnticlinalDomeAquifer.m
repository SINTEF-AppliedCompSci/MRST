%% Cyclic Storage of H₂ in Anticlinal Dome Aquifer
% This example demonstrates a numerical simulation of hydrogen storage 
% within a dome-shaped reservoir. It represents the third case discussed 
% in the paper: "Phase behavior and black-oil simulations of hydrogen 
% storage in saline aquifers" (Ahmed, E., et al., 2024).
%
% This simulation examines hydrogen behavior in saline aquifers under varying depth and pressure scenarios. 
% The current setup models a low-depth scenario at 600 meters with an approximate pressure of 10 bar. 
% Users can switch to a high-depth, high-pressure case at 1200 meters and 120 bar by adjusting the depth and pressure in the options.
%
% Available input files cover specific case studies: (for lwo depth case of the paper)
% - **DOME_RS.DATA**: Dissolved hydrogen in brine, with evaporated water interactions disabled.
% - **DOME_RSRV.DATA**: Dissolved hydrogen in brine, with evaporated water interactions enabled.
% - **DOME_RS_SALT.DAT**: Hydrogen interaction with brine at a molality of 4 mol/kg.
% - **DOME_IMMISC.DAT**: Hydrogen  in immiscible case.
% - **DOME_RS_SALT.DAT**: Hydrogen interaction with brine at a molality of 4 mol/kg.

%
% For high-depth scenarios, change with the PVT files ending in `_HDEPTH.inc.
%
% Readers are encouraged to examine these files to understand the model configurations for different aquifer scenarios.

% The PVT (Pressure-Volume-Temperature) data used in the simulations 
% was generated using the RK-EoS (Redlich-Kwong Equation of State) 
% model, as referenced in the same paper. The code to generate PVT 
% data is also available and can create data files for varying 
% salinity conditions.

mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui vemmech test-suite
%% Notice on Computational Cost
warning('ComputationalCost:High', ...
       ['Please be advised that this example takes a long time ', ...
        'to run']);

name = 'H2_STORAGE_RSRV_DOME';
dataFile = fullfile(ROOTDIR,'..','modules','h2store','examples',  'data','AnticlinalDomeAquifer', 'DOME_RSRV.DATA');
deck = readEclipseDeck(dataFile);
%% Initialize the dome case simulation with specified parameters
[~, options, state0, model, schedule, deck,~] = modelForDome3DAquifer(deck);

%% Plot rock parameters
figure;
% porosity
subplot(1, 2, 1);
plotCellData(model.G, model.rock.poro, model.G.cells.centroids(:,2) > 2500, 'LineStyle', 'none');
colorbar;
title('Porosity Distribution');
xlabel('X [m]');
ylabel('Y [m]');
view(3);
% permeability
subplot(1, 2, 2);
plotCellData(model.G, log10(model.rock.perm(:,1)), model.G.cells.centroids(:,2) > 2500, 'LineStyle', 'none');
colorbar;
title('Permeability Distribution  (log10)');
xlabel('X [m]');
ylabel('Y [m]');
view(3);
sgtitle('Porosity and Permeability Distribution in the Dome Aquifer');

%% Add custom output functions to the model for analysis
model.OutputStateFunctions{end+1} = 'CapillaryPressure';
model.OutputStateFunctions{end+1} = 'SurfaceDensity';
model.OutputStateFunctions{end+1} = 'ShrinkageFactors';
model.outputFluxes = false;

%% Setup the non-linear solver and select the linear solver
nls = NonLinearSolver();
lsolve = selectLinearSolverAD(model);
nls.LinearSolver = lsolve;

%% Pack the simulation problem for execution
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);

%% Run the simulation
simulatePackedProblem(problem, 'restartStep',1);
%% gGet reservoir and well states
[ws,states] = getPackedSimulatorOutput(problem);
%% Plot well states
plotWellSols(ws);
%%
figure; plotToolbar(model.G, states); title(['MRST ', name]);
%% Compute total gas volume in the reservoir
nstep = numel(schedule.step.val);
dt = schedule.step.val;
time = cumsum(dt);
timeSteps = (1:nstep)';
[totGV_all, vH2InW_all, vH2InG_all] = computeGasVolumes(model, states, timeSteps);

%% Initialize variables for hydrogen volume, pressure, and rates
hydrogenLossPercentage = zeros(1, nstep);
totalGasVolume = zeros(1, nstep);
avPressure_res = zeros(1, nstep);
avSat_res = zeros(1, nstep);
bhp_res = zeros(1, nstep);
qGs_res = zeros(1, nstep);
qWs_res = zeros(1, nstep);

%% Compute cumulative hydrogen volume and averaged pressure over time
for step = 1:nstep
    totalGasVolume(step) = sum(totGV_all{step});
    dissolvedHydrogenVolume = sum(vH2InW_all{step});
    avPressure_res(step) = mean(states{step}.pressure);
    avSat_res(step) = mean(states{step}.s(:,2));
    
    if totalGasVolume(step) > 0
        hydrogenLossPercentage(step) = (dissolvedHydrogenVolume / totalGasVolume(step)) * 100;
    else
        hydrogenLossPercentage(step) = 0;
    end
    
    bhp_res(step) = ws{step}.bhp;
    qGs_res(step) = ws{step}.qGs;
    qWs_res(step) = ws{step}.qOs;
end

%% Plotting results over time (in months)
timeMonths = convertTo(time, 30.*day);

figure;
subplot(4,1,1);
plot(timeMonths, avSat_res, 'LineWidth', 1.5);
ylabel('Average Saturation');
title('Average Saturation Over Time');

subplot(4,1,2);
plot(timeMonths, avPressure_res, 'LineWidth', 1.5);
ylabel('Average Pressure (Pa)');
title('Average Pressure Over Time');

subplot(4,1,3);
plot(timeMonths, totalGasVolume, 'LineWidth', 1.5);
ylabel('Total Gas Volume');
title('Total Gas Volume Over Time');

subplot(4,1,4);
plot(timeMonths, cellfun(@sum, vH2InW_all), 'LineWidth', 1.5);
ylabel('Dissolved H2 Volume');
title('Dissolved H2 Volume Over Time');
xlabel('Time (months)');

figure;
subplot(4,1,1);
plot(timeMonths, bhp_res, 'LineWidth', 1.5);
ylabel('BHP (Pa)');
title('Bottom Hole Pressure Over Time');

subplot(4,1,2);
plot(timeMonths, qGs_res, 'LineWidth', 1.5);
ylabel('Gas Rate');
title('Gas Rate Over Time');

subplot(4,1,3);
plot(timeMonths, qWs_res, 'LineWidth', 1.5);
ylabel('Water Rate');
title('Water Rate Over Time');

subplot(4,1,4);
plot(timeMonths, hydrogenLossPercentage, 'LineWidth', 1.5);
ylabel('H2 Loss (%)');
title('Hydrogen Loss Percentage Over Time');
xlabel('Time (months)');
%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.
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
