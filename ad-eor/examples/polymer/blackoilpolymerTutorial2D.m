% 2D Three-Phase Polymer Injection Case
%
% This example contains a simple 4000 m-by-200 m-by-125 m reservoir
% described on 20-by-1-by-5 uniform Cartesian grid. One injection well is
% located in the bottom two layers and one production well is located in
% the top two layers. Hydrostatic equilibration is used for initialization.
%
% The polymer injection schedule follows a typical polymer waterflooding
% strategy. The flooding process begins with primary waterflooding for 1260
% days, followed by a polymer slug injected over 1700 days, and then
% switching back to water injection. The injection well is under rate
% control with target rate 1000 m3/day and upper limit of 450 bar on the
% bottom-hole pressure (bhp), whereas the production well is under pressure
% control with target bottom-home pressure 260 bar.

mrstModule add ad-core ad-blackoil ad-eor ad-props ...
               deckformat mrst-gui

%% Set up model and initial conditions
% The data required for the example
% If the data does not exist locally, download it automatically
% The following are all the files needed for this tutorial
% The first two files are the data for a simulation with shear-thinning
% effect. The second two fils are the data for a simulation without shear
% effect. The last two are the reference results from Eclipse.
fname = {'BOPOLYMER.DATA', ...
         'POLY.inc', ...
         'BOPOLYMER_NOSHEAR.DATA', ...
         'POLY_NOSHEAR.inc', ...
         'smry.mat', ...
         'smry_noshear.mat'};
files = fullfile(getDatasetPath('BlackoilPolymer2D', 'download', true),...
                                fname);

% check to make sure the files are complete
e = cellfun(@(pth) exist(pth, 'file') == 2, files);

if ~all(e)
    pl = ''; if sum(e) ~= 1, pl = 's'; end
    msg = sprintf('Missing data file%s\n', pl);
    msg = [msg, sprintf('  * %s\n', fname{~e})];
    error('Dataset:Incomplete', msg);
end
gravity reset on;

% Parsing the data file with shear-thinning effect.
deck = readEclipseDeck(files{1});
% The deck is using metric system, MRST uses SI unit internally
deck = convertDeckUnits(deck);

% Construct physical model, initial state and dynamic well controls.
[state0, model, schedule] = ...
   initEclipseProblemAD(deck, 'UseLegacyModels', true);

% Add polymer concentration
state0.cp   = zeros([model.G.cells.num, 1]);

% maximum polymer concentration, used to handle the polymer adsorption
state0.cpmax= zeros([model.G.cells.num, 1]);

%% Select nonlinear and linear solvers

% Using physically normalized residuals for non-linear convergence
% calcuation.
model.useCNVConvergence = true;

% Setting up the non-linear solver.
nonlinearsolver = NonLinearSolver();
nonlinearsolver.useRelaxation = true;

%% Visualize the properties of the black-oil fluid model
% We launch the interactive viewer for the black-oil fluid model, and
% specify the pressure range for which the tables are given. Note that
% extrapolation beyond the specified values for black-oil properties can
% result in non-physical curves, depending on how the input was given.
inspectFluidModel(model, 'pressureRange', (50:10:600)*barsa)
example_name = 'blackoil2D';
vizPolymerModel();

%% Run the schedule with plotting function
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

% The AD-solvers allow for dyanmic plotting during the simulation process.
% We set up the following function to plot the evolution of the related
% variables (s:2 means oil saturation by default), the change of the well
% curves, and the a panel showing simulation progress. You can customize
% the function based on your own preference.
close all
fn = getPlotAfterStep(state0, model, schedule, ...
    'plotWell', true, 'plotReservoir', true, 'view', [20, 8], ...
    'field', 's:2');
[wellSols, states, reports] = ...
    simulateScheduleAD(state0, model, schedule, ...
                    'NonLinearSolver', nonlinearsolver, 'afterStepFn', fn);

%% Comparing the result with reference result from commercial simualtor.
% loading the reference result smary.mat
load (files{5});
% the time for the reference result
T_ref = smry.get(':+:+:+:+', 'TIME', ':');
% the time for the MRST result
T_mrst = convertTo(cumsum(schedule.step.val), day);

% generate a color map for plotting use.
color_map = lines(10);
color_mrst = color_map(1, :);
color_ref  = color_map(2, :);

mrstplot = @(T_mrst, data, color) plot(T_mrst, data, '-', ...
                                       'linewidth', 2, 'color', color);
referenceplot = @(T_ref, data, color) plot(T_ref, data, '--',...
                                       'linewidth', 4, 'color', color);

% Plotting the water injection rate
h = figure(); clf;
set(gca,'FontSize',20);
set(h, 'Position', [100, 100, 900, 600]);
well_name = 'INJE01';
reference = smry.get(well_name, 'WWIR', ':');
% the first value of the result of the commerical simualtor is always zero.
reference(1) = nan;
mrst = convertTo(abs(getWellOutput(wellSols, 'qWs', well_name)), ...
                                                              meter^3/day);
hold on;
mrstplot(T_mrst, mrst, color_mrst);
referenceplot(T_ref, reference, color_ref);
title(['Water injection rate for ', well_name]);
xlabel('Time (days)');
ylabel('Water injection rate (m^3/day)');
axis tight;
legend({'MRST', 'reference'})
pause(0.1);

% Plotting the bhp for the injection well
h = figure(); clf;
set(gca,'FontSize',20);
set(h, 'Position', [100, 100, 900, 600]);
well_name = 'INJE01';
reference = smry.get(well_name, 'WBHP', ':');
% the first value of the result of the commerical simualtor is always zero.
reference(1) = nan;
mrst = convertTo(abs(getWellOutput(wellSols, 'bhp', well_name)), barsa);
hold on;
mrstplot(T_mrst, mrst, color_mrst);
referenceplot(T_ref, reference, color_ref);
title(['Bottom hole pressure for ', well_name]);
xlabel('Time (days)');
ylabel('Bottom hole pressure (Bar)');
axis tight;
legend({'MRST', 'reference'})
pause(0.1);

% Plotting the oil production rate
h = figure(); clf;
set(gca,'FontSize',20);
set(h, 'Position', [100, 100, 900, 600]);
well_name = 'PROD01';
reference = smry.get(well_name, 'WOPR', ':');
% the first value of the result of the commerical simualtor is always zero.
reference(1) = nan;
mrst = convertTo(abs(getWellOutput(wellSols, 'qOs', well_name)), ...
                                                              meter^3/day);
hold on;
mrstplot(T_mrst, mrst, color_mrst);
referenceplot(T_ref, reference, color_ref);
title(['Oil production rate for ', well_name]);
xlabel('Time (days)');
ylabel('Oil production rate (m^3/day)');
axis tight;
legend({'MRST', 'reference'})
pause(0.1);

% Plotting the oil production rate
h = figure(); clf;
set(gca,'FontSize',20);
set(h, 'Position', [100, 100, 900, 600]);
well_name = 'PROD01';
reference = smry.get(well_name, 'WWPR', ':');
% the first value of the result of the commerical simualtor is always zero.
reference(1) = nan;
mrst = convertTo(abs(getWellOutput(wellSols, 'qWs', well_name)), ...
                                                             meter^3/day);
hold on;
mrstplot(T_mrst, mrst, color_mrst);
referenceplot(T_ref, reference, color_ref);
title(['Water production rate for ', well_name]);
xlabel('Time (days)');
ylabel('Water production rate (m^3/day)');
axis tight;
legend({'MRST', 'reference'})
pause(0.1);


%% Run the simulation without shear effect.
% You can load the files{3} to run the simulation.
% Here we just modify the model directly to disable the shear effect.
close all
model.usingShear = false;

fn = getPlotAfterStep(state0, model, schedule, ...
    'plotWell', true, 'plotReservoir', true, 'view', [20, 8], ...
    'field', 's:2');
[wellSolsNoShear, statesNoShear, reportsNoShear] = ...
    simulateScheduleAD(state0, model, schedule, ...
                    'NonLinearSolver', nonlinearsolver, 'afterStepFn', fn);

%% Plotting the results from two simulations and their reference results.
% load the reference results for the non-shear case.
load (files{6});
% the time for the reference result
% since the commercial software might cut the time steps, the actually used
% schedule can be different from the previous running with shear-thinning.
T_ref_noshear = smry.get(':+:+:+:+', 'TIME', ':');

color_mrst_noshear = color_map(3, :);
color_ref_noshear = color_map(4,:);

% Plotting the water injection rate
h = figure(); clf;
set(gca,'FontSize',20);
set(h, 'Position', [100, 100, 900, 600]);
well_name = 'INJE01';
reference = smry.get(well_name, 'WWIR', ':');
reference_noshear = smry_noshear.get(well_name, 'WWIR', ':');
% the first value of the result of the commerical simualtor is always zero.
reference(1) = nan; reference_noshear(1) = nan;

mrst = convertTo(abs(getWellOutput(wellSols, 'qWs', well_name)),...
                                                        meter^3/day);
mrst_noshear = convertTo(abs(getWellOutput(wellSolsNoShear, 'qWs', ...
                                          well_name)), meter^3/day);
hold on;
mrstplot(T_mrst, mrst, color_mrst);
referenceplot(T_ref, reference, color_ref);
mrstplot(T_mrst, mrst_noshear, color_mrst_noshear);
referenceplot(T_ref_noshear, reference_noshear, color_ref_noshear);
title(['Water injection rate for ', well_name]);
xlabel('Time (days)');
ylabel('Water injection rate (m^3/day)');
axis tight;
legend({'MRST', 'reference', 'MRST no shear', 'reference no shear'})
pause(0.1);


% Plotting the bhp for the injection well
h = figure(); clf;
set(gca,'FontSize',20);
set(h, 'Position', [100, 100, 900, 600]);
well_name = 'INJE01';
reference = smry.get(well_name, 'WBHP', ':');
reference_noshear = smry_noshear.get(well_name, 'WBHP', ':');
% the first value of the result of the commerical simualtor is always zero.
reference(1) = nan; reference_noshear(1) = nan;
mrst = convertTo(abs(getWellOutput(wellSols, 'bhp', well_name)), barsa);
mrst_noshear = convertTo(abs(getWellOutput(wellSolsNoShear, 'bhp', ...
                                           well_name)), barsa);
hold on;
mrstplot(T_mrst, mrst, color_mrst);
referenceplot(T_ref, reference, color_ref);
mrstplot(T_mrst, mrst_noshear, color_mrst_noshear);
referenceplot(T_ref_noshear, reference_noshear, color_ref_noshear);
title(['Bottom hole pressure for ', well_name]);
xlabel('Time (days)');
ylabel('Bottom hole pressure (Bar)');
axis tight;
legend({'MRST', 'reference', 'MRST no shear', 'reference no shear'})
pause(0.1);

% Plotting the oil production rate
h = figure(); clf;
set(gca,'FontSize',20);
set(h, 'Position', [100, 100, 900, 600]);
well_name = 'PROD01';
reference = smry.get(well_name, 'WOPR', ':');
reference_noshear = smry_noshear.get(well_name, 'WOPR', ':');
% the first value of the result of the commerical simualtor is always zero.
reference(1) = nan; reference_noshear(1) = nan;
mrst = convertTo(abs(getWellOutput(wellSols, 'qOs', well_name)), ...
                                                          meter^3/day);
mrst_noshear = convertTo(abs(getWellOutput(wellSolsNoShear, 'qOs', ...
                                          well_name)), meter^3/day);
hold on;
mrstplot(T_mrst, mrst, color_mrst);
referenceplot(T_ref, reference, color_ref);
mrstplot(T_mrst, mrst_noshear, color_mrst_noshear);
referenceplot(T_ref_noshear, reference_noshear, color_ref_noshear);
title(['Oil production rate for ', well_name]);
xlabel('Time (days)');
ylabel('Oil production rate (m^3/day)');
axis tight;
legend({'MRST', 'reference'})
pause(0.1);

% Plotting the water production rate
h = figure(); clf;
set(gca,'FontSize',20);
set(h, 'Position', [100, 100, 900, 600]);
well_name = 'PROD01';
reference = smry.get(well_name, 'WWPR', ':');
reference_noshear = smry_noshear.get(well_name, 'WWPR', ':');
% the first value of the result of the commerical simualtor is always zero.
reference(1) = nan; reference_noshear(1) = nan;
mrst = convertTo(abs(getWellOutput(wellSols, 'qWs', well_name)), ...
                                                              meter^3/day);
mrst_noshear = convertTo(abs(getWellOutput(wellSolsNoShear, 'qWs', ...
                                          well_name)), meter^3/day);
hold on;
mrstplot(T_mrst, mrst, color_mrst);
referenceplot(T_ref, reference, color_ref);
mrstplot(T_mrst, mrst_noshear, color_mrst_noshear);
referenceplot(T_ref_noshear, reference_noshear, color_ref_noshear);
title(['Water production rate for ', well_name]);
xlabel('Time (days)');
ylabel('Water production rate (m^3/day)');
axis tight;
legend({'MRST', 'reference', 'MRST no shear', 'reference no shear'})
pause(0.1);


%%
% save resMRSTPolymer wellSols states schedule;
fprintf('The simulation has been finished! \n');

%% Copyright notice

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
