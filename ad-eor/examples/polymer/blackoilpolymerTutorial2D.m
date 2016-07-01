% 2D Three-Phase Polymer Injection Case
% This example contains a simple 4000 m × 200 m × 125 m reservoir described
% on 20 × 1 × 5 uniform Cartesian grid. One injection well is located in
% the bottom two layers and one production well is located in the top two
% layers. Hydrostatic equilibration is used for initialization.

% The polymer injection schedule follows a typical polymer waterflooding
% strategy. The flooding process begins with primary waterflooding for 1260
% days, followed by a polymer slug injected over 1700 days, and then
% switching back to water injection. The injection well is under rate
% control with target rate 1000 m3 /day and upper limit of 450 bar on the
% bottom-hole pressure (bhp), whereas the production well is under pressure
% control with target bottom-home pressure 260 bar.

mrstModule add ad-core ad-blackoil ad-eor ad-fi ad-props deckformat
close all;

%% Set up model and initial conditions
% The data required for the example
% If the data does not exist locally, download it automatically
% The following are all the files needed for this tutorial
% The first two files are the data for a simulation with shear-thinning
% effects. The second two fils are the data for a simulation without shear
% effects. The last two are the reference results from Eclipse.
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

if ~all(e),
    pl = ''; if sum(e) ~= 1, pl = 's'; end
    msg = sprintf('Missing data file%s\n', pl);
    msg = [msg, sprintf('  * %s\n', fname{~e})];
    error('Dataset:Incomplete', msg);
end

% Parsing the data file with shear-thinning effects.
deck = readEclipseDeck(files{1});
% The deck is using metric system, MRST uses SI unit internally
deck = convertDeckUnits(deck);

% Process the grid
G = initEclipseGrid(deck);
G = computeGeometry(G);

% rock properties, including permeability, porosity, compressibility
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% fluid properties, such as densities, viscosities, relative
% permeability, etc.
fluid          = initDeckADIFluid(deck);

% Constructing the physical model used for this simulation
% Here, we use three phase blackoil-polymer model
model = ThreePhaseBlackOilPolymerModel(G, rock, fluid, 'inputdata', deck);

% Set initial equilibrium saturation and polymer related two variables
gravity on;
state0         = initEclipseState(G, deck, initEclipseFluid(deck));
% polymer concentration
state0.c       = zeros([G.cells.num, 1]);
% maximum polymer concentration, used to handle the polymer adsorption
state0.cmax    = zeros([G.cells.num, 1]);

% Convert the deck schedule into a MRST schedule by parsing the wells
schedule = convertDeckScheduleToMRST(model, deck);

%% Select nonlinear and linear solvers
% Use the AGMG algebraic multigrid solver if this is present, and if not,
% use the default MATLAB solver, which is a direct solver based on LU
% factorization. For this case, which is small in term of numbers of cells
% and unknowns, the direct solver is efficient enough, while for larger
% cases, the AGMG will improve the solution speed and be required.

% Using physically normalized residuals
model.useCNVConvergence = true;

% Choose the linear solver for pressure equation
if ~isempty(mrstPath('agmg'))
    mrstModule add agmg
    pSolver = AGMGSolverAD();
else
    pSolver = BackslashSolverAD();
end

% set up the CPR linear solver
linsolve = CPRSolverAD('ellipticSolver', pSolver);
nonlinearsolver = getNonLinearSolver(model, 'DynamicTimesteps', false, ...
                                     'LinearSolver', linsolve);
nonlinearsolver.useRelaxation = true;

%% To check the fluid properties interactively

fluidPlotPanelAD(model);


%% TODO: plotting the polymer properties

%% Plotting the grid, wells, and porosity and permeability
h=figure(); clf
set(h, 'Position', [100, 100, 900, 600]);

subplot(2,2,1);
plotGrid(G, 'FaceColor', 'none');
% plot wells, blue for injector, red for producer.
W = schedule.control(1).W; sgn = [W.sign];
plotWell(G, W(sgn>0), 'color', 'b')
plotWell(G, W(sgn<0), 'color', 'r')
axis tight off; view(20,8);

prop = {'PERMX', 'PERMZ', 'porosity'};
for i = 1:3
    if strcmp(prop{i},'PERMX')
        s_plot = rock.perm(:,1)/(milli*darcy);
    end
    if strcmp(prop{i},'PERMZ')
        s_plot = rock.perm(:,3)/(milli*darcy);
    end
    if strcmp(prop{i},'porosity')
        s_plot = rock.poro;
    end
    subplot(2, 2, i+1);
    plotCellData(G, s_plot, 'EdgeColor','k');
    title(prop{i});
    axis tight off; view(20,8); colorbar;
end
pause(0.1);

%% Plotting the inital saturations and pressures
h=figure(); clf
set(h, 'Position', [100, 100, 900, 600]);
% plotting the initial saturations
phases = {'water','oil','gas'};
for i=1:3
    subplot(2,2,i+1);
    plotCellData(G, state0.s(:,i), 'EdgeColor','k');
    title(['Initial ', phases{i}, ' saturation'])
    axis tight off; view(20,8); caxis([0, 1]); colorbar;
end
% plotting the inital pressure
subplot(2, 2, 1);
plotCellData(G, state0.pressure/barsa, 'EdgeColor','k');
title('Inital pressure (Bar)');
axis tight off; view(20,8); colorbar;
pause(0.1);

%% Run the schedule with plotting function
% Once a system has been created it is trivial to run the schedule. Any
% options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

% The AD-solvers allow for dyanmic plotting during the simulation process.
% We set up the following function to plot the evolution of the related
% variables (s:2 means oil saturation by default), the change of the well
% curves, and the a panel showing simulation progress. You can customize
% the function based on your own preference.

fn = getPlotAfterStep(state0, model, schedule, ...
    'plotWell', true, 'plotReservoir', true, 'view', [20, 8], ...
    'field', 's:2');

[wellSols, states, reports] = ...
    simulateScheduleAD(state0, model, schedule, ...
                    'NonLinearSolver', nonlinearsolver, 'afterStepFn', fn);

%% Comparing results with the reference results from commercial simualtor.
% loading the reference result smary.mat
load (files{5});
% the time for the reference result
T_ref = smry.get(':+:+:+:+', 'TIME', ':');
% the time for the MRST result
T_mrst = convertTo(cumsum(schedule.step.val), day);

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
% the first result of the commerical simualtor is always zero.
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
% the first result of the commerical simualtor is always zero.
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
% the first result of the commerical simualtor is always zero.
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
% the first result of the commerical simualtor is always zero.
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


%% Run the simulation without shear effects.
% You can load the files{3} to run the simulation.
% Here we just modify the model to disable the shear effects.

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
% the first result of the commerical simualtor is always zero.
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
% the first result of the commerical simualtor is always zero.
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
% the first result of the commerical simualtor is always zero.
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
% the first result of the commerical simualtor is always zero.
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
