%% Calibrate Sleipner model 
% This example calibrates the Sleipner model (top surface elevations, rock
% property, CO2 entry rates, CO2 density) by minimizing the misfit between
% simulated and observed CO2 plume thicknesses.
%
% The calibration result obtained by running this script is only one
% possible solution that is found within a "family of optimals".
%
% For more details, see Nilsen et al. "Using sensitivities and
% vertical-equilibrium models for parameter estimation of CO2 injection
% models with application to Sleipner data", Energy Procedia, 2017.
% (Presented at the 13th International Conference on Greenhouse Gas Control
% and Technologies, GHGT-13, 14-18 November 2016, Lausanne, Switzerland)

mrstModule add co2lab-common co2lab-ve optimization
gravity('reset',[0 0 9.8])
gravity on

%% Set up model with a schedule of 12 years (1999-2010).
% An initial simulation will be run.
num_years = 12;
[ model, schedule, initState, states, wellSols ] = ...
    calibrateSleipnerSetup( 'num_years',num_years );

Gt = model.G;

schedule_base = schedule;
states_base   = states;
wellSols_base = wellSols;

%% Get 2001, 2004, 2006, 2010 plume heights data
assert(num_years == 12);
assert(numel(states) == 12);
clear newplumes
newplumes=cell(numel(states),1);
for ty=[3,6,8,12]
    newplumes{ty}.year = ty + 1998;
    newplumes{ty}.FH = getSleipnerPlumeHeights('year',newplumes{ty}.year);
    newplumes{ty}.h = newplumes{ty}.FH(Gt.cells.centroids(:,1), Gt.cells.centroids(:,2));
end

% Initial simulated heights vs. observed plume thicknesses
plotObsAndSim(Gt, newplumes, states_base)

%% Define parameters to optimize.
% The limits for grid elevation change (dzLims) refer to the
% dimensional unit of meters for which the initial grid can change. The
% limits for CO2 density, permeability, and porosity are relative 
% multipliers applied to the initial variable values. The limits for CO2
% entry rates refer to the (constant) multiplier applied to the
% initial 1999-2010 entry rates.

setup = struct('model', model, 'schedule', schedule, 'state0', initState);

params = addParameter([], setup, 'name', 'permx', ...
                      'type', 'multiplier', ...
                      'lumping', ones(Gt.cells.num, 1), ...
                      'relativeLimits', [0.5 12], 'scaling', 'linear');
params = addParameter(params, setup, 'name', 'porevolume', ...
                      'type', 'multiplier', ...
                      'lumping', ones(Gt.cells.num, 1), ...
                      'relativeLimits', [0.5 1.4])
params = addParameter(params, setup, 'name', 'rhoG', ...
                      'belongsTo', 'model', ...
                      'location', {'fluid', 'rhoG'}, ...
                      'setfun', @setrhoGfun, ...
                      'relativeLimits' , [0.5, 2], 'scaling', 'linear');
params = addParameter(params, setup, 'name', 'rate', ...
                      'type', 'multiplier', ...
                      'relativeLimits', [0.1 1.5]); % @@@
params = addParameter(params, setup, 'name', 'dz', ...
                      'belongsTo', 'model', ...
                      'location', {'G', 'dz'}, ...
                      'boxLims', [-5, 5], 'uniformLimits', true);

%% Define objective function and compute initial value

% Determine scaling by comparing with base result
v_base = matchToPlumeData(model, states_base, newplumes);
scaling = max(sum(abs(cell2mat(v_base))), 1);

% Initial parameter values
p0 = getScaledParameterVector(setup, params);

% Define an objective function p -> obj(p) that can be used in the 
% nonlinear optimization routine
objh = @(p) evaluatePlumeMatch(p, @matchToPlumeData, setup, params, newplumes, scaling);
                               

% Run optimization (parameter calibration)
[v, p_opt, history] = unitBoxBFGS(p0, objh);

% Recompute optimal solution

% Plot final simulated heights vs. observed plume thicknesses

%% Store results

%% Some final output

% Initial versus calibrated grid