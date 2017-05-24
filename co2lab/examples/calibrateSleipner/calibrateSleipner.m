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

mrstModule add co2lab
gravity('reset',[0 0 9.8])
gravity on


%% Set up model with a schedule of 12 years (1999-2010).
% An initial simulation will be run.
num_years = 12;
[ smodel, schedule, initState, states, wellSols ] = calibrateSleipnerSetup( 'num_years',num_years );
Gt = smodel.G;

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


%% Set limits of variables to optimize (will be stored in scaling.boxLims)
% We calibrate specific variables of the benchmark model within certain
% limits. The limits for grid elevation change (dzLims) refer to the
% dimensional unit of meters for which the initial grid can change. The
% limits for CO2 density, permeability, and porosity refer to the
% multipliers applied to the initial variable values. The limits for CO2
% entry rates (qLims) refer to the (constant) multiplier applied to the
% initial 1999-2010 entry rates.
q_mean = mean(arrayfun(@(x)x.W.val, schedule.control));
dzLims  = [-5 5];           % elevations
qLims   = [0.1 1.5]*q_mean; % rates (same multiplier applied to each rate)
rhoLims = [0.4 2];          % CO2 density
permLims = [0.5 12];        % permeability
poroLims = [0.5 1.4];       % porosity


%% Define objective function and compute initial value
obj = @(wellSols, states, schedule, varargin)...
        matchToDataSens(smodel, wellSols, states, schedule, newplumes, varargin{:});

v_base = obj(wellSols, states, schedule);
v_base = cell2mat(v_base);


%% Set up box limits for scaling
scaling.boxLims = [ones(smodel.G.cells.num, 1)*dzLims; qLims; rhoLims; permLims; poroLims];
scaling.obj     = max(sum(abs(v_base)), 1);

% Convert schedule-params to control vector (assume all control-steps have same)
q_mean = mean(arrayfun(@(x)x.W.val, schedule_base.control));
c  = schedule.control(1);
u  = [c.dz; q_mean; c.rhofac; c.permfac; c.porofac];
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
u  = (u-umin)./(umax-umin);


%% Define function handle for objective evaluation
smodel.nonlinearTolerance = 1e-9;
f = @(u)evalObjectiveAndSens(u, obj, initState, smodel, schedule, scaling);


%% Calibration
[~, u_opt_s, history] = unitBoxBFGS(u, f);

% recompute optimal
[v_opt_s, der_opt, wellSols_opt, states_opt] = evalObjectiveAndSens(u_opt_s, obj, initState, smodel, schedule, scaling);
% scale back:
v_opt = -v_opt_s*scaling.obj;

[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
u_opt = u_opt_s.*(umax-umin)+umin;

fprintf('rate multiplier: %4.2f\n', u_opt(end-3)/q_mean);
fprintf('rho mutiplier: %4.2f, perm mutiplier: %4.2f, poro mutiplier: %4.2f\n', u_opt(end-2), u_opt(end-1), u_opt(end)) 

% Final simulated heights vs. observed plume thicknesses
plotObsAndSim(Gt, newplumes, states_opt)


%% Store results
% We store the convergence details of the calibration (history), the
% multipliers of the variables (mult), as well as the calibrated variables
% in their proper units (perm, poro, rho, rate, dz). The calibrated
% variables of perm and poro refer to their averaged values. The rate
% multiplier is applied to each annual rate provided in the benchmark.
history.hess = history.hess{end};   % NB: store only converged H to reduce datasize
mult.poro = u_opt(end);             % NB: full u stored in history.u{end}
mult.perm = u_opt(end-1);
mult.rho  = u_opt(end-2);
mult.rate = u_opt(end-3)/q_mean;
res = struct('multipliers', mult, ...
             'perm',    mult.perm * convertTo(mean(smodel.rock.perm),darcy), ...
             'poro',    mult.poro * mean(smodel.rock.poro), ...
             'rho',     mult.rho * smodel.fluid.rhoGS, ...
             'rate',    mult.rate .* (arrayfun(@(x)x.W.val, schedule_base.control)), ...
             'dz',      u_opt(1:end-4), ...
             'v_base',  sum(v_base), ...
             'v_opt',   v_opt, ...
             'scaling', scaling, ...
             'states_opt', {states_opt}, ...
             'history', history);

save('calibrateSleipner_results.mat', 'res', '-v7.3')


%% Some final output:
% This is one possible solution that is found within a "family of optimals".
fprintf('\nVariable             Initial    Calibrated    Unit \n')
fprintf('----------------------------------------------------\n')
fprintf(' avg. porosity      |  %1.3f   |  %1.3f      |      \n', mean(smodel.rock.poro), mult.poro*mean(smodel.rock.poro))
fprintf(' avg. permeability  |  %2.3f   |  %2.3f      | darcy\n', convertTo(mean(smodel.rock.perm),darcy), mult.perm*convertTo(mean(smodel.rock.perm),darcy))
fprintf(' CO2 density        |  %4.1f   |  %4.1f      | kg/m3\n', smodel.fluid.rhoGS, mult.rho*smodel.fluid.rhoGS)
fprintf(' avg. CO2 entry rate| %0.4f   | %0.4f      | m3/s \n', q_mean, mult.rate*q_mean)

% Initial versus calibrated grid:
figure
subplot(1,2,1)
title({'Initial';'top-surface elevations'})
plotCellData(smodel.G, smodel.G.cells.z, 'edgealpha',0.1); colorbar; axis equal tight off
subplot(1,2,2)
title({'Calibrated';'top-surface elevations'})
plotCellData(smodel.G, smodel.G.cells.z + res.dz, 'edgealpha',0.1); colorbar; axis equal tight off

