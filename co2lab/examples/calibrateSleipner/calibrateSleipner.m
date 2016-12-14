%% Calibrate Sleipner model 
% This example calibrates the Sleipner model (top surface elevations, rock
% property, CO2 entry rates, CO2 density) by minimizing the misfit between
% simulated and observed CO2 plume thicknesses.

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
q_mean = mean(arrayfun(@(x)x.W.val, schedule.control));
dzLims  = [-5 5];
qLims   = [0.1 1.5]*q_mean;
rhoLims = [0.1 2];
permLims = [0.5 12];
poroLims = [0.5 1.4];


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
% @@
plotObsAndSim(Gt, newplumes, states_opt)


%% Store results
history.hess = history.hess{end};       % NB: store only converged H to reduce datasize
res = struct('perm', u_opt(end-1), ...  % NB: full u stored in history.u{end}
             'poro', u_opt(end), ...
             'rho', u_opt(end-2), ...
             'rate', u_opt(end-3)/q_mean, ...
             'dz', u_opt(1:end-4), ...
             'v_base', sum(v_base), ...
             'v_opt', v_opt, ...
             'scaling', scaling, ...
             'states_opt', {states_opt}, ...
             'f', f, ...
             'history', history);

save(['calibrateSleipnerResults.mat'], 'res', '-v7.3')