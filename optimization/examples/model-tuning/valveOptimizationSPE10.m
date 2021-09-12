%% Parameter tuning of a very coarse upscaling of the Egg model   
mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid...
               mrst-gui ad-props incomp optimization...
               example-suite linearsolvers spe10

%% Setup reference model
example  = MRSTExample('spe10_wo', 'layers', 1);
problem  = example.getPackedSimulationProblem();
model = problem.SimulatorSetup.model;
% Consider 1 injector (south) and one producer (north)
W = addWell([], model.G, model.rock, 1:60:220*60, 'Name', 'I', 'Type', 'rate', 'Val', 30/day, ...
            'Sign', 1, 'compi', [1 0]);
W = addWell(W, model.G, model.rock, 60:60:220*60, 'Name', 'P', 'Type', 'bhp', 'Val', 150*barsa, ...
            'Sign', -1', 'compi', [0 1]);
schedule = problem.SimulatorSetup.schedule;
schedule.control.W = W;
schedule.step.control = schedule.step.control(1:36);
schedule.step.val = schedule.step.val(1:36);
problem.SimulatorSetup.schedule = schedule;
problem.SimulatorSetup.NonLinearSolver.useLinesearch = true;
%problem.SimulatorSetup.NonLinearSolver = NonLinearSolver();
problem.SimulatorSetup.model.fluid = ...
    initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                       'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                       'n',     [2, 2, 0]);

problem.SimulatorSetup.state0.s = ones(model.G.cells.num,1)*[0 1];                   

%problem.SimulatorSetup.model.verbose = true
simulatePackedProblem(problem);

[wellSolsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef    = problem.SimulatorSetup.model;



% Plot
%example.plot(statesRef, 'step_index', numel(statesRef))
%plotWellSols(wellSolsRef)

%% Specify WI-multipliers for injector and producer imitate valves
steps = schedule.step.val;
% valves will be allowed to change every 60 days, so we refomulate the
% shedule to allow for this
controlstep = ceil(cumsum(steps)/(60*day));
schedule.step.control = controlstep;
schedule.control = repmat(schedule.control(1), [max(controlstep), 1]);

setup = problem.SimulatorSetup;
setup.schedule = schedule;
% Wach of the two wells has 220 connections, we set one valve for evere 22 connections,
% such that each well has 10 valves. This is set up by the lumping-vector
lumping = ceil((1:440)'/20);
% Add parameters for each control step (period where valves are constant)
params = [];
for k = 1:max(controlstep)
    params = addParameter(params, setup, 'name', 'conntrans', 'type', 'multiplier', ...
                          'lumping', lumping, 'controlSteps', k, ...
                          'boxLims', [0.0001 1], 'scaling', 'log');
end

%%

%% Define the mismatch function
% The mismatch function is defined as a function handle to a library
% function from the optimization module that computes the mismatch between
% a given simulation and a reference state. For an oil-water system, the
% match is computed based on three quantities (water/oil rate and bhp) and
% these should be given an associated weight. Weights on the order the
% reciprocal of rate magnitudes/pressure variations result in a properly 
% scaled mismatch function 
options  = {'OilPrice',             50 , ...
            'WaterProductionCost',   5 , ...
            'WaterInjectionCost',    3 , ...
            'DiscountFactor',       0.1};
        
npvFn = @(model, states, schedule, tmp, cp, tstep, state)...
    NPVOW(model, states, schedule, 'computePartials', cp, 'tstep', tstep, ...
          'state', state, options{:}, 'from_states', false);
% mismatchFn = @(model, states, schedule, states_ref, tt, tstep, state) ...
%     matchObservedOW(model, states, schedule, states_ref,...
%                    'computePartials', tt, 'tstep', tstep, weighting{:},...
%                    'state', state, 'from_states', false);

val = NPVOW(model, statesRef, schedule, options{:});
val = sum(vertcat(val{:}));

%objhh = @(tstep,model,state) npvFn(setup.model, statesRef, setup.schedule, [], true, tstep, state);
%gradient = computeSensitivitiesAdjointAD(setup, statesRef, params, objhh);        




%% Model calibration
pvec = getScaledParameterVector(setup, params);
objh = @(p) evaluateMatch(p, npvFn, setup, params, [], 'objScaling', -val);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v, p_opt, history] = unitBoxBFGS(pvec/100, objh, 'objChangeTol', 1e-5, ...
    'gradTol', 1e-4, 'maxIt', 15, 'lbfgsStrategy', 'dynamic', 'lbfgsNum', 5);

save opt_res p_opt history

%% Create a new coarse model setup with the optimal parameters, 
%% and rerun the simulation for full time horizon
setup_opt = updateSetupFromScaledParameters(setup, params, p_opt); 
% reset time-steps to full schedule
% setup_opt.schedule.step = schedule.step;
[wellSols_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, setup_opt.schedule);
return

%%
% arrange mulitpliers for imspection
nc = numel(params);
[injMult, prodMult] = deal(nan(nc, 10));
for k = 1:nc
    tmp = params{k}.getParameter(setup_opt);
    injMult(k,:) = tmp(1:10);
    prodMult(k,:) = tmp(11:20);
end
gp = cartGrid([size(injMult)]);
for k = 1:2
    figure
    if k == 1
        plotCellData(gp, log(injMult(:)), 'EdgeAlpha', .2);
        title('Injector valves');
    else
        plotCellData(gp, log(prodMult(:)), 'EdgeAlpha', .2);
        title('Producer valves');
    end
    ax = gca;
    ax.XTick = 0:2:14;
    ax.XTickLabel= (0:2*60:14*60)
    ax.YTick = .5:10;
    ax.YTickLabel = 1:10;
    xlabel('Time [days]')
    ylabel('Valve No')
    ax.FontSize = 14;
end

figure
plotCellData(model.G, log(model.rock.perm(:,1)), 'EdgeAlpha', .1);
axis off
figure
plotCellData(model.G, statesRef{end}.s(:,2), 'EdgeAlpha', .1);
axis off
figure
plotCellData(model.G, states_opt{end}.s(:,2), 'EdgeAlpha', .1);
axis off

