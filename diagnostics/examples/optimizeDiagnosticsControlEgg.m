%% Control optimization: mean of diagnostics proxy for NPV over ensemble of models
% This example computes well controls (rate/BHP) that optimize an efficient
% NPV-proxy over a set of realizations of the Egg model, and optionally runs
% simulations for initial and optimized controls (more time consuming)

% Run simulations for initial and optimized controls (see end of script):
runSimulations = false;

mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat optimization
mrstModule add linearsolvers agmg sequential incomp diagnostics ensemble test-suite

%% define optimization problem
% select ensemble members and create problems
realizations = [1, 2, 3];
for k = 1:numel(realizations)
    samples.problem{k} = setupEggForDiagnosticsOptimization(realizations(k));
    samples.problem{k}.seed    = k;
end
samples.num = numel(realizations);

% select objective options
objOpts = {'ro' ,        50/stb, ...
    'rwi',        -5/stb, ...
    'rwp',        -3/stb, ...
    'discount',    0.1,   ...
    'timeHorizon', 2.5};  % default time-unit for discount/horizon is year
obj = struct('function',    @(model, varargin)DiagnosticsNPV(model, objOpts{:}, varargin{:}), ...
             'scaling', 7.5e7);

% Impose bounds on all wells
scheduleInit = samples.problem{1}.SimulatorSetup.schedule;
bnds  = processBounds(scheduleInit(1).control.W, ...
    'rate(inj)',   [1 150]/day, ...      % target rate bounds
    'bhp(inj)',  [400 420]*barsa, ...    % upper bhp limit bounds
    'bhp(prod)', [380 400]*barsa);       % target bhp bounds
% create mappings from/to control vector to/from wells
maps = setupSimulationControlMappings(scheduleInit, bnds);

% setup optimization problem
p = OptimizationProblem(samples, ...
    'name', 'DiagnosticsControlOptimizationEgg', ...
    'objective',     obj, ...
    'maps',          maps,...
    'setupType',     'diagnostics');

%% run optimization
% fetch initial scaled control vector
uOpt = p.optimize(samples.problem{1}, 'maxIt', 15, 'objChangeTol', 1e-4);
% plot all objective values
figure, p.plotObjectiveValues();
% uncomment line below to delete optimization output:
% p.reset('prompt', false)

%% Simulations ------------------------------------------------------------
% run full simulations for base and optimal cases (if opted)
if runSimulations
    problems = samples.problem;
    Wopt    = problem.SimulatorSetup.schedule.control.W;
    [statesInit, statesOpt] = deal(cell(1, numel(realizations)));
    
    for k = 1:numel(realizations)
        state0 = problems{k}.SimulatorSetup.state0;
        model  = problems{k}.SimulatorSetup.model.getReservoirModel();
        nls = NonLinearSolver;
        nls.LinearSolver = selectLinearSolverAD(model);
        
        % setup 3 year schedule with initial controls
        W = problems{k}.SimulatorSetup.schedule.control.W;
        scheduleInit = simpleSchedule(rampupTimesteps(3*year, 30*day), 'W', W);
        % simulate
        [~, statesInit{k}] = simulateScheduleAD(state0, model, scheduleInit, 'NonLinearSolver', nls);
        
        % update problem to optimal controls
        problems{k} = p.updateProblemFun(problems{k}, uOpt);
        % setup 3 year schedule with optimized controls
        W = problems{k}.SimulatorSetup.schedule.control.W;
        scheduleOpt = simpleSchedule(rampupTimesteps(3*year, 30*day), 'W', W);
        % simulate
        [~, statesOpt{k}] = simulateScheduleAD(state0, model, scheduleOpt, 'NonLinearSolver', nls);
    end
    
    % Compute evolution of NPV for all realizations
    [npvInit, npvOpt] = deal(cell(1, numel(realizations)));
    npvOpts = {'OilPrice',               50,  ...
        'WaterProductionCost',     3, ...
        'WaterInjectionCost',      5, ...
        'DiscountFactor',         .1};
    
    for k = 1:numel(realizations)
        tmp = NPVOW(model, statesInit{k}, scheduleInit, npvOpts{:});
        npvInit{k} = cumsum(vertcat(tmp{:}));
        tmp = NPVOW(model, statesOpt{k}, scheduleOpt, npvOpts{:});
        npvOpt{k} = cumsum(vertcat(tmp{:}));
    end
    
    % plot NPV-evolutions
    t = cumsum(scheduleInit.step.val)/year;
    figure, hold on
    plot(t, mean(horzcat(npvInit{:}),2), '-k', 'LineWidth', 2);
    plot(t, mean(horzcat(npvOpt{:}),2), '--k', 'LineWidth', 2);
    if numel(realizations)==1
        leg = {'Initial', 'Optimized'};
    else
        ax = gca; ax.ColorOrderIndex = 1;
        plot(t, horzcat(npvInit{:}), '-');
        ax = gca; ax.ColorOrderIndex = 1;
        plot(t, horzcat(npvOpt{:}), '--');
        leg = [{'Initial mean', 'Optimized mean'}, ...
            arrayfun(@(r)sprintf('Initial %d', r), realizations, 'UniformOutput', false), ...
            arrayfun(@(r)sprintf('Optimized %d', r), realizations, 'UniformOutput', false)];
    end
    legend(leg, 'Location', 'southeast')
    axis([2 3 7.3e7 7.8e7])
    xlabel('Time [year]'), ylabel('NPV [$]')
    
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
