%% Run an expected NPV-optimization over a subset of the Egg-ensemble 
% NOTE: this example requires ~40 simulations and adjoint runs of
% realizations from the Egg ensamble. The optimization can be cancelled 
% (or killed) at any stage, and be continued at a later time by rerunning 
% the script. Required simulation-, gradient- and objective values are saved 
% in mrstOutputDirectory(). 
%
% The example optimizes the mean (expected) NPV for the selected
% realizations. The optimization procedure is based on a version of the
% generalized reduced gradient (GRG) method, which computes gradients both
% for well target controls and well upper/lower simulation limits. Hence,
% both targets and limits (that become active during at least one
% simulation) are adjusted to optimize the objective.

mrstModule add ad-core ad-blackoil ad-props mrst-gui optimization example-suite
mrstVerbose false

%% setup optmization-problem
% Select objective function handle. If simulations/adjoints are to be run in 
% seperate matlab-sessions, anonymous functions should be avoided. Extra objective 
% arguments can be added through 'objectiveOpts'.
objective      = @NPVOW;
objectiveOpts  = {'OilPrice',           50/stb, ...
                  'WaterInjectionCost'   3/stb, ...
                  'WaterProductionCost', 3/stb};
              
% Provide handle to function that sets up problem-structure as a function
% of realization number. The function setupFnEggEnsemble also sets up a
% control-logic for each simulation such that wells are shut whenever
%   * water-cut exceeds 0.95
%   * producers become injectors and vice versa
%   * well-rates drop below 0.1 m^3/day
setupFn = @setupFnEggEnsemble;

% Provide initial guess in terms of a schedule. Only control-targets and
% limits are used, so this can be the schedule for any of the simulation 
% problems
problem_tmp = setupFn(1);
initialGuess = problem_tmp.SimulatorSetup.schedule;

% Provide indices to the realizations that should be included in the
% optimization. Here we use four
realizations = [10 20 30 40];

% Setup target/limit upper and lower bounds (initial schedule should lie 
% within bounds). All upper/lower limits that are finite in schedule, and 
% defined below, will become active if they become active during at least 
% one of the simulations. Bounds can be set up individually for any 
% target/limit of any well. The function processBounds simplifies the setup, 
% when bounds are equal across all injectors/producers.
W_tmp = initialGuess.control(1).W;
bnds = processBounds(W_tmp, 'rate(inj)', [1 150]/day, ...       % target rate bounds
                            'bhp(inj)', [400 420]*barsa, ...    % upper bhp limit bounds
                            'bhp(prod)', [380 400]*barsa);      % target bhp bounds                            
                                                
% It can be advantageous to provide an objective value guess (at least order 
% of magnitude), to help the optimizer in setting an initial step-length. 
% This value should be set as
objectiveScaling = 1e8;

% Function-evaluations can be run in seperate matlab-sessions. Depending on
% problem size, this number should not be gratear than maxNumCompThreads().
% Here, we run evaluations in the same session as the optimizer:
nWorkers = maxNumCompThreads() -1;
if nWorkers < 2
    % don't bother for just a single thread 
    nWorkers = 0;
end

% Finally, set up the optimization-problem
p = OptimizationProblem('EggEnsembleOptimization', ...      % problem name                                
                        'objective',        objective, ...                
                        'objectiveOpts',    objectiveOpts, ...                    
                        'setupFn',          setupFn, ...
                        'initialGuess',     initialGuess, ...
                        'realizations',     realizations, ...
                        'bounds',           bnds, ...
                        'objectiveScaling', objectiveScaling, ...
                        'nWorkers',         nWorkers)       %#ok
                    
%% Run optimization 
% note that this requires aprox 10 times number of realizations simulations
% and adjoints
runBackground = true;
if ~runBackground
    scheduleOpt = p.optimize('stepInit', .5, 'maxIt', 7);
else
    p.optimizeBackground('stepInit', .5, 'maxIt', 7);
end

%% While optimization is running, progress can be checked with 
figure, p.plotObjectiveValues();

%% or for continous plotting
figure, p.monitorProgress();
 
%% Once the optimization is done, we can inspect the optimization results
% get objective value for each iteration 
p = p.refreshControlList();
values = p.loadAllFiles('NPVOW.mat', 'value');
isfin   = cellfun(@isfinite, values);
if any(isfin)
    values = cell2mat(values(isfin));
    % plot well-solutions for initial schedule
    p.plotWellSols(1);
    % plot well-solution for best performing schedule
    [~, ix] = max(values);
    if ix ~= 1
        p.plotWellSols(ix);
    end
end

% to delete all output: p.cleanup()

%% Copyright Notice
%
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
