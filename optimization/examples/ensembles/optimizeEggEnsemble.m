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

mrstModule add ad-core ad-blackoil ad-props mrst-gui optimization test-suite ensemble
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

realizations = [10 20 30 40];
samples = struct(...
    'getSampleProblem', @(dd, seed)setupFnEggEnsemble(seed, dd, realizations), ...
    'num', numel(realizations));

% Provide initial guess in terms of a schedule. Only control-targets and
% limits are used, so this can be the schedule for any of the simulation 
% problems
problem_tmp = setupFn(realizations(1));
schedule = problem_tmp.SimulatorSetup.schedule;


% Setup target/limit upper and lower bounds (initial schedule should lie 
% within bounds). All upper/lower limits that are finite in schedule, and 
% defined below, will contribute to gradient if they become active during at least 
% one of the simulations. Bounds can be set up individually for any 
% target/limit of any well. The function processBounds simplifies the setup, 
% when bounds are equal across all injectors/producers.
W_tmp = schedule.control(1).W;
bnds = processBounds(W_tmp, 'rate(inj)', [1 150]/day, ...       % target rate bounds
                            'bhp(inj)', [400 420]*barsa, ...    % upper bhp limit bounds
                            'bhp(prod)', [380 400]*barsa);      % target bhp bounds                            
                                                
maps = setupSimulationControlMappings(schedule, bnds);                        
% It can be advantageous to provide an objective value guess (at least order 
% of magnitude), to help the optimizer in setting an initial step-length. 
% This value should be set as
objectiveScaling = 5e8;

objFun = @(model, states, schedule, varargin)NPVOW(model, states, schedule, varargin{:}, objectiveOpts{:});
objStruct = struct('function', objFun, ...
                   'scaling', objectiveScaling);

% Function-evaluations can be run in seperate matlab-sessions. 
if maxNumCompThreads() > 2
    simStrategy = 'background';
else
    simStrategy = 'serial';
end

% Finally, set up the optimization-problem
p = OptimizationProblem(samples, ...
                        'name', 'EggEnsembleOptimization', ...                                   
                        'objective',        objStruct, ...                
                        'maps',             maps,...
                        'setupType',     'simulation', ...
                        'verboseSimulation',     true, ...
                        'simulationStrategy',  simStrategy)   %#ok
                          
%% Run optimization 
% note that this requires aprox 10 times number of realizations simulations
% and adjoints

%p.reset('prompt', false); % un-comment to rerun optimization

% optimization can be run in seperate session:
runBackground = true;
u = p.maximizeObjective(problem_tmp, 'maxIt', 6, 'background', runBackground);


%% While optimization is running, progress can be checked with 
figure, p.plotObjectiveValues();

 
%% Once the optimization is done, we can inspect the optimization results

ids = p.iterationObjectiveValues.getValidIds;
if ~isempty(ids)
    % initial schedule well solutions
    p.plotWellSols(problem_tmp, 1);
    if max(ids) > 1
        % optimized schedule well solutions
        p.plotWellSols(problem_tmp, max(ids));
    end
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
