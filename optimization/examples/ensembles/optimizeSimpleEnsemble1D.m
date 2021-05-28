%% Run an NPV-optimization problem over an ensemble of simple 1D-models
mrstModule add ad-core ad-blackoil ad-props mrst-gui optimization


%% setup optmization-problem
% Select objective function handle. If simulations/adjoints are to be run in 
% seperate matlab-sessions, anonymous functions should be avoided. Extra objective 
% arguments can be added through 'objectiveOpts'.
objective      = @NPVOW;
objectiveOpts  = {'OilPrice',           40/stb, ...
                  'WaterInjectionCost'   5/stb, ...
                  'WaterProductionCost', 4/stb};
              
% Provide handle to function that sets up problem-structure as a function
% of realization number.
setupFn = @setupFnSimpleProblem1D;

% Provide initial guess in terms of a schedule. Only control-targets and
% limits are used, so this can be the schedule for any of the simulation 
% problems
problem_tmp = setupFn(1);
initialGuess = problem_tmp.SimulatorSetup.schedule;

% Provide indices to the realizations that should be included in the
% optimization. Here we use three
realizations = (1:3);

% Setup target and limit bounds (initial schedule should lie within bounds). 
% All upper/lower limits that are finite in shedule, and defined below, will
% become active controls if they become active during at least one of the
% simulations.
W_tmp = initialGuess.control(1).W;
irate = W_tmp(1).val;
bnds  = processBounds(W_tmp, 'rate(inj)',  [.1*irate, 3*irate], ... % allowed target rates
                             'bhp(inj)',   [120, 500]*barsa, ...    % allowed upper bhp-limits
                             'bhp(prod)',  [50, 120]*barsa, ...     % allowed target bhp
                             'lrat(prod)', [-3*irate, -.1*irate]);  % allowed lower lrat-limit                                              
                                                
% It can be advantageous to provide an objective value guess (at least order 
% of magnitude), to help the optimizer in setting an initial step-length. 
% This value is set as
G = problem_tmp.SimulatorSetup.model.G;
objectiveScaling = .9*40*sum(G.cells.volumes)/stb; 

% Function-evaluations can be run in seperate matlab-sessions. Depending on
% problem size, this number should not be gratear than maxNumCompThreads().
% Here, we run evaluations in the same session as the optimizer:
nWorkers = 0;

% Finally, set up the optimization-problem
p = OptimizationProblem('SimpleEnsembleOptimization1D', ...      % problem name                                
                        'objective',        objective, ...                
                        'objectiveOpts',    objectiveOpts, ...                    
                        'setupFn',          setupFn, ...
                        'initialGuess',     initialGuess, ...
                        'realizations',     realizations, ...
                        'bounds',           bnds, ...
                        'objectiveScaling', objectiveScaling, ...
                        'nWorkers',         nWorkers)   %#ok
                    
%% Run optimization 
runBackground = false;
if ~runBackground
    % running in current session will trigger plotting from the optimizer
    % if it is set to default (unitBoxBFGS)
    scheduleOpt = p.optimize('stepInit', .1, 'maxIt', 6);
else
    p.optimizeBackground('stepInit', .1, 'maxIt', 6);
end

%% Once the optimization is done, we can inspect the optimization results
p = p.refreshControlList();
% plot all the computed NPV-values
figure, p.plotObjectiveValues();
% get objective value for each iteration 
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
