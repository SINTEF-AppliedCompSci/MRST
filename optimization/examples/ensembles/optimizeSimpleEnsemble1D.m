%% Run an NPV-optimization problem over an ensemble of simple 1D-models
mrstModule add ad-core ad-blackoil ad-props mrst-gui optimization ensemble


%% setup optmization-problem
% Select objective function handle. If simulations/adjoints are to be run in 
% seperate matlab-sessions, anonymous functions should be avoided. Extra objective 
% arguments can be added through 'objectiveOpts'.
objective      = @NPVOW;
objectiveOpts  = {'OilPrice',           40/stb, ...
                  'WaterInjectionCost'   5/stb, ...
                  'WaterProductionCost', 4/stb};
              
% Setup base problems
numRealizations = 2;
problems = cell(1, numRealizations);
for k =1:numRealizations
    problems{k} = setupFnSimpleProblem1D(k);
    problems{k}.seed = k;
end
samples = struct('problem', {problems}, ...
                  'num', numRealizations);

% Provide indices to the realizations that should be included in the
% optimization. Here we use three
schedule = problems{1}.SimulatorSetup.schedule;

%%
% Setup target and limit bounds (initial schedule should lie within bounds). 
% All upper/lower limits that are finite in shedule, and defined below, will
% become active controls if they become active during at least one of the
% simulations.
W = schedule.control(1).W;
irate = W(1).val;
bnds  = processBounds(W, 'rate(inj)',  [.1*irate, 3*irate], ... % allowed target rates
                             'bhp(inj)',   [120, 500]*barsa, ...    % allowed upper bhp-limits
                             'bhp(prod)',  [50, 120]*barsa, ...     % allowed target bhp
                             'lrat(prod)', [-3*irate, -.1*irate]);  % allowed lower lrat-limit                                              
[maps, u] = setupSimulationControlMappings(schedule, bnds);


% It can be advantageous to provide an objective value guess (at least order 
% of magnitude), to help the optimizer in setting an initial step-length. 
% This value is set as
G = problems{1}.SimulatorSetup.model.G;
objectiveScaling = .9*40*sum(G.cells.volumes)/stb; 

objFun = @(model, states, schedule, varargin)NPVOW(model, states, schedule, varargin{:}, objectiveOpts{:});
objStruct = struct('function', objFun, ...
                   'scaling', objectiveScaling);
% Function-evaluations can be run in seperate matlab-sessions. Depending on
% problem size, this number should not be gratear than maxNumCompThreads().
% Here, we run evaluations in the same session as the optimizer:
nWorkers = 0;

% Finally, set up the optimization-problem
p = OptimizationProblem(samples, ...
                        'name', 'SimpleEnsembleOptimization1D', ...      % problem name                                
                        'objective',        objStruct, ...                
                        'maps',             maps,...
                        'setupType',     'simulation', ...
                        'verboseSimulation',     true)   %#ok
                    
%% Run optimization 
p.reset('prompt', false)
u = p.optimize(problems{1}, 'stepInit', .1, 'maxIt', 6);


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
