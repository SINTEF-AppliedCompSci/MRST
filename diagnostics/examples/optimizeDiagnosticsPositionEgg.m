%% Trajectory optimization: mean of diagnostics proxy for NPV over ensemble of models
% This example finds well trajectories that optimize an efficient
% NPV-proxy over a set of realizations of the Egg model, and optionally runs
% simulations for initial and optimized controls (more time consuming)

mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat optimization test-suite
mrstModule add linearsolvers agmg sequential incomp diagnostics wellpaths ensemble test-suite
%% define optimization problem
% select which ensemble members to optimize over (here we select a single realization)
realizations = [1];
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
obj          = struct('function',    @(model, varargin)DiagnosticsNPV(model, objOpts{:}, varargin{:}), ...
                      'scaling', 7.5e7);
                  
% select wells to optimize, here we select producer 1 and 2 (wells 9 and 10)
W = samples.problem{1}.SimulatorSetup.schedule.control(1).W; 
wno = [9, 10];
% setup instances of WellPositionControl-class 
for k = 1:numel(wno)
    W(wno(k)).posControl= WellPositionControl(samples.problem{1}.SimulatorSetup.model.G, 'w', W(wno(k)), ...
                                'perturbationSize', [16 16 8], ... % perturbation used for residual perturbations 
                                'maxUpdateWell',    [16 16 4], ... % max well update per iteration
                                'maxUpdatePoint',   [16 16 4], ... % max trajectory point update per iteration
                                'nPoints',                 2);     % number of trajectory coordinates (2 => straight line)
end
                                 
% setup optimization problem      
p = OptimizationProblem(samples, ...
                        'name', 'DiagnosticsPositionOptimizationEgg', ...
                        'objective',     obj, ...
                        'setupType',     'diagnostics', ...
                        'W',             W);   
                    
%% run optimization
uOpt = optimPlaceSimple(p, samples.problem{1}, 'stepInit', .02, 'plotTrajectories', true);
% plot all objective values
figure, p.plotObjectiveValues();
% uncomment line below to delete optimization output:
% p.reset('prompt', false)

%
%% plot initial and optimal well-positions 
problem = samples.problem{1};
G = problem.SimulatorSetup.model.G;
WInit = problem.SimulatorSetup.schedule.control.W;
WInit = addTrajectories(WInit, G, 2);
% update problem with optimized positions
problem = p.updateProblemFun(problem, uOpt);
WOpt = problem.SimulatorSetup.schedule.control.W;
for k = 1:numel(WOpt)
    WOpt(k).name = [WOpt(k).name, '(opt)'];
end
WOpt = addTrajectories(WOpt, G, 2);

figure, plotGrid(G, 'FaceColor', 'none', 'EdgeColor', [.6 .6 .6], 'EdgeAlpha', .4)
h1 = WellPlotHandle(G, WInit(setdiff((1:12), wno)));
h1.injectorColor = [.5 .5 .5]; h1.producerColor = [.5 .5 .5];
h2 = WellPlotHandle(G, WInit(wno));
h3 = WellPlotHandle(G, WOpt(wno));
set([h3.producers.body], 'LineStyle', ':')

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
