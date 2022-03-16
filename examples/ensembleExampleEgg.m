%% Egg ensemble example
% In this example, we show how to set up an ensemble from a collection of
% deck files, exemplified by the Egg model [1]

%% Add modules
mrstModule add ad-core ad-props ad-blackoil
mrstModule add test-suite ensemble
mrstModule add mrst-gui
mrstVerbose on

%% Set up base problem
baseCase = TestCase('egg_wo');
% Extract interesting part of the schedule
steps = cumsum(baseCase.schedule.step.val) <= 1500*day;
baseCase.schedule.step.val     = baseCase.schedule.step.val(steps);
baseCase.schedule.step.control = baseCase.schedule.step.val(steps);
% Plot setup
baseCase.plot(baseCase.model.rock, 'log10', true); colormap(pink);
material dull

%% Set up samples
% We use the function getDeckEGG as a generatorFn to set up a sample. The
% realizations are numbered from 0 to 100, hence seed-1.
generatorFn = @(problem, seed) getDeckEGG('realization', seed-1);
% The class DeckSamples implements what we need for setting up samples from
% a deck file. The class uses `initEclipseProblemAD` to set up the problem,
% and extra input arguments to this function can be provided in the class
% property `initArgs`. The property `gridFromDeck` is used to determine if
% the simulation grid should be constructed from the deck, or taken from
% the base problem. We are only interested in simulating the first 1500
% days. Since `initEclipseProblemAD` will set up the schedule from deck, we
% must postprocess the problem after setting the sample in order to extract
% the subschedule. We do this with the optional input argument
% `processProblemFn`.
processProblemFn = @(problem) getSubSchedule(problem, steps);
samples = DeckSamples('generatorFn'     , generatorFn     , ... % Generator function
                      'processProblemFn', processProblemFn, ... % Get subschedule
                      'num'             , 101             );    % Number of samples
disp(samples)

%% Set up QoI
% For our QoI, we choose the total oil production rate
is_prod = vertcat(baseCase.schedule.control(1).W.sign) < 0;
qoi = WellQoI('wellIndices', is_prod, 'names', 'qOs');
    
%% Set up ensemble
ensemble = MRSTEnsemble(baseCase, samples, qoi, ...
               'simulationStrategy', 'background'); % Run in the background

%% Simulate the ensemble members
% We simulate 8 samples. Each time the code in this block is called, 8 new
% samples will be simulated until all ensemble memebers have been run. To
% reset the ensemble data, call ensemble.reset();
ensemble.simulateEnsembleMembers('batchSize', 8, 'plotProgress', true);

%% Plot the QoI
color = lines(7); color = color(end,:);
close all, ensemble.plotQoI('subplots', true, 'subplotDir', 'vertical', 'color', color);
f = gcf; f.Position(4) = f.Position(4)*2;

%%
mrstModule add sequential diagnostics
solver = @(problem) pressureSolverAD(problem);

ensembleFD = MRSTEnsemble(baseCase, samples, qoi, ...
               'solve'             , solver  , ...
               'directory'         , [ensemble.directory, '_fd'], ...
               'simulationStrategy', 'background'); % Run in the background

%%
%ensembleFD.qoi.diagnosticsType = 'tof';
ensembleFD.simulateEnsembleMembers('batchSize', 8, 'plotProgress', true);
         
%% References
% [1] Jansen, J. D., et al., "The egg model–a geological ensemble for
% reservoir simulation." 'Geoscience Data Journal 1.2 (2014): 192-195.'

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
