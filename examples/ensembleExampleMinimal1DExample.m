%% Ensemble of 1D reservoirs
% This is a minimal example for creating an ensemble of 1D reservoirs with
% two phases (oil and water) driven by an injector and a producer well
% located at each end of the reservoir.

mrstModule add ad-core ad-blackoil mrst-gui ad-props ...
    example-suite incomp ensemble

mrstVerbose off

%% Choose a base problem and ensemble size
% The base problem contains all properties that are common throughout the
% ensemble, and here we have implemented our own super simple example (see 
% example_template.m for the MRSTExample template).

baseProblemName = 'ensemble_base_problem_1d_reservoir';
numCells = 30;
baseProblemOptions = {'ncells', numCells};

baseExample = MRSTExample(baseProblemName, baseProblemOptions{:});

% Change these flags to investigate the baseExample
simulateExample = false;
plotExample = false;
rerunBaseProblemFromScratch = false;

if simulateExample
    problem = baseExample.getPackedSimulationProblem();
    if rerunBaseProblemFromScratch
        clearPackedSimulatorOutput(problem);
    end
    simulatePackedProblem(problem);

    [wellSols, states, reports] = getPackedSimulatorOutput(problem);
    if plotExample
        baseExample.plot(states);
    end
end


%% Create samples that represent the stochastic parameters in our ensemble

ensembleSize = 20;

configData = cell(ensembleSize, 1);
for i = 1:ensembleSize
    configData{i}.poro = gaussianField(baseExample.model.G.cartDims, [0.2 0.4]); 
    configData{i}.perm = configData{i}.poro.^3.*(1e-5)^2./(0.81*72*(1-configData{i}.poro).^2);
end

samples = RockSamples('data', configData);
disp(samples);

%% Select quantity of interest class

qoi = WellQoI('wellNames', {'P1'}, 'cumulative', true, ...
    'fldname', {'qOs', 'qWs'});
disp(qoi);

%% Create the ensemble

ensemble = MRSTEnsemble(baseExample, samples, qoi, ... 
    'simulationStrategy', 'parallel', ...
    'maxWorkers', 8, ...
    'verbose', true, ...
    'reset', true ...
    );
disp(ensemble);

%% Run ensemble
ensemble.simulateEnsembleMembers('plotProgress', true);

%% Plot results
close all, ensemble.plotQoI();

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
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
