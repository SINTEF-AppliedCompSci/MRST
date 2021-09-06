%% Diagnostics GUI for plotting 3D flow diagnostics fields for multiple ensemble
% members
%
% This example shows how to load an ensemble of models into the interactive
% diagnostics GUI. Each ensemble member must have the same grid and well
% configuration. The diagnosticsViewer class will calculate simple flow
% diagnostics for each ensemble member and allow you to compare both grid based
% and other properties for each model. These flow diagnostics are based on
% results from a single phase incompressible simulation.
%
% For more information see Chapter 13 of the MRST book "An Introduction to
% Reservoir Simuation in MATLAB / GNU Octave" Lie 2019.
%
% Here we run flow diagnostics on multiple realizations of the Egg model and
% load the results into the interactive GUI. We setup an intial state with an
% an oil cap above the water zone. This setup is relatively arbitrary but
% demonstrates the use of basic diagnostic multi-phase recovery estimates.
%

%% Load modules

mrstModule add diagnostics ad-core ad-blackoil example-suite
%% Get the first 10 realisations of the Egg model

% Here we setup the ensemble using the ModelEnsemble class. 
% Each ensemble member has a different permeability but all other properties
% are the same.

m = ModelEnsemble('egg', 'nMembers', 50, 'setupFn', @setupEggFn);
% after initial setup, m = ModelEnsemble('egg') is sufficient

% Input to the diagnostics viewer is generated using the ensemble class
% setup function.
realizations = 1:10;

[models, wells] = deal(cell(1, numel(realizations)));
for k = 1: numel(realizations)
    tmp = m.setupFn(k);
    [models{k}, wells{k}] = deal(tmp.model, tmp.W);
end


%% Setup the two-phase intial condition
% This has no bearing on the calculation of basic time-of-flight and tracer
% diagnostics but allows us to make rough estimates of recovery and water cut
% for specific models / wells.

gcz = models{1}.G.cells.centroids(:, 3);
height = max(gcz) - min(gcz);
pos = 1 - (gcz - min(gcz))./height;

oil = pos > 0.1;
wat = pos < 0.1;

[state0] = deal(cell(1, numel(realizations)));

for i = 1:numel(realizations)
    
    state0{i} = initResSol(models{1}.G, 200*barsa);
    state0{i}.s(wat, 1) = 1;
    state0{i}.s(oil, 2) = 1;
    state0{i}.s = bsxfun(@rdivide, state0{i}.s, sum(state0{i}.s, 2));

end

%% Launch the diagnostics GUI.
% Assigning it to d allows us to inspect and access properties of the GUI after
% it has been launched. 

[d] = DiagnosticsViewer(models,wells,'state0',state0);


%%
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


