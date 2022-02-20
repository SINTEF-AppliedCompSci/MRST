%% Lightweight diagnostics for egg ensemble
%
% This example shows how to setup an ensemble, calculate flow diagnostics based
% on residence time distributions and visualise results with the lightweight
% EnsembleGUI. Here we use the Egg ensemble.
%
% For more information see Chapter 13 of the MRST book "An Introduction to
% Reservoir Simuation in MATLAB / GNU Octave" Lie 2019.
%

%% First, we create an instance of the ModelEnsemble class for the Egg ensemble.
% This requires a name, number of members and a corresponding setup function.
% Once setup, diagnostics and full simulations can be run directly using
% m.computeDiagnostics() and m.runSimulations(). Results are easy to access
% using mrst ResultHandlers.
%

mrstModule add diagnostics ad-core ad-props ad-blackoil incomp mrst-gui example-suite

m = ModelEnsemble('egg', 'nMembers', 50, 'setupFn', @setupEggFn);
% after initial setup, m = ModelEnsemble('egg') is sufficient

diagnNo = 1;
if numel(m.diagnostics) > 1
    fprintf('Ensemble contains multiple diagnostics scenarios, selecting first\n');
end

%% Compute diagnostics
m.computeDiagnostics();


%% Finally we pass the ModelEnsemble object to the EnsembleGUI for interactive
% plotting and inspection of results.
g = EnsembleGUI(m);

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

