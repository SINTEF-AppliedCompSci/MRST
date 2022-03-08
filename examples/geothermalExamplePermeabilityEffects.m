%% Geothermal example illustrating p/T-dependent permeability
mrstModule add ad-core ad-props
mrstModule add geothermal compositional
mrstModule add spe10 upr
mrstModule add test-suite
mrstModule add mrst-gui

mrstVerbose on

%%
test = TestCase('permeability_effects_geothermal');

%%
problem = test.getPackedSimulationProblem('Name', 'static');
simulatePackedProblem(problem, 'restartStep', 1);

%%
K0 = 273.15*Kelvin;
Tmin = K0 + 100*Kelvin;
Tmax = K0 + 200*Kelvin;
permMult = 1e-6;
tau = @(T) (min(max(T, Tmin), Tmax) - Tmin)./(Tmax - Tmin);
perm = @(p,T) test.model.rock.perm.*((1-tau(T)) + tau(T)*permMult);

testPerm = test;
testPerm.model.rock.perm = perm;
testPerm.model = testPerm.model.setupOperators();

%%
problemPerm = testPerm.getPackedSimulationProblem('Name', 'dynamic');
simulatePackedProblem(problemPerm, 'restartStep', 1);

%%
[~, states, reports] = getPackedSimulatorOutput(problem);
[~, statesPerm, reportsPerm] = getPackedSimulatorOutput(problemPerm);
model = testPerm.model.validateModel();
for i = 1:numel(statesPerm)
    statesPerm{i}.perm = model.getProp(statesPerm{i}, 'Permeability');
end

%%
close all
test.plot(states);
testPerm.plot(statesPerm);

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