mrstModule add ad-unittest ad-core ad-blackoil

%% Run SPE 1
tests = matlab.unittest.TestSuite.fromClass(?TestSPE1);
result_spe1 = run(tests);

%% Run SPE 3
tests = matlab.unittest.TestSuite.fromClass(?TestSPE3);
result_spe3 = run(tests);

%% Run SPE 9
tests = matlab.unittest.TestSuite.fromClass(?TestSPE9);
result_spe9 = run(tests);

%% Run oil/water test
tests = matlab.unittest.TestSuite.fromClass(?TestSimpleOW);
result_ow = run(tests);

%% Run oil/water test
tests = matlab.unittest.TestSuite.fromClass(?TestEGG);
result_egg = run(tests);

%% Run tests on the nonlinear solver
tests = matlab.unittest.TestSuite.fromClass(?TestSimulateScheduleAD);
result_sched = run(tests);

%%
clear
tests = matlab.unittest.TestSuite.fromClass(?ResultHandlerTest);
res = run(tests);

%%
tests = matlab.unittest.TestSuite.fromClass(?TestAdjoints);
result_adjoints = run(tests);

%% Run all simulation tests
folder = mrstPath('query', 'ad-unittest');

tests = matlab.unittest.TestSuite.fromFolder(fullfile(folder, 'test_sim'));
results = run(tests);

%%
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
