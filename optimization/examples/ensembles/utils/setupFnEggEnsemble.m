function problem = setupFnEggEnsemble(memberIx)
%Undocumented Utility Function

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

[G, ~, ~, deck] = setupEGG('realization', memberIx);
[state0, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none', 'useMex', true);

% The schedule only has a single control-step lasting approx 10 years. We 
% reduce the simulation period to ~five years and devide it into 10 
% initially equal control steps 

t  = cumsum( schedule.step.val );
ix = t <= 5*year;
t  = t(ix);

controlNo = ceil( 10*t/t(end)); 
control   = schedule.control(1);
schedule.step.control = controlNo;
schedule.step.val     = schedule.step.val(ix);
schedule.control = repmat(control, [1, max(controlNo)]);


nonlinear.useLinesearch = true;

% We set up a control-logic function below such that wells are shut down
% whenever
%   * water-cut exceeds 0.95
%   * producers become injectors and vice versa
%   * well-rates drop below 0.1 m^3/day

problem = packSimulationProblem(state0, model, schedule, 'tmp', ...
                                'Name',             num2str(memberIx), ...
                                'NonLinearSolver',  nonlinear, ...
                                'ExtraArguments',   {'controlLogicFn', @exampleLogicFunc});

end


function [schedule, report, isAltered] = exampleLogicFunc(state, schedule, report, cnum)
[schedule, report, isAltered] = controlLogicFunc(state, schedule, report, cnum, ...
                                                 'wcutLim', .95, ...
                                                 'closeWellOnSignChange', true, ...
                                                 'rateLim', .1/day);
end
