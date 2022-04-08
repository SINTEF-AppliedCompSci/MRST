function problem = setupFnEggEnsemble(seed, directory, subset)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
if nargin < 3
    memberIx = seed;
else
    memberIx = subset(seed);
end
if nargin < 2
    directory = '';
end

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

nm = ['egg', num2str(seed)];

problem = packSimulationProblem(state0, model, schedule, nm , ...
                                'Name',             num2str(seed), ...
                                'NonLinearSolver',  nonlinear, ...
                                'Directory',        directory, ...
                                'ExtraArguments',   {'controlLogicFn', @exampleLogicFunc});
problem.seed = seed;                            

end


function [schedule, report, isAltered] = exampleLogicFunc(state, schedule, report, cnum)
[schedule, report, isAltered] = controlLogicFunc(state, schedule, report, cnum, ...
                                                 'wcutLim', .95, ...
                                                 'closeWellOnSignChange', true, ...
                                                 'rateLim', .1/day);
end
