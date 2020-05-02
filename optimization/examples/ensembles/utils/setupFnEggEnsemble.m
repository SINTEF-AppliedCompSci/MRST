function problem = setupFnEggEnsemble(memberIx)

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
