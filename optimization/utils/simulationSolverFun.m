function varargout = simulationSolverFun(problem, objective, varargin)
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

opt = struct('computeGradient',          false, ...
             'objectiveHandler',            [], ...
             'parameters',                  [], ...
             'adjointLinearSolver',         [], ...
             'clearStatesAfterAdjoint', false,  ...
             'maps',                        [], ...
             'nSteps',                      []);
opt = merge_options(opt, varargin{:});         
         
if ~isempty(opt.nSteps)
    problem.SimulatorSetup.schedule.step.val = ...
         problem.SimulatorSetup.schedule.step.val(1:opt.nSteps);
    problem.SimulatorSetup.schedule.step.control = ...
         problem.SimulatorSetup.schedule.step.control(1:opt.nSteps);
end

% simulate       
if isfield(problem, 'OutputHandlers')
    simulatePackedProblem(problem); 
    states  = problem.OutputHandlers.states;
else
    setup = problem.SimulatorSetup;
    [~, states]  = simulateScheduleAD(setup.state0, setup.model, setup.schedule);
end

% compute objective
setup   = problem.SimulatorSetup;
if iscell(objective)
    % multiple objectives
    objective = objective{problem.seed};
end
vals = objective(setup.model, states, setup.schedule);
obj.value = sum(vertcat(vals{:}));

% compute gradient
if opt.computeGradient || nargout == 2
    objh = @(tstep, model, state)objective(model, states, setup.schedule, ...
           'ComputePartials', true, 'tStep', tstep,'state', state, 'from_states', false);
    if isempty(opt.parameters) % assume standard type
        gradient = computeGradientAdjointAD(setup.state0, states, setup.model, ...
                   setup.schedule, objh, 'OutputPerTimestep', true);
        if isfield(problem, 'OutputHandlers')
            ws  = problem.OutputHandlers.wellSols;
        else
            ws = applyFunction(@(s)s.wellSol, states);
        end
        obj.gradient = processAdjointGradients(gradient, ws, 'controlIx', ...
                                                setup.schedule.step.control);
        if ~isempty(opt.maps)
            obj.gradient = mapGradient(obj.gradient, opt.maps);
        end
    else
        params = opt.parameters;
        gradient = computeSensitivitiesAdjointAD(setup, states, opt.parameters, objh,...
                                                     'LinearSolver',opt.adjointLinearSolver);
        %scaling
        nms = applyFunction(@(x)x.name, params);
        scaledGradient = cell(numel(nms), 1);
        for k = 1:numel(nms)
            pval = [];
            if ~strcmp(params{k}, 'linear')
                pval = params{k}.getParameter(setup);
            end
            scaledGradient{k} = params{k}.scaleGradient( gradient.(nms{k}), pval);
        end
        obj.gradient = vertcat(scaledGradient{:});
    end
    if opt.clearStatesAfterAdjoint && isa(states, 'ResultHandler')
        states.resetData();
    end
end

% return values
if ~isempty(opt.objectiveHandler)
    if isfield(problem, 'seed')
        n = problem.seed;
    else
        warning('Problem is lacking a ''seed''-field, assuming single realization scenario.')
        n = 1;
    end
    opt.objectiveHandler{n} = obj;
end
if nargout > 0
    varargout{1} = obj.value;
end
if nargout > 1
    varargout{1} = obj.gradient;
end
end

function g = mapGradient(grad, maps)
g = zeros(numel(maps.type), 1);
for k = 1:numel(g)
    [sno, wno, tp] = deal(maps.stepNo(k), maps.wellNo(k), maps.type{k});
    g(k) = grad.(tp)(wno, sno);
end
end
