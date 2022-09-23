function varargout = simulationSolverFun(problem, objective, varargin)
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

opt = struct('computeGradient',          false, ...
             'objectiveHandler',            [], ...
             'parameters',                  [], ...
             'adjointLinearSolver',         [], ...
             'clearStatesAfterAdjoint', false,  ...
             'maps',                        [], ...
             'nSteps',                      [], ...
             'scalarObjective',           true);
[opt, extra] = merge_options(opt, varargin{:});         
         
if ~isempty(opt.nSteps)
    problem.SimulatorSetup.schedule.step.val = ...
         problem.SimulatorSetup.schedule.step.val(1:opt.nSteps);
    problem.SimulatorSetup.schedule.step.control = ...
         problem.SimulatorSetup.schedule.step.control(1:opt.nSteps);
end

setup   = problem.SimulatorSetup;
% simulate       
if isfield(problem, 'OutputHandlers')
    simulatePackedProblem(problem); 
    states = problem.OutputHandlers.states;
    report = problem.OutputHandlers.reports;
    nStates = states.numelData;
else
    [~, states, report]  = simulateScheduleAD(setup.state0, setup.model, setup.schedule);
    nStates = numel(states);
end
hasMinisteps = isfield(setup, 'OutputMinisteps') && setup.OutputMinisteps;

% compute objective
if iscell(objective)
    % multiple objectives
    objective = objective{problem.seed};
end
% Check for 1-1 between states and schedule and handle mini-steps
objArg = {};
if nStates ~= numel(setup.schedule.step.val) 
    if hasMinisteps
        [setup.schedule, stepMap] = includeMiniSteps(setup.schedule, report, states.numelData);
        % For parameter calibration, provide a map from steps to reportstep.
        % It's up to the objective whether this is treated adequately
        if ~isempty(opt.parameters)
            objArg = {'reportStepMap', stepMap};
        end
    else
        error('Unable to resolve incompatibility between states end schedule');
    end
end

vals = objective(setup.model, states, setup.schedule, objArg{:});

if opt.scalarObjective
    obj.value = sum(vertcat(vals{:}));
else
    % assumed to be of ls-type, so return non-squared
    obj.value = (vertcat(vals{:})).^(1/2);
end

% compute gradient
if opt.computeGradient || nargout == 2
    objh = @(tstep, model, state)objective(model, states, setup.schedule, ...
           'ComputePartials', true, 'tStep', tstep,'state', state, ...
           'from_states', false, objArg{:});
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
                                                     'LinearSolver',opt.adjointLinearSolver, ...
                                                     'isScalar', opt.scalarObjective, extra{:});
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
        if opt.scalarObjective
            obj.gradient = vertcat(scaledGradient{:});
        else
            J = vertcat(scaledGradient{:})';
            % adjust jacobian for sum r^2 -> sum r
            nzIx  = abs(obj.value) > eps*norm(obj.value);
            numnz = nnz(nzIx);
            J(nzIx,:) = spdiags((2*obj.value(nzIx)), 0, numnz, numnz)\J(nzIx,:);
            obj.gradient = J;
        end
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

function [schedule, stepMap] = includeMiniSteps(schedule, report, nStates)
if isa(report, 'ResultHandler')
    tmp  = cell(report.numelData, 1);
    for k = 1:report.numelData
        tmp{k} = report{k};
    end
    report = tmp;
end
[schedule, ts, nSteps] = convertReportToSchedule(report, schedule);
assert(nStates == sum(nSteps), ...
    'Number of ministeps in report did not match number of states');
stepMap = nan(numel(ts),1);
stepMap(cumsum(nSteps)) = (1:numel(nSteps))';
end

