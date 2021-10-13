function varargout = simulationSolverFun(problem, objective, varargin)
opt = struct('computeGradient',          false, ...
             'objectiveHandler',            [], ...
             'parameters',                  [], ...
             'adjointLinearSolver',         [], ...
             'clearStatesAfterAdjoint', false,  ...
             'maps',                        []);
opt = merge_options(opt, varargin{:});         
         
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
obj.value = -sum(vertcat(vals{:}));

% compute gradient
if opt.computeGradient || nargout == 2
    objh = @(tstep, model, state)objective(setup.model, states, setup.schedule, ...
           'ComputePartials', true, 'tStep', tstep,'state', state, 'from_states', false);
    if isempty(opt.parameters) % assume standard type
        gradient = computeGradientAdjointAD(setup.state0, states, setup.model, ...
                   setup.schedule, objh, 'OutputPerTimestep', true);
        ws  = problem.OutputHandlers.wellSols;
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
    if opt.clearStatesAfterAdjoint
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