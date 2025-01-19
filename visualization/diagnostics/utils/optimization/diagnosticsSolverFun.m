function varargout = diagnosticsSolverFun(problem, objective, varargin)
opt = struct('computeGradient',          false, ...
             'objectiveHandler',            [], ...
             'parameters',                  [], ...
             'adjointLinearSolver',         [], ...
             'maps',                        [], ...
             'clearStatesAfterAdjoint', false, ...
             'W',                          []);
opt = merge_options(opt, varargin{:});         
         
computeGradient = opt.computeGradient || nargout == 2;

% compute
W   = problem.SimulatorSetup.schedule.control(1).W;
objModel = objective(problem.SimulatorSetup.model);

res = objModel.compute(W, 'computeGradient', computeGradient);
obj.value = res.value;

if computeGradient
    g = res.gradient;
    if isfield(g, 'well') && ~isempty(g.well) && ~isempty(opt.maps)
        gw = processAdjointGradients({g.well}, {res.state.wellSol});
        if ~isempty(opt.maps)
            gw = -mapGradient(gw, opt.maps);
        end
    else
        gw = [];
    end
    if isfield(g, 'position') 
        gp = -vertcat(g.position{:});
    else
        gp = [];
    end
    obj.gradient = [gw; gp];
end

% return values
if ~isempty(opt.objectiveHandler)
    if isfield(problem, 'seed')
        n = problem.seed;
    else
        n = str2double(opt.objectiveHandler.dataFolder);
        if ~isfinite(n)
            n = 1;
        end
    end
    opt.objectiveHandler{n} = obj;
end

if nargout > 0
    varargout{1} = obj.value;
end
if nargout == 2
    varargout{2} = obj.gradient;
end
end

function g = mapGradient(grad, maps)
g = zeros(numel(maps.type), 1);
for k = 1:numel(g)
    [sno, wno, tp] = deal(maps.stepNo(k), maps.wellNo(k), maps.type{k});
    g(k) = grad.(tp)(wno, sno);
end
end
