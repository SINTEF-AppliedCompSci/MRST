function sens = computeSensitivitiesTwinModelsAdjointAD(setup_1, states_1, setup_2, states_2, param, getObjective, varargin)
% Compute parameter sensitivities for twin models using adjoint simulation
%
% SYNOPSIS:
%   sens = computeSensitivitiesTwinModelsAdjointAD(setup_1, states_1, setup_2, states_2, param, getObjective)
%
% PARAMETERS:
%   setup_1, setup_2 - Simulation setups for both models
%   states_1, states_2 - State histories for both models
%   param - Cell array of ModelParameter objects
%   getObjective - Function handle for objective evaluation
%
% RETURNS:
%   sens - Structure containing sensitivity values for each parameter

% Copyright information
opt = struct('LinearSolver',          [], ...
             'isScalar',            true, ...
             'accumulateResiduals',   [], ...
             'matchMap',              []);
opt = merge_options(opt, varargin{:});

% Initialize sensitivity structure
sens = struct();
for k = 1:numel(param)
    sens.(param{k}.name) = 0;
end

% Setup linear solver
if ~isempty(opt.LinearSolver)
    linsolve = opt.LinearSolver;
elseif isfield(setup_2, 'AdjointLinearSolver')
    linsolve = setup_2.AdjointLinearSolver;
else
    linsolve = BackslashSolverAD();
end

% Split initial state parameters
isInitParam = cellfun(@(p) strcmp(p.belongsTo, 'state0'), param);
if any(isInitParam)
    [initparam, param] = deal(param(isInitParam), param(~isInitParam));
end

% Validate models and prepare AD parameters
[setup_1.model, setup_1.schedule] = validateModelSchedule(setup_1);
[setup_2.model, setup_2.schedule] = validateModelSchedule(setup_2);

[modelParam_1, scheduleParam_1] = initModelParametersADI(setup_1, param);
[modelParam_2, scheduleParam_2] = initModelParametersADI(setup_2, param);
modelParam_1 = modelParam_1.removeStateFunctionGroupings();
modelParam_1 = validateModel(modelParam_1);
modelParam_2 = modelParam_2.removeStateFunctionGroupings();
modelParam_2 = validateModel(modelParam_2);
% Prepare state access functions
getState_1 = @(i) getStateFromInput(setup_1.schedule, states_1, setup_1.state0, i);
getState_2 = @(i) getStateFromInput(setup_2.schedule, states_2, setup_2.state0, i);

% Initialize adjoint variables
if ~isempty(opt.accumulateResiduals) || ~opt.isScalar
    [colIx, nrow, ncol] = getColumnIndex(opt.accumulateResiduals, opt.matchMap, setup_1, states_1{end});
    lambda1 = zeros(nrow, ncol);
    lambda2 = zeros(nrow, ncol);
else
    colIx = repmat({nan}, [numel(setup_1.schedule.step.val), 1]);
    lambda1 = [];
    lambda2 = [];
end

% Run adjoint simulation
nstep_1 = numel(setup_1.schedule.step.val);
for step = nstep_1:-1:1
    fprintf('Solving reverse mode step %d of %d\n', nstep_1 - step + 1, nstep_1);
    
    % Solve adjoint equations for both models
    [lami1, lambda1, ~, lami2, lambda2, ~] = ...
        setup_1.model.solveAdjointTwinModels(setup_1.model, linsolve, ...
        getState_1, getState_2, getObjective, setup_1.schedule, ...
        setup_2.schedule, lambda1, lambda2, step, 'colIx', colIx{step});

    % Compute partial derivatives
    [eqdth1, modelParam_1] = partialWRTparam(modelParam_1, getState_1, scheduleParam_1, step, param);
    [eqdth2, modelParam_2] = partialWRTparam(modelParam_2, getState_2, scheduleParam_2, step, param);

    % Accumulate sensitivities
    result = 0;
    for k = 1:numel(lami1)
        result = result + lami1{k}'*eqdth1{k} + lami2{k}'*eqdth2{k};
    end
    
    result = result.jac;
    if numel(result) ~= numel(param)
        result = horzcat(result{:});
        cumn = reshape(cumsum(cellfun(@(p) p.nParam, param)), [], 1);
        [lo, hi] = deal([1; cumn(1:end-1)+1], cumn);
        result = applyFunction(@(i1, i2) result(i1:i2), lo, hi);
    end
    
    for k = 1:numel(param)
        sens.(param{k}.name) = sens.(param{k}.name) + result{k}.';
    end
end

% Handle initial state parameters
if any(isInitParam)
    sens = computeInitialStateSensitivities(setup_1, states_1, initparam, sens, lami1);
end
end

%% Helper Functions
function [model, schedule] = validateModelSchedule(setup)
    model = validateModel(setup.model);
    schedule = setup.model.validateSchedule(setup.schedule);
end

function sens = computeInitialStateSensitivities(setup, states, initparam, sens, lami)
    schedule = setup.schedule;
    forces = setup.model.getDrivingForces(schedule.control(schedule.step.control(1)));
    forces = merge_options(setup.model.getValidDrivingForces(), forces{:});
    model = setup.model.validateModel(forces);
    
    state0 = model.validateState(setup.state0);
    state0.wellSol = states{1}.wellSol;
    state0 = model.getStateAD(state0);
    dt = schedule.step.val(1);
    
    linProblem = model.getAdjointEquations(state0, states{1}, dt, forces,...
        'iteration', inf, 'reverseMode', true);
    
    pNames = fieldnames(sens);
    nms = applyFunction(@lower, pNames(cellfun(@(p) strcmp(p.belongsTo, 'state0'), initparam)));
    varNms = applyFunction(@lower, linProblem.primaryVariables);
    
    for k = 1:numel(nms)
        kn = find(strcmp(nms{k}, varNms));
        for nl = 1:numel(lami)
            if isa(linProblem.equations{nl}, 'ADI')
                sens.(nms{k}) = sens.(nms{k}) + linProblem.equations{nl}.jac{kn}'*lami{nl};
            end
        end
        if strcmp(initparam{k}.type, 'multiplier')
            sens.(nms{k}) = sens.(nms{k}).*initparam{k}.referenceValue;
        end
        sens.(nms{k}) = initparam{k}.collapseGradient(sens.(nms{k}));
    end
end
function [eqdth, model] = partialWRTparam(model, getState, schedule, step, params)
validforces = model.getValidDrivingForces();
current = getState(step);
before  = getState(step - 1);
dt_steps = schedule.step.val;
dt       = dt_steps(step);
cNo      = schedule.step.control(step);

control = schedule.control(cNo);
forces  = model.getDrivingForces(control);
forces  = merge_options(validforces, forces{:});

reValidate = step == numel(schedule.step.val) || ...
             cNo ~= schedule.step.control(step+1);
if reValidate
    if isa(model, 'ReservoirModel') && ~isempty(model.FacilityModel)
        model.FacilityModel = model.FacilityModel.validateModel(forces);
    else
        model = model.validateModel(forces);
    end
end
% Initial state typically lacks wellSol-field, so add if needed
if step == 1
    before = model.validateState(before);
end
% initialize before-state in case it contains cached properties
before  = model.getStateAD(before, false);
problem = model.getEquations(before, current, dt, forces, 'iteration', inf, 'resOnly', true);

% We need special treatment of well-control parameters (if present) due to 
% non-diff misc logic in well equations
isPolicy      = cellfun(@(p)strcmp(p.controlType, 'policy'), params);
isWellControl = cellfun(@(p)any(strcmp(p.controlType, ...
                {'bhp', 'rate', 'wrat', 'orat', 'grat'})), params);
if any(isPolicy | isWellControl)
    ws     = current.wellSol;
    %nw     = numel(ws);
    isOpen = vertcat(ws.status);
    nOpen  = nnz(isOpen);
    eqNo   = strcmp('closureWells', problem.equationNames);
    ceq    = problem.equations{eqNo};
    assert(~isa(ceq, 'ADI'), 'Expected equation to be of class ADI');
    if any(isPolicy) && step > 1
        % experimental support for policies: 
        assert(isfield(control.misc, 'policy'));
        tmp = control.extra.policy.function(before, schedule, [], step-1);
        ceq = -vertcat(tmp.control(cNo).W(isOpen).val);
    end
    if any(isWellControl)
        % set to AD if not already
        ceq = setEqToADI(ceq, params);
        ix = find(isWellControl);
        for k = 1:numel(ix)
            pnum = ix(k);
            if ismember(cNo, params{pnum}.controlSteps)
                assert(strcmp(params{pnum}.type, 'value'), ...
                    'Well control parameter of non-value type not supported');
                cntr     = params{pnum}.controlType;
                % non-zero gradient only if control is active
                isActive = isOpen & reshape(strcmp(cntr, {ws.type}), [], 1);
                jac = - spdiags(double(isActive(isOpen)), 0, nOpen, nOpen);
                jac = jac(:, params{pnum}.subset);
                assert(all(size(ceq.jac{pnum}) == size(jac)), 'Problem requires debugging')
                ceq.jac{pnum} = jac;
            end
        end
    end
    % scaling (identical to scaling of well eqs)
    scale = getControlEqScaling({ws.type}, model.FacilityModel);
    scale = scale(isOpen);
    problem.equations{eqNo} = scale.*ceq;
end             
eqdth = problem.equations;
end

%--------------------------------------------------------------------------

function state = getStateFromInput(schedule, states, state0, i)
if i == 0
    state = state0;
elseif i > numel(schedule.step.val)
    state = [];
else
    state = states{i};
end
end
%--------------------------------------------------------------------------

function [modelParam, scheduleParam] = initModelParametersADI(setup, param)
v  = applyFunction(@(p)p.getParameter(setup), param);
if isprop(setup.model.AutoDiffBackend, 'useMex') && setup.model.AutoDiffBackend.useMex
    % mex-version of operators are not compatible
    setup.model.AutoDiffBackend.useMex = false;
    % operators are updated in subsequent call to validateModel
end
% use same backend as setup.model
if isfield(setup, 'model') && isprop(setup.model, 'AutoDiffBackend')
    [v{:}] = setup.model.AutoDiffBackend.initVariablesAD(v{:});
else
    [v{:}] = initVariablesADI(v{:});
end
for k = 1:numel(v)
    % don't set for well controls
    if ~any(strcmp(param{k}.controlType, {'bhp', 'rate', 'wrat', 'orat', 'grat'}))
        setup = param{k}.setParameter(setup, v{k});
    end
end
[modelParam, scheduleParam] = deal(setup.model, setup.schedule);
end
%--------------------------------------------------------------------------
function eq = setEqToADI(eq, p)
if ~isa(eq, 'ADI')
    % we don't necessarily have sample, so create from scratch 
    np  = cellfun(@(p)p.nParam, p);
    neq = numel(eq);
    jac = applyFunction(@(npk)sparse([],[],[],neq, npk ), np);
    % use default ADI
    eq  = ADI(eq, jac);
end
end

%--------------------------------------------------------------------------
function sc = getControlEqScaling(tp, facility)
% In order to produce correct gradients, this scaling must be precisely the
% same as in the control equation. This is rather shaky!
sc    = ones(numel(tp), 1);
is_bhp = strcmp('bhp', tp);
if any(is_bhp)
    if isfinite(facility.toleranceWellRate) && isfinite(facility.toleranceWellBHP)
        eqScale = facility.toleranceWellRate/facility.toleranceWellBHP;
    else
        eqScale = 1/(day()*barsa());
    end
    sc(is_bhp) = eqScale;
end
end

%--------------------------------------------------------------------------
function [colIx, nrow, ncol] = getColumnIndex(accum, matchMap, setup, state)
if ~isempty(accum) && ~isempty(matchMap)
    warning('Accumulation of residuals will be ignored (not compatible with ''matchMap''-option)');
end
ns  = numel(setup.schedule.step.val);
nw  = numel(setup.schedule.control(1).W);
nt  = setup.model.water + setup.model.oil + setup.model.gas + 1;
if isprop(setup.model, 'thermal') && setup.model.thermal
    nt = nt +1;
end
% we need a lambda corresponding to the last step
%nwo = nnz(vertcat(setup.schedule.control(end).W.status));
nwo = nnz(vertcat(state.wellSol.status));
nrow = (nt-1)*setup.model.G.cells.num + nt*nwo; % shaky
if isempty(matchMap)
    if isempty(accum) || isempty(accum.wells)
        accum.wells = (1:nw)';
    end
    assert(numel(accum.wells)==nw, 'Mismatch: number of wells');
    if ~isfield(accum, 'types') || isempty(accum.types)
        accum.types = (1:nt)';
    end
    assert(numel(accum.types)==nt, 'Mismatch: number of types (rates/bhp)');
    if ~isfield(accum, 'steps') || isempty(accum.steps)
        accum.steps = (1:ns)';
    end
    assert(numel(accum.steps)==ns, 'Mismatch: number of time steps');
    % build collumn index
    [mw, mt] = deal(max(accum.wells), max(accum.types));
    colIx = applyFunction(@plus, repmat({(1:2*mw*(mt)+ mw + 3)'}, [ns, 1]), ...
                          num2cell((2*mw*(mt) +mw+3)*(accum.steps-1)));
    [colIx{accum.steps<=0}] = deal({});
    ncol = max(colIx{end});
else
    nMatch = zeros(ns,1);
    fn = fieldnames(matchMap.steps);
    assert(ns == numel(matchMap.steps), ...
        'Dimension mismatch between schedule and matchMap');
    for k = 1:numel(fn)
        nMatch = nMatch +  reshape(cellfun(@(x)nnz(isfinite(x)), {matchMap.steps.(fn{k})}), [], 1);
    end
    cum = cumsum(nMatch);
    colIx = applyFunction(@(x,y)x:y, [0; cum(1:end-1)]+1, cum);
    ncol  = cum(end);
end
end