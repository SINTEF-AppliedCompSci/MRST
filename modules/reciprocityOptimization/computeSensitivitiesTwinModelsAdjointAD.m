function sens = computeSensitivitiesTwinModelsAdjointAD(setup_1, states_1,setup_2, states_2, param, getObjective, varargin)
%Compute parameter sensitivities using adjoint simulation
%
% SYNOPSIS:
%   sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, getObjective, varargin)
%
% DESCRIPTION:
%   For a given schedule, compute senistivities with regards to parameters
%
% REQUIRED PARAMETERS:
%
%   SimulatorSetup - structure containing:
%       state0   - Physical model state at `t = 0`
%       model    - Subclass of PhysicalModel class such as
%                  `GenericBlackOilModel` that models the physical
%                  effects we want to study.
%
%       schedule - Schedule suitable for `simulateScheduleAD`.
%
%   states       - All previous states. Must support the syntax
%                  `state = states{i}`. If the problem is too large to fit in
%                  memory, it can use `ResultHandler` class to retrieve files
%                  from the disk.
%
%   param        - array of parameters of class ModelParameter
%
%   getObjective - Function handle for getting objective function value
%                  for a given timestep with derivatives. Format: @(tstep)
%
% OPTIONAL PARAMETERS:
%
%   'LinearSolver'           - Subclass of `LinearSolverAD` suitable for 
%                              solving the adjoint systems.
%   'isScalar'               - boolean indicating whether objective is scalar 
%                             (default true)
%   'accumulateResiduals'    - experimental/undocumented feature
%
% RETURNS:
%   sens - Structure with parameter sensitivites of the form
%          sens.(paramName) = paramValue
%
% SEE ALSO:
%   `computeGradientAD`

%{
Copyright 2009-2021 SINTEF Digital, Applied Mathematics.

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

if isa(setup_1.model, 'ReservoirModel')
    assert(isa(setup_1.model,'GenericReservoirModel'),... 
       'The model must be derived from GenericReservoirModel.')
end

opt = struct('LinearSolver',          [], ...
             'isScalar',            true, ...
             'accumulateResiduals',   [], ...
             'matchMap',              []);       
opt = merge_options(opt, varargin{:});

if ~isempty(opt.LinearSolver)
    linsolve = opt.LinearSolver;
elseif isfield(setup_2, 'AdjointLinearSolver')
    linsolve = setup_2.AdjointLinearSolver;
else
    linsolve = BackslashSolverAD();
end
    
sens = struct;
for k = 1:numel(param)
    sens.(param{k}.name) = 0;
end
pNames = fieldnames(sens);
% split parameters in inital state/non-initial state due to different
% handling
isInitParam = cellfun(@(p)strcmp(p.belongsTo, 'state0'), param);
if any(isInitParam)
    [initparam, param] = deal(param(isInitParam), param(~isInitParam));
end
% validate simulation model:
setup_1.model = validateModel(setup_1.model);
setup_1.schedule = setup_1.model.validateSchedule(setup_1.schedule);

% validate simulation model:
setup_2.model = validateModel(setup_2.model);
setup_2.schedule = setup_2.model.validateSchedule(setup_2.schedule);

% inititialize parameters to ADI
[modelParam_1, scheduleParam_1] = initModelParametersADI(setup_1, param);
% reset discretization/flow functions to account for AD-parameters
modelParam_1 = modelParam_1.removeStateFunctionGroupings();
        
modelParam_1 = validateModel(modelParam_1);

% inititialize parameters to ADI
[modelParam_2, scheduleParam_2] = initModelParametersADI(setup_2, param);
% reset discretization/flow functions to account for AD-parameters
modelParam_2 = modelParam_2.removeStateFunctionGroupings();
        
modelParam_2 = validateModel(modelParam_2);

nstep_1    = numel(setup_1.schedule.step.val);
getState_1 = @(i) getStateFromInput(setup_1.schedule, states_1, setup_1.state0, i);

nstep_2    = numel(setup_2.schedule.step.val);
getState_2 = @(i) getStateFromInput(setup_2.schedule, states_2, setup_2.state0, i);

if ~isempty(opt.accumulateResiduals) || ~opt.isScalar
    % allocate correct size of Lagrange multiplier matrix
    [colIx, nrow, ncol] = getColumnIndex(opt.accumulateResiduals, opt.matchMap, setup_1, states_1{end});
    lambda1 = zeros(nrow, ncol);
    [colIx, nrow, ncol] = getColumnIndex(opt.accumulateResiduals, opt.matchMap, setup_2, states_2{end});
    lambda2 = zeros(nrow, ncol);
else
    colIx = repmat({nan}, [nstep_1, 1]);
    lambda1 = [];

    lambda2 = [];

end

% run adjoint
for step = nstep_1:-1:1
    fprintf('Solving reverse mode step %d of %d\n', nstep_1 - step + 1, nstep_1);
    [lami1, lambda1,rep1,lami2,lambda2,rep2]= setup_1.model.solveAdjointTwinModels(setup_2.model,linsolve, getState_1,getState_2, ...
        getObjective, setup_1.schedule,setup_2.schedule, lambda1,lambda2, step, 'colIx', colIx{step});
    [eqdth1, modelParam_1] = partialWRTparam(modelParam_1, getState_1, scheduleParam_1, step, param);

    [eqdth2, modelParam_2] = partialWRTparam(modelParam_2, getState_2, scheduleParam_2, step, param);
%%% I'm not sure but I think equations order may be different 
% cc=lami2;
% lami2{3}=lami2{4};
% lami2{4}=cc{3};

    result = 0;
    for k = 1:numel(lami1)
        result = result + lami1{k}'*eqdth1{k}+lami2{k}'*eqdth2{k};
    end
    result = result.jac;
    if numel(result) ~= numel(param) % might be the case for e.g., GenericAD
        result = horzcat(result{:});
        cumn   = reshape(cumsum(cellfun(@(p)p.nParam, param)), [], 1);
        [lo, hi] = deal([1; cumn(1:end-1)+1], cumn);
        result = applyFunction(@(i1, i2)result(i1:i2), lo, hi);
    end
    for k = 1:numel(param)
        sens.(param{k}.name) = sens.(param{k}.name) + result{k}.';
    end
end

% compute partial derivative of first eq wrt init state and compute
% initial state sensitivities
if any(isInitParam)
    schedule = setup_1.schedule;
    forces = setup_1.model.getDrivingForces(schedule.control(schedule.step.control(1)));
    forces = merge_options(setup_1.model.getValidDrivingForces(), forces{:});
    model = setup_1.model.validateModel(forces);
    
    state0 = model.validateState(setup_1.state0);
    % set wellSols just to make subsequent function-calls happy, sensitivities wrt wellSols doesn't make sense anyway
    state0.wellSol = states_1{1}.wellSol;
    state0 = model.getStateAD(state0);
    dt = schedule.step.val(1);
    
    linProblem = model.getAdjointEquations(state0, states_1{1}, dt, forces,...
        'iteration', inf, 'reverseMode', true);
    nms    = applyFunction(@lower, pNames(isInitParam));
    varNms = applyFunction(@lower, linProblem.primaryVariables);
    for k = 1:numel(nms)
        kn = find(strcmp(nms{k}, varNms));
        assert(numel(kn)==1, 'Unable to match initial state parameter name %s\n', nms{k});
        for nl = 1:numel(lami)
            if isa(linProblem.equations{nl}, 'ADI')
                sens.(nms{k}) = sens.(nms{k}) + linProblem.equations{nl}.jac{kn}'*lami{nl};
            end
        end
        if strcmp(initparam{k}.type, 'multiplier')
            sens.(nms{k}) = sens.(nms{k}).*initparam{k}.referenceValue;
        end
        sens.(nms{k}) =  initparam{k}.collapseGradient(sens.(nms{k}));
    end
end       

end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
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
    colIx = applyFunction(@plus, repmat({(1:mw*(mt+2)+6)'}, [ns, 1]), ...
                          num2cell((mw*(mt+2)+6)*(accum.steps-1)));
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
