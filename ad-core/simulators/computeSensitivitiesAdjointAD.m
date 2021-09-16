function sens = computeSensitivitiesAdjointAD(setup, states, param, getObjective, varargin)
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
%   'LinearSolver'   - Subclass of `LinearSolverAD` suitable for solving the
%                      adjoint systems.
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


assert(isa(setup.model,'GenericReservoirModel'),... 
       'The model must be derived from GenericReservoirModel.')
assert(isa(param{1}, 'ModelParameter'), ...
        'Parameters must be initialized using ''ModelParameter''.')

opt = struct('LinearSolver', []);       
opt = merge_options(opt, varargin{:});
if mrstVerbose && setup.model.toleranceCNV >= 1e-3
   fprintf(['The accuracy in the gradient depend on the',...
            'acuracy on the CNV tolerance.\n',...
            'For better accuracy set a lower value for '...
            'model.toleranceCNV.'] )
end

getState = @(i) getStateFromInput(setup.schedule, states, setup.state0, i);

if isempty(opt.LinearSolver)
    linsolve = BackslashSolverAD();
else
    linsolve = opt.LinearSolver;
end

sens = struct;
for k = 1:numel(param)
    sens.(param{k}.name) = 0;
end
pNames = fieldnames(sens);
isInitStateParam = cellfun(@(p)strcmp(p.belongsTo, 'state0'), param);
% inititialize parameters to ADI
[modelParam, scheduleParam] = initModelParametersADI(setup, param(~isInitStateParam));
% reset discretization/flow functions to account for AD-parameters
modelParam.FlowDiscretization = [];
modelParam.FlowPropertyFunctions = [];
modelParam = validateModel(modelParam);

nstep = numel(setup.schedule.step.val);
lambda = [];
nt = nstep;
for step = nt:-1:1
    fprintf('Solving reverse mode step %d of %d\n', nt - step + 1, nt);
    [lami, lambda]= setup.model.solveAdjoint(linsolve, getState, ...
        getObjective, setup.schedule, lambda, step);
    eqdth = partialWRTparam(modelParam, getState, scheduleParam, step);
    for kp = 1:numel(param)
        if ~isInitStateParam(kp)
            nm = param{kp}.name;
            for nl = 1:numel(lami)
                if isa(eqdth{nl}, 'ADI')
                    sens.(nm) = sens.(nm) + eqdth{nl}.jac{kp}'*lami{nl};
                end
            end
        end        
    end
end

% compute partial derivative of first eq wrt init state and compute
% initial state sensitivities
if any(isInitStateParam)
    schedule = setup.schedule;
    forces = setup.model.getDrivingForces(schedule.control(schedule.step.control(1)));
    forces = merge_options(setup.model.getValidDrivingForces(), forces{:});
    model = setup.model.validateModel(forces);
    
    state0 = model.validateState(setup.state0);
    % set wellSols just to make subsequent function-calls happy, sensitivities wrt wellSols doesn't make sense anyway
    state0.wellSol = states{1}.wellSol;
    dt = schedule.step.val(1);
    
    linProblem = model.getAdjointEquations(state0, states{1}, dt, forces,...
        'iteration', inf, 'reverseMode', true);
    nms    = applyFunction(@lower, pNames(isInitStateParam));
    varNms = applyFunction(@lower, linProblem.primaryVariables);
    for k = 1:numel(nms)
        kn = find(strcmp(nms{k}, varNms));
        assert(numel(kn)==1, 'Unable to match initial state parameter name %s\n', nms{k});
        for nl = 1:numel(lami)
            if isa(linProblem.equations{nl}, 'ADI')
                sens.(nms{k}) = sens.(nms{k}) + linProblem.equations{nl}.jac{kn}'*lami{nl};
            end
        end
        kp = strcmp(nms{k}, pNames);
        sens.(nms{k}) =  param{kp}.collapseGradient(sens.(nms{k}));
    end
end       

end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function eqdth = partialWRTparam(model, getState, schedule, step)
validforces = model.getValidDrivingForces();
current = getState(step);
before  = getState(step - 1);
dt_steps = schedule.step.val;
dt       = dt_steps(step);
cNo      = schedule.step.control(step);

control = schedule.control(cNo);
forces  = model.getDrivingForces(control);
forces  = merge_options(validforces, forces{:});
model   = model.validateModel(forces);
% Initial state typically lacks wellSol-field, so add if needed
if step == 1
    before = model.validateState(before);
end
% initialize before-state in case it contains cached properties
before  = model.getStateAD(before, false);
problem = model.getEquations(before, current, dt, forces, 'iteration', inf, 'resOnly', true);
% experimental support for policies:
if isfield(control, 'policy') && step > 1
    % direct-set jacobian to avoid misc logic in well equations
    % NOTE: fix scaling of bhp-controls
    eqNo  = strcmp('closureWells', problem.equationNames);
    isAct = vertcat(current.wellSol.status);
    tmp   = control.policy.function(before, schedule, [], step-1);
    problem.equations{eqNo} = -vertcat(tmp.control(cNo).W(isAct).val);
end
eqdth   = problem.equations;
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
% use same backend as problem.model
if isfield(setup, 'model') && isprop(setup.model, 'AutoDiffBackend')
    [v{:}] = setup.model.AutoDiffBackend.initVariablesAD(v{:});
else
    [v{:}] = initVariablesADI(v{:});
end
for k = 1:numel(v)
    setup = param{k}.setParameter(setup, v{k});
end
[modelParam, scheduleParam] = deal(setup.model, setup.schedule);
end
