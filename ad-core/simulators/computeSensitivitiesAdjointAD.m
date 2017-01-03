function sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, getObjective, varargin)
%Compute parameter sensitivities using adjoint simulation 
%
% SYNOPSIS:
%   sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, getObjective, varargin)
%
% DESCRIPTION:
%   For a given schedule, compute senistivities with regards to parameters
%  
%
% REQUIRED PARAMETERS:
%
%   state0       - Physical model state at t = 0
%
%   states       - All previous states. Must support the syntax 
%                  state = states{i}. If the problem is too large to fit in
%                  memory, it can use ResultHandler class to retrieve files
%                  from the disk.
%                  
%   model        - Subclass of PhysicalModel class such as
%                 ThreePhaseBlackOilModel that models the physical effects
%                 we want to study.
%
%   schedule     - Schedule suitable for simulateScheduleAD.
%
%   getObjective - Function handle for getting objective function value 
%                  for a given timestep with derivatives. Format: @(tstep)
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   
% 
%
% RETURNS:
%
%  'Parameters'     - Defaults to {'transmissibility', 'porevolumes'} which currenly are the only two options. 
% 
%  'ParameterTypes' - Defaults to 'value' for all. The other option is 'multiplyer'
%
%  'Regularization' - One (sparse) matrix for each parameter-set, requires also modelBase 
%
%  'LinearSolver'   - Subclass of 'LinearSolverAD' suitable for solving the
%                     adjoint systems.
%
% SEE ALSO:
%   computeGradientAD


    opt = struct('Parameters', {{'transmissibility', 'porevolume'}}, ...
                 'ParameterTypes', {{'value', 'value'}}, ...
                 'Regularization', [], ...
                 'LinearSolver',     [], ...
                 'modelBase', []);
    opt = merge_options(opt, varargin{:});
    
    getState = @(i) getStateFromInput(schedule, states, state0, i);
    
    if isempty(opt.LinearSolver)
        linsolve = BackslashSolverAD();
    else
        linsolve = opt.LinearSolver;
    end
    
    param  = opt.Parameters;
    pTypes = opt.ParameterTypes;
    sens = cell(1, numel(param));
    [sens{:}] = deal(0);
    % inititialize parameters to ADI
    modelParam = initModelParametersADI(model, param, pTypes);
    
    nstep = numel(schedule.step.val);
    lambda = [];
    nt = nstep;
    for step = nt:-1:1
        fprintf('Solving reverse mode step %d of %d\n', nt - step + 1, nt);
        [~, lambda] = model.solveAdjoint(linsolve, getState, ...
                                         getObjective, schedule, lambda, step);
        dFdth = partialWRTparam(modelParam, getState, schedule, step);
        for kp = 1:numel(sens)
            sens{kp} = sens{kp} + dFdth{kp}'*lambda;
        end
    end
    
    % set ADI back to double
    opnms = {'T', 'pv'};
    for k  =1:numel(opnms)
        model.operators.(opnms{k}) = double(model.operators.(opnms{k}));
    end
    
   
    % add regularization term
    if ~isempty(opt.Regularization)
        assert(numel(sens)==numel(opt.Regularization), ...
               'Expected one regularization matrix for each sensitivity-vector');
        assert(~isempty(opt.modelBase), ...
               'Base-model needed to compute regularization' )
        sens = addRegularizationSens(sens, param, model, opt.modelBase, opt.Regularization);
    end
end

function dFdth = partialWRTparam(model, getState, schedule, step)
    eqs = blabla
end

function state = getStateFromInput(schedule, states, state0, i)
    if i == 0
        state = state0;
    elseif i > numel(schedule.step.val)
        state = [];
    else
        state = states{i};
    end
end

function model = initModelParametersADI(model, param, pTypes)
np = numel(param);
m = cell(1, np);
for k = 1:np
    switch param{k}
        case 'transmissibility'
            if strcmp(pTypes{k}, 'multiplyer')
                m{k} = ones(size(model.operators.T));
            else
                m{k} = model.operators.T;
            end
        case 'porevolume'
            if strcmp(pTypes{k}, 'multiplyer')
                m{k} = ones(size(model.operators.pv));
            else
                m{k} = model.operators.pv;
            end
        otherwise
            error(['Unknown parameter: ', param{k}]);
    end
end
[m{:}] = initVariablesADI(m{:});
for k = 1:np
    switch param{k}
        case 'transmissibility'
            if strcmp(pTypes{k}, 'multiplyer')
                model.operators.T = m{k}.*model.operators.T;
            else
                model.operators.T = m{k};
            end
        case 'porevolume'
            if strcmp(pTypes{k}, 'multiplyer')
                model.operators.pv = m{k}.*model.operators.pv;
            else
                model.operators.pv = m{k};
            end
    end
end
end

function sens = addRegularizationSens(sens, param, pTypes, model, modelBase, reg)
np = numel(param);
for k = 1:np
    switch param{k}
        case 'transmissibility'
            dparam = model.operators.T - modelBase.operators.T; 
            if strcmp(pTypes{k}, 'multiplyer')
                dparam = dparam./model.operators.T;
            end
        case 'porevolume'
            dparam = model.operators.pv - modelBase.operators.pv; 
            if strcmp(pTypes{k}, 'multiplyer')
                dparam = dparam./model.operators.pv;
            end
    end
    if any(~isfinite(dparam))
        dparam(~isfinite(dparam)) = 0;
    end
    
    sens{k} = sens{k} + reg{k}*dparam;
end
end
