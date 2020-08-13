function sens = computeSensitivitiesAdjointAD(state0, states, model, schedule, getObjective, varargin)
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
%   state0       - Physical model state at `t = 0`
%
%   states       - All previous states. Must support the syntax 
%                  `state = states{i}`. If the problem is too large to fit in
%                  memory, it can use `ResultHandler` class to retrieve files
%                  from the disk.
%                  
%   model        - Subclass of PhysicalModel class such as
%                  `ThreePhaseBlackOilModel` that models the physical 
%                  effects we want to study.
%
%   schedule     - Schedule suitable for `simulateScheduleAD`.
%
%   getObjective - Function handle for getting objective function value 
%                  for a given timestep with derivatives. Format: @(tstep)
%
% OPTIONAL PARAMETERS:
%
%   'Parameters'     - Defaults to `{'transmissibility', 'porevolumes'}`. 
% 
%   'ParameterTypes' - Defaults to 'value' for all. The other option is
%                      'multiplier'
%
%   'Regularization' - One (sparse) matrix for each parameter-set, requires
%                      also modelBase
%
%   'LinearSolver'   - Subclass of `LinearSolverAD` suitable for solving the
%                      adjoint systems.
%
% RETURNS:
%   sens - Sensitivites.
%
% SEE ALSO:
%   `computeGradientAD`

%{
Copyright 2009-2017 SINTEF Digital, Applied Mathematics.

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

    opt = struct('Parameters', {{'transmissibility', 'porevolume'}}, ...
                 'ParameterTypes', [], ...
                 'Regularization', [], ...
                 'LinearSolver',     [], ...
                 'initStateSensitivity', false, ...
                 'modelBase', []);
    opt = merge_options(opt, varargin{:});
    
    % check that OutputStateFunctions is empty, give warning if not 
    if ~isempty(model.OutputStateFunctions)
        warning('Model property ''OutputStateFunctions'' is non-empty, this may result in incorrect sensitivities.');
    end
    
    getState = @(i) getStateFromInput(schedule, states, state0, i);
    
    if isempty(opt.LinearSolver)
        linsolve = BackslashSolverAD();
    else
        linsolve = opt.LinearSolver;
    end
    
    [param, pTypes]  = checkParam(opt.Parameters, opt.ParameterTypes);
         
    sens = struct;
    for k = 1:numel(param)
        sens.(param{k}) = 0;
    end
    
    % inititialize parameters to ADI
    [modelParam, scheduleParam, paramValues] = initModelParametersADI(model, schedule, param);
    % reset discretization/flow functions to account for AD-parameters
    modelParam.FlowDiscretization = [];
    modelParam.FlowPropertyFunctions = [];
    modelParam = validateModel(modelParam);
    
    
    nstep = numel(schedule.step.val);
    lambda = [];
    nt = nstep;
    for step = nt:-1:1
        fprintf('Solving reverse mode step %d of %d\n', nt - step + 1, nt);
        [lami, lambda]= model.solveAdjoint(linsolve, getState, ...
                                         getObjective, schedule, lambda, step);
        eqdth = partialWRTparam(modelParam, getState, scheduleParam, step);
        for kp = 1:numel(param)
            for nl = 1:numel(lami)
                if isa(eqdth{nl}, 'ADI')
                    sens.(param{kp}) = sens.(param{kp}) + eqdth{nl}.jac{kp}'*lami{nl};
                end
            end
        end 
    end
    
    % compute partial derivative of first eq wrt init state and compute
    % initial state sensitivities
    if opt.initStateSensitivity
        forces = model.getDrivingForces(schedule.control(schedule.step.control(1)));
        forces = merge_options(model.getValidDrivingForces(), forces{:});
        model = model.validateModel(forces);
        
        state0 = model.validateState(state0);
        % set wellSols just to make subsequent function-calls happy, sensitivities wrt wellSols doesn't make sense anyway
        state0.wellSol = states{1}.wellSol;
        dt = schedule.step.val(1);
        
        problem = model.getAdjointEquations(state0, states{1}, dt, forces,...
            'iteration', inf, 'reverseMode', true);
        pnm = problem.primaryVariables;
        for kn = 1:numel(pnm)
            sens.init.(pnm{kn}) = 0;
            for nl = 1:numel(lami)
                if isa(problem.equations{nl}, 'ADI')
                    sens.init.(pnm{kn}) = sens.init.(pnm{kn}) + problem.equations{nl}.jac{kn}'*lami{nl};
                end
            end
        end
    end
                            
        
        
    
    % multiply by value if type is multiplier:
    multIx = find(strcmp(pTypes, 'multiplier'));
    for k = 1:numel(multIx)
        sens.(param{multIx(k)}) = sens.(param{multIx(k)}).*paramValues{multIx(k)};
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

function eqdth = partialWRTparam(model, getState, schedule, step)
    validforces = model.getValidDrivingForces();
    current = getState(step);
    before  = getState(step - 1);
    dt_steps = schedule.step.val;
    dt = dt_steps(step);
    lookupCtrl = @(step) schedule.control(schedule.step.control(step));
    % get forces and merge with valid forces
    forces = model.getDrivingForces(lookupCtrl(step));
    forces = merge_options(validforces, forces{:});
    model = model.validateModel(forces);
    % Initial state typically lacks wellSol-field, so add if needed
    if step == 1
        before = model.validateState(before);
    end
    
    problem = model.getEquations(before, current, dt, forces, 'iteration', inf, 'resOnly', true);
    eqdth   = problem.equations;
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

function [model, schedule, paramValues] = initModelParametersADI(model, schedule, param)
np = numel(param);
m = cell(1, np);
% in stage 1, set values, in stage 2, initialize ADI
for stage = 1:2
    for k = 1:np
        switch param{k}
            case 'transmissibility'
                if stage == 1
                    m{k} = model.operators.T; 
                else
                    model.operators.T = m{k};
                end
            case 'porevolume'
                if stage == 1
                    m{k} = model.operators.pv;
                else
                    model.operators.pv = m{k};
                end
            case 'permx'
                % always compute sensitivity for perm in X,Y and Z
                % -directions, input must be organized such that permy and
                % permz appear at index k+1 and k+2, resp.
                if stage == 1
                    nd = size(model.rock.perm, 2);
                    if nd == 1
                        ix = [1 1 1];
                    elseif nd ==3
                        ix = [1 2 3];
                    else
                        error('Sensitivity code assumes either isotropic or diagnoal (3-elem) perm');
                    end
                    for ii = 1:3
                        m{k+ii-1} = model.rock.perm(:,ix(ii));
                    end
                else
                    th = cell(1,3);
                    for ii = 1:3
                        % half -transmissibilities projected along each
                        % coordinate direction
                        th{ii} = perm2directionalTrans(model, m{k+ii-1}, ii);
                    end
                    cf = model.G.cells.faces(:,1);
                    nf = model.G.faces.num;
                    % mapping from from cell-face to face
                    M = sparse(cf, (1:numel(cf))', 1, nf, numel(cf));
                    % consider only internal faces
                    ie = model.operators.internalConn;
                    model.operators.T = 1./(M(ie,:)*(1./(th{1}+th{2}+th{3})));
                end
            case 'conntrans'
                % assume wells are the same throughout sim and
                if stage == 1
                    %m{k} = collectConnectionFactors(schedule.control);
                    m{k} = vertcat(schedule.control(end).W.WI);
                else
                    W    = schedule.control(1).W;
                    ncon = arrayfun(@(x)numel(x.WI), W);
                    for step = 1:numel(schedule.control)
                        ix = 0;
                        for wn = 1:numel(W) 
                            schedule.control(step).W(wn).WI = m{k}(ix + (1:ncon(wn))');
                            ix = ix + ncon(wn);
                        end
                    end
                end
            case {'swl', 'swcr', 'swu', 'sowcr'}  
                if stage == 1
                    m{k} = model.fluid.(param{k});
                else
                    % reset fluid-functions
                    model.fluid = initSimpleScaledADIFluid(model.fluid, param{k}, m{k});
                end
        end
    end
    if stage == 1
        paramValues = m;
        % after setting values, initialize ADI:
        [m{:}] = initVariablesADI(m{:});
    end
end
if isnumeric(model.operators.pv)
    model.operators.pv = double2ADI(model.operators.pv, m{1});
end
end


function sens = addRegularizationSens(sens, param, pTypes, model, modelBase, reg)
warning('Regularization not tested ...')
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

function ti = perm2directionalTrans(model, p, cdir)
% special utility function for calculating transmissibility along coordinate direction cdir
% In particular:
% trans = t1+t2+t2, where ti = perm2dirtrans(model, perm(:, i), i);
assert(size(p,2)==1, ...
       'Input p should be single column representing permeability in direction cdir');
% make perm represent diag perm tensor with zero perm orthogonal to cdir
dp = value(p);
r.perm = zeros(numel(dp), 3);
r.perm(:, cdir) = dp;
ti = computeTrans(model.G, r);
if isa(p, 'ADI')
    % make ti ADI (note ti is linear function of p)
    p = p./dp;
    cellno = rldecode(1 : model.G.cells.num, diff(model.G.cells.facePos), 2).';
    ti = ti.*p(cellno);
end
end

function [param_proc, types_proc] = checkParam(param, types)
if isempty(types)
   types = repmat({'value'}, 1, numel(param));
end
assert(numel(param)==numel(types));
param_proc = param;
types_proc = types;

pix = 0;
for k = 1:numel(param)
    p = param{k};
    pix = pix+1;
    switch p
        case {'transmissibility', 'porevolume', 'conntrans', 'swl', 'swcr', 'swu', 'sowcr'}
            param_proc{pix} = param{k};
            types_proc{pix} = types{k};
        case 'permeability'
            param_proc(pix:(pix+2)) = {'permx', 'permy', 'permz'};
            types_proc(pix:(pix+2)) = repmat(types(k), 1, 3);
            pix = pix+2;
        otherwise
            error('Unknown parameter: %s\n', param{k});
    end
end
if any(strcmp(param_proc, 'transmissibility')) && any(strcmp(param_proc, 'permx'))
    warning('Cannot compute sensitivities for transmissibility and permeability in the same run ... ')
    ix = strcmp(param_proc, 'transmissibility');
    param_proc = param_proc(~ix);
    types_proc = types_proc(~ix);
end 
end

