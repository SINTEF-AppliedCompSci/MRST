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
    
    param  = checkParam(opt.Parameters);
    pTypes = opt.ParameterTypes;
    sens = cell(1, numel(param));
    [sens{:}] = deal(0);
  
    
    % inititialize parameters to ADI
    modelParam = initModelParametersADI(model, schedule, param);
    
    nstep = numel(schedule.step.val);
    lambda = [];
    nt = nstep;
    for step = nt:-1:1
        fprintf('Solving reverse mode step %d of %d\n', nt - step + 1, nt);
        [lami, lambda]= model.solveAdjoint(linsolve, getState, ...
                                         getObjective, schedule, lambda, step);
        eqdth = partialWRTparam(modelParam, getState, schedule, step);
        for kp = 1:numel(sens)
            for nl = 1:numel(lami)
                if isa(eqdth{nl}, 'ADI')
                    sens{kp} = sens{kp} - eqdth{nl}.jac{kp}'*lami{nl};
                end
            end
        end
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

function [model, schedule] = initModelParametersADI(model, schedule, param)
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
            case 'welltrans'
                % assume wells are the same throughout sim
                if stage == 1
                    m{k} = vertcat(schedule.control(1).W.WI);
                else
                    schedule = updateWellTrans(schedule, m{k});
                end
            case {'swcr', 'swu', 'sowcr'}  
                if stage == 1
                    m{k} = model.fluid.(param{k});
                else
                    model.fluid.(param{k}) = m{k};
                end
        end
    end
    % after setting values, initialize ADI:
    [m{:}] = initVariablesADI(m{:});
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

function ti = perm2directionalTrans(model, p, cdir)
% special utility function (typically for sensitivity calculations)
% for calculating transmissibility along coordinate direction cdir
% In particular:
% trans = t1+t2+t2, where ti = perm2dirtrans(model, perm(:, i), i);
assert(size(p,2)==1, ...
       'Input p should be single column representing permeability in direction cdir');
   
% make perm represent diag perm tensor with zero perm orthogonal to cdir

dp = double(p);
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

function param_proc = checkParam(param)
param_proc = param;
pix = 0;
for k = 1:numel(param)
    p = param{k};
    pix = pix+1;
    switch p
        case {'transmissibility', 'porevolume'}
            param_proc{pix} = param{k};
        case 'permeability'
            param_proc(pix:(pix+2)) = {'permx', 'permy', 'permz'};
            pix = pix+2;
        otherwise
            error('Unknown parameter: %s\n', param{k});
    end
end
end

