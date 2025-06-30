function [misfitVals ,varargout] = evaluateMatchTwinModelsSummands(pvec, obj, setup_model1, setup_model2, parameters, states_ref, prior, varargin)
% Utility function (for optimization) that simulates two models with parameters 
% obtained from the vector 'pvec' (scaled parameters) and computes vector 
% of residuals/mismatch with respect to a reference output state 'states_ref'
%
% SYNOPSIS:
%   misfitVals = evaluateMatchSummands(p, obj, setup_model1, setup_model2, parameters, states_ref, ['pn', pv, ...]) 
%
% REQUIRED PARAMETERS:
%   pvec         - An array containing the parameters' values scaled in unit-interval [0 ,1]
%   obj          - Objective function that evaluates the residuals between states
%                  evaluated at parameters p for both models and the reference state states_ref
%   setup_model1 - Simulation setup structure for model 1 containing: state0, model, and schedule.
%   setup_model2 - Simulation setup structure for model 2 containing: state0, model, and schedule.
%   parameters   - cell-array of parameters of class ModelParameter
%   states_ref   - Physical model states corresponding to the reference.
%
% OPTIONAL PARAMETERS:
%   'objScaling'   - scaling value for the objective function obj/objScaling.
%   'NonlinearSolver'- Subclass of `NonLinearSolver` suitable for solving the
%                      nonlinear systems of the forward model.
%   'Verbose'         - Indicate if extra output is to be printed.
%
% RETURNS:
%   misfitVals       - Difference between states(p) from both models and states_ref
%   sensitivities    - Jacobian of misfitVals with respect to p
%   wellSols         - Well solution at each control step for both models
%   states           - State at each control step for both models
%
% SEE ALSO:
% `evalObjective`, `computeSensitivitiesAdjointAD`, `unitBoxBFGS`

opt = struct('Verbose',           mrstVerbose(),...
             'NonlinearSolver', [],...
             'objScaling',1, ...
             'enforceBounds',  true, ...
             'accumulateResiduals', []);

[opt, extra] = merge_options(opt, varargin{:});
accum = opt.accumulateResiduals;
if isempty(accum)
    accum = struct('wells', [], 'types', [], 'steps', []);
end
nparam = cellfun(@(x)x.nParam, parameters);
p_org = pvec;
if opt.enforceBounds
    pvec = max(0, min(1, pvec));
end
pvec = mat2cell(pvec, nparam, 1);

% Set parameter values for model 1
pval_model1 = cell(size(parameters));
setupNew_model1 = setup_model1;
setupNew_model1.model.FlowDiscretization = [];
setupNew_model1.model.FlowPropertyFunctions = [];
for k = 1:numel(parameters)
    pval_model1{k}  = parameters{k}.unscale(pvec{k});
    setupNew_model1 = parameters{k}.setParameter(setupNew_model1, pval_model1{k});
end
[wellSols_model1,states_model1] = simulateScheduleAD(setupNew_model1.state0, setupNew_model1.model, setupNew_model1.schedule,...
                                       'NonLinearSolver',opt.NonlinearSolver,...
                                       'Verbose',opt.Verbose, extra{:});

% Set parameter values for model 2
pval_model2 = cell(size(parameters));
setupNew_model2 = setup_model2;
setupNew_model2.model.FlowDiscretization = [];
setupNew_model2.model.FlowPropertyFunctions = [];
for k = 1:numel(parameters)
    pval_model2{k}  = parameters{k}.unscale(pvec{k});
    if strcmp(prior.name, parameters{k}.name)
        scaled_prior = parameters{k}.scale(prior.value);
        dist_prior{k} =  pvec{k}(prior.location) - scaled_prior;
    end
    setupNew_model2 = parameters{k}.setParameter(setupNew_model2, pval_model2{k});
end
[wellSols_model2,states_model2] = simulateScheduleAD(setupNew_model2.state0, setupNew_model2.model, setupNew_model2.schedule,...
                                       'NonLinearSolver',opt.NonlinearSolver,...
                                       'Verbose',opt.Verbose, extra{:});

% Compute objective for both models
misfitVals = obj(setupNew_model1.model, states_model1, setupNew_model1.schedule, ...
                 setupNew_model2.model, states_model2, setupNew_model2.schedule, states_ref, dist_prior,false, [], [], []);

% Aggregate residuals if necessary
if ~isempty(accum.steps)
    tmp = repmat({0}, [max(accum.steps), 1]);
    for k = 1:numel(misfitVals)
        if accum.steps(k) > 0
            tmp{accum.steps(k)} = tmp{accum.steps(k)} + misfitVals{k};
        end
    end
    misfitVals = tmp;
end
misfitVals = (vertcat(misfitVals{:})).^(1/2);

% If additional outputs are requested
if nargout > 1
    objh = @(tstep,model1,model2,state1,state2) obj(setupNew_model1.model, states_model1, setupNew_model1.schedule, ...
                                                    setupNew_model2.model, states_model2, setupNew_model2.schedule, states_ref, dist_prior,true, tstep, state1, state2);
    nms = applyFunction(@(x)x.name, parameters);
    scaledGradient = cell(numel(nms), 1);

    setupNew_model1.model = setupNew_model1.model.validateModel();
    setupNew_model2.model = setupNew_model2.model.validateModel();
    
    gradient = computeSensitivitiesTwinModelsAdjointAD(setupNew_model1, states_model1, setupNew_model2, states_model2, parameters, objh, ...
                                                       'accumulateResiduals', accum, 'isScalar', false);
    for k = 1:numel(nms)
        scaledGradient{k} = parameters{k}.scaleGradient(gradient.(nms{k}), pval_model1{k});
    end
    J = vertcat(scaledGradient{:})' / opt.objScaling;
    
    nzIx  = abs(misfitVals) > eps*norm(misfitVals);
    numnz = nnz(nzIx);
    J(nzIx,:) = spdiags((2*misfitVals(nzIx)).^(-1), 0, numnz, numnz) * J(nzIx,:);
    
    varargout{1} = J;
    
    if nargout > 2
        [varargout{2:5}] = deal(wellSols_model1, states_model1, wellSols_model2, states_model2);
    end
end
end
    