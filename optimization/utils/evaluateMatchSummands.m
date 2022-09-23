function [misfitVals ,varargout] = evaluateMatchSummands(pvec, obj, setup, parameters, states_ref, varargin)
% Utility function (for optimization) that simulates a model with parameters obtained 
% from the vector 'pvec' (scaled parameters) and computes vector of residuals/mismatch 
% with respect a reference output state 'states_ref'
%
% SYNOPSIS:
%   misfitVals = evaluateMatchSummands(p, obj, setup, parameters, states_ref, ['pn', pv, ...]) 
%
%
% REQUIRED PARAMETERS:
%   pvec         - An array containing the parameters' values scaled in unit-interval [0 ,1]
%   obj          - Objective function that evaluates the residuals between
%                  states evaluated at parameters p (states(p)) and the reference state states_ref
%   setup        - Simulation setup structure containing: state0, model, and schedule.
%   parameters   - cell-array of parameters of class ModelParameter
%   states_ref   - Physical model states corresponding to the reference.
%
% OPTIONAL PARAMETERS:
%   'objScaling'   - scaling value for the objective function obj/objScaling.
%   'NonlinearSolver'- Subclass of `NonLinearSolver` suitable for solving the
%                      non linear systems of the forward model.
%   'Verbose'         - Indicate if extra output is to be printed such as
%                       detailed convergence reports and so on.
% RETURNS:
%   misfitVals       - Diference between states(p) and states_ref
%   sesitivities     - Jacobian of misfitVals with respect p
%   wellSols         - Well solution at each control step (or timestep if
%                      'OutputMinisteps' is enabled.)
%   states           - State at each control step (or timestep if
%                      'OutputMinisteps' is enabled.)
%
% SEE ALSO:
% `evalObjective`, `computeSensitivitiesAdjointAD`, `unitBoxBFGS` 

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

% Create new setup, and set parameter values
pval = cell(size(parameters));
setupNew = setup;
setupNew.model.FlowDiscretization = [];
setupNew.model.FlowPropertyFunctions = [];
for k = 1:numel(parameters)
    pval{k}  = parameters{k}.unscale(pvec{k});
    setupNew = parameters{k}.setParameter(setupNew, pval{k});
end
[wellSols,states] = simulateScheduleAD(setupNew.state0, setupNew.model, setupNew.schedule,...
                                       'NonLinearSolver',opt.NonlinearSolver,...
                                       'Verbose',opt.Verbose, extra{:});
                                   
misfitVals = obj(setupNew.model, states, setupNew.schedule, states_ref, false, [],[]);
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

if nargout > 1
    objh = @(tstep,model,state) obj(setupNew.model, states, setupNew.schedule, states_ref, true, tstep, state);
    nms = applyFunction(@(x)x.name, parameters);
    scaledGradient = cell(numel(nms), 1);
    
    setupNew.model = setupNew.model.validateModel();
    gradient = computeSensitivitiesAdjointAD(setupNew, states, parameters, objh,...
                         'accumulateResiduals', accum, 'isScalar', false);
    for k = 1:numel(nms)
        scaledGradient{k} = parameters{k}.scaleGradient( gradient.(nms{k}), pval{k});
    end
    J = vertcat(scaledGradient{:})'/opt.objScaling;
    % adjust jacobian for sum r^2 -> sum r
    nzIx  = abs(misfitVals) > eps*norm(misfitVals);
    numnz = nnz(nzIx);
    J(nzIx,:) = spdiags((2*misfitVals(nzIx)).^(-1), 0, numnz, numnz)*J(nzIx,:);
    varargout{1} = J;
    if nargout > 2
        [varargout{2:3}] = deal(wellSols, states);
    end
    if nargout > 4
        varargout{4} =  setupNew;
    end
end
end

    