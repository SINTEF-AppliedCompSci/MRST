function [misfitVal,varargout] = evaluateMatch(pvec, obj, setup, parameters, states_ref, varargin)
% Utility function (for optimization) that simulates a model with parameters obtained 
% from the vector 'pvec' (scaled parameters) and computes mismatch with respect a 
% reference output state 'states_ref'
%
% SYNOPSIS:
%   misfitVal = evaluateMatch(p, obj, setup, parameters, states_ref, ['pn', pv, ...]) 
%  [misfitVal, sesitivities] = evaluateMatch(...) 
%  [misfitVal, sesitivities, wellSols, states] = evaluateMatch(...) 
%  [misfitVal, ~, wellSols, states] = evaluateMatch(...,'Gradient','none')
%
% DESCRIPTION:
%   For a given parameter array p, compute mistmach and sesitivities with regards to parameter p
%
% REQUIRED PARAMETERS:
%   pvec         - An array containing the parameters' values scaled in unit-interval [0 ,1]
%   obj          - Objective function that evaluates the diferences/match between
%                  states evaluated at parameters p (states(p)) and the reference state states_ref
%   setup        - Simulation setup structure containing: state0, model, and schedule.
%   parameters   - cell-array of parameters of class ModelParameter
%   states_ref   - Physical model states corresponding to the reference.
%
% OPTIONAL PARAMETERS:
%   'Gradient'       - Method to calculate the sensitivities/gradient:
%                      'AdjointAD':        Compute parameter sensitivities
%                                          using adjoint simulation  (default)
%                      'PerturbationADNUM': Compute parameter sensitivities
%                                           using perturbations (first-order forward finite diferences)
%                      'none':              Avoind computing parameters sensitivities
%   'PerturbationSize'- small value <<1 to perturb parameter p (default = 1e-7)
%   'objScaling'   - scaling value for the objective function obj/objScaling.
%   'AdjointLinearSolver' - Subclass of `LinearSolverAD` suitable for solving the
%                      adjoint linear systems.
%   'NonlinearSolver'- Subclass of `NonLinearSolver` suitable for solving the
%                      non linear systems of the forward model.
%   'Verbose'         - Indicate if extra output is to be printed such as
%                       detailed convergence reports and so on.
% RETURNS:
%   misfitVal        - Diference between states(p) and states_ref
%   sesitivities     - Gradient of misfitVal with respect p
%   wellSols         - Well solution at each control step (or timestep if
%                      'OutputMinisteps' is enabled.)
%   states           - State at each control step (or timestep if
%                      'OutputMinisteps' is enabled.)
%
% SEE ALSO:
% `evalObjective`, `computeSensitivitiesAdjointAD`, `unitBoxBFGS` 

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
             'Gradient',   'AdjointAD',...
             'NonlinearSolver', [],...
             'AdjointLinearSolver',[],...
             'PerturbationSize',1e-7,...
             'objScaling',1, ...
             'enforceBounds',  true);

[opt, extra] = merge_options(opt, varargin{:});

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
misfitVal  = - sum(vertcat(misfitVals{:}))/opt.objScaling ;

if nargout > 1
    objh = @(tstep,model,state) obj(setupNew.model, states, setupNew.schedule, states_ref, true, tstep, state);
    nms = applyFunction(@(x)x.name, parameters);
    scaledGradient = cell(numel(nms), 1);
    
    switch opt.Gradient
        case 'none'
            if nargout > 2
                [varargout{2:3}] = deal(wellSols, states);
            end
            return
        case 'AdjointAD'
            gradient = computeSensitivitiesAdjointAD(setupNew, states, parameters, objh,...
                                                        'LinearSolver',opt.AdjointLinearSolver);            
            % do scaling of gradient
            for k = 1:numel(nms)
               scaledGradient{k} = parameters{k}.scaleGradient( gradient.(nms{k}), pval{k});
            end
        case 'PerturbationADNUM'
            % do manual pertubuation of the defiend control variabels
            eps_pert = opt.PerturbationSize;            
            val=nan(size(pvec));
            try  % Try parallel loop
                parfor i=1:numel(p_org)
                    val(i) = evaluateMatch(perturb(p_org,i,eps_pert),...
                         obj,setup,parameters, states_ref,...
                        'Gradient', 'none',...
                        'NonlinearSolver',opt.NonlinearSolver,...
                        'objScaling',opt.objScaling,'enforceBounds',  false);
                end
            catch % Try serial loop instead
                for i=1:numel(p_org)
                    val(i) = evaluateMatch(perturb(p_org,i,eps_pert),...
                         obj,state0_org,model_org,schedule_org,parameters, states_ref,...
                        'Gradient', 'none',...
                        'NonlinearSolver',opt.NonlinearSolver,...
                        'objScaling',opt.objScaling, 'enforceBounds',  false);
                end
            end 
            gradient= (val-misfitVal)./eps_pert;
            
            scaledGradient = mat2cell(gradient, nparam, 1);            
        otherwise
            error('Greadient method %s is not implemented',opt.Gradient)
    end

        varargout{1} = vertcat(scaledGradient{:})/opt.objScaling;
end
if nargout > 2
    [varargout{2:3}] = deal(wellSols, states);
end
if nargout > 4
   varargout{4} =  setupNew;
end
end


% Utility function to perturb the parameter array in coordinate i with
% eps_pert
function p_pert = perturb(p_org,i,eps_pert) 
    p_pert = p_org;
    p_pert(i) = p_pert(i) + eps_pert;
end


    