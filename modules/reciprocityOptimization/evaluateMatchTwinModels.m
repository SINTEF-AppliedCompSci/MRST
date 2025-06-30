function [misfitVal,varargout] = evaluateMatchTwinModels(pvec, obj, setup_model1,setup_model2, parameters, states_ref,prior, varargin)
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
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

% Create new setup, and set parameter values for model1
pval_model1 = cell(size(parameters));
setupNew_model1 = setup_model1;
setupNew_model1.model.FlowDiscretization = [];
setupNew_model1.model.FlowPropertyFunctions = [];
for k = 1:numel(parameters)
    pval_model1{k}  = parameters{k}.unscale(pvec{k});
    if strcmp(prior.name, parameters{k}.name)
        scaled_prior = parameters{k}.scale(prior.value);
        dist_prior{k} =  pvec{k}(prior.location) - scaled_prior;
    end
    setupNew_model1 = parameters{k}.setParameter(setupNew_model1, pval_model1{k});
end
[wellSols_model1,states_model1] = simulateScheduleAD(setupNew_model1.state0, setupNew_model1.model, setupNew_model1.schedule,...
                                       'NonLinearSolver',opt.NonlinearSolver,...
                                       'Verbose',opt.Verbose, extra{:});

% Create new setup, and set parameter values for model2
pval_model2 = cell(size(parameters));
setupNew_model2 = setup_model2;
setupNew_model2.model.FlowDiscretization = [];
setupNew_model2.model.FlowPropertyFunctions = [];
for k = 1:numel(parameters)
    pval_model2{k}  = parameters{k}.unscale(pvec{k});
    setupNew_model2 = parameters{k}.setParameter(setupNew_model2, pval_model2{k});
end
[wellSols_model2,states_model2] = simulateScheduleAD(setupNew_model2.state0, setupNew_model2.model, setupNew_model2.schedule,...
                                       'NonLinearSolver',opt.NonlinearSolver,...
                                       'Verbose',opt.Verbose, extra{:});

misfitVals = obj(setupNew_model1.model, states_model1, setupNew_model1.schedule, setupNew_model2.model, states_model2, setupNew_model2.schedule,states_ref, dist_prior, false, [],[],[]);
misfitVal  = - sum(vertcat(misfitVals{:}))/opt.objScaling ;

if nargout > 1
    objh = @(tstep,model1,model2,state1,state2) obj(setupNew_model1.model, states_model1, setupNew_model1.schedule, setupNew_model2.model, states_model2, setupNew_model2.schedule,states_ref, dist_prior, true, tstep, state1,state2);
    nms = applyFunction(@(x)x.name, parameters);
    scaledGradient_1 = cell(numel(nms), 1);
    
    switch opt.Gradient
        case 'none'
            if nargout > 2
                [varargout{2:5}] = deal(wellSols_model1, states_model1,wellSols_model2, states_model2);
            end
            return
        case 'AdjointAD'
            setupNew_model1.model = setupNew_model1.model.validateModel();
            setupNew_model2.model = setupNew_model2.model.validateModel();

            gradient_1 = computeSensitivitiesTwinModelsAdjointAD(setupNew_model1, states_model1, setupNew_model2, states_model2, parameters, objh,...
                                                        'LinearSolver',opt.AdjointLinearSolver);
            % do scaling of gradient
            for k = 1:numel(nms)
               scaledGradient_1{k} = -parameters{k}.scaleGradient( gradient_1.(nms{k}), pval_model1{k});
%                scaledGradient_1{k}  = smoothdata(scaledGradient_1{k},'movmean',SmoothingFactor=1.0);
            end
        case 'PerturbationADNUM'
            % do manual pertubuation of the defiend control variabels
            eps_pert = opt.PerturbationSize;            
            val=nan(size(pvec));
%             try  % Try parallel loop
%                 parfor i=1:numel(p_org)
%                     val(i) = evaluateMatchTwinModels(perturb(p_org,i,eps_pert),...
%                          obj,setup_model1,setup_model2,parameters, states_ref,...
%                         'Gradient', 'none',...
%                         'NonlinearSolver',opt.NonlinearSolver,...
%                         'objScaling',opt.objScaling,'enforceBounds',  false);
%                 end
%             catch % Try serial loop instead
                for i=1:1
                    vall(i) = evaluateMatchTwinModels(perturb(p_org,i,eps_pert),...
                         obj,setupNew_model1,setupNew_model2,parameters, states_ref,...
                        'Gradient', 'none',...
                         'NonlinearSolver',opt.NonlinearSolver,...
                         'objScaling',opt.objScaling,'enforceBounds',  false);
                end
% 
%                 for i=2:2
%                     vall(i) = evaluateMatchTwinModels(perturb(p_org,145,eps_pert),...
%                          obj,setupNew_model1,setupNew_model2,parameters, states_ref,...
%                         'Gradient', 'none',...
%                          'NonlinearSolver',opt.NonlinearSolver,...
%                          'objScaling',opt.objScaling,'enforceBounds',  false);
%                 end
% 
%                     
%                 for i=3:3
%                     vall(i) = evaluateMatchTwinModels(perturb(p_org,289,eps_pert),...
%                          obj,setupNew_model1,setupNew_model2,parameters, states_ref,...
%                         'Gradient', 'none',...
%                          'NonlinearSolver',opt.NonlinearSolver,...
%                          'objScaling',opt.objScaling,'enforceBounds',  false);
%                 end
% 
%                                
%                 for i=4:4
%                     vall(i) = evaluateMatchTwinModels(perturb(p_org,946,eps_pert),...
%                          obj,setupNew_model1,setupNew_model2,parameters, states_ref,...
%                         'Gradient', 'none',...
%                          'NonlinearSolver',opt.NonlinearSolver,...
%                          'objScaling',opt.objScaling,'enforceBounds',  false);
%                 end
% 
%                  
%                 for i=5:5
%                     vall(i) = evaluateMatchTwinModels(perturb(p_org,1108,eps_pert),...
%                          obj,setupNew_model1,setupNew_model2,parameters, states_ref,...
%                         'Gradient', 'none',...
%                          'NonlinearSolver',opt.NonlinearSolver,...
%                          'objScaling',opt.objScaling,'enforceBounds',  false);
%                 end
% %             end 
% val(1:144,1)=vall(1);
% val(144*3+1:144*4,1)=vall(1);
% val(144*6+1:945,1)=vall(1);
% 
% val(1+144*1:288,1)=vall(2);
% val(577:720,1)=vall(2);
% val(1027:1107,1)=vall(2);
% 
% val(144*2+1:144*3,1)=vall(3);
% val(144*5+1:144*6,1)=vall(3);
% val(1189:1269,1)=vall(3);
% 
% val(946:1026,1)=vall(4);
% val(1108:1188,1)=vall(5);
val(1:length(p_org))=vall(1);

            gradient_1= (val'-misfitVal)./eps_pert;
            
            scaledGradient_1 = mat2cell(gradient_1, nparam, 1);            
        otherwise
            error('Greadient method %s is not implemented',opt.Gradient)
    end

        varargout{1} = vertcat(scaledGradient_1{:})/opt.objScaling;
end
if nargout > 2
    [varargout{2:5}] = deal(wellSols_model1, states_model1,wellSols_model2, states_model2);
end
if nargout > 6
   varargout{6} =  setupNew_model1;
end
end


% Utility function to perturb the parameter array in coordinate i with
% eps_pert
function p_pert = perturb(p_org,i,eps_pert) 
    p_pert = p_org;
    p_pert(i) = p_pert(i) + eps_pert;
end


    