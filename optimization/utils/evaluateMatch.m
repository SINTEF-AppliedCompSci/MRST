function [misfitVal,varargout] = evaluateMatch(p, obj, state0_org ,model_org,schedule_org,objScaling,parameters, states_ref,varargin)
%Undocumented Utility Function

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
    'stepchange',false);

opt = merge_options(opt, varargin{:});


nparam = cellfun(@(x)x.n, parameters);
p_org=p;
p = max(0, min(1, p));
p = mat2cell(p, nparam, 1);

SimulatorSetup = struct('model', model_org, 'schedule', schedule_org, 'state0', state0_org);
initStateSensitivity = false;
ignore_these = zeros(numel(parameters),1);
for k = 1:numel(parameters)
    pval    = parameters{k}.unscale(p{k});
    SimulatorSetup = parameters{k}.setParameterValue(SimulatorSetup, pval);
    if (strcmp(parameters{k}.name,'initSw')||strcmp(parameters{k}.name,'p0'))
        initStateSensitivity = true;
        ignore_these(k)=1;
    end    
end
% skipping state0 (not yet part of ModelParameter)
%SimulatorSetup.model.toleranceCNV = 1e-6;
% figure(1), plot(SimulatorSetup.model.operators.pv)
% figure(2), plot(SimulatorSetup.model.operators.T)
% figure(3), plot(vertcat(SimulatorSetup.schedule.control.W.WI))
[wellSols,states] = simulateScheduleAD(SimulatorSetup.state0, SimulatorSetup.model, SimulatorSetup.schedule,...
    'NonLinearSolver',opt.NonlinearSolver);

misfitVals = obj(SimulatorSetup.model, states, SimulatorSetup.schedule, states_ref, false, [],[]);
misfitVal  = - sum(vertcat(misfitVals{:}))/objScaling ;

if nargout > 1
    objh = @(tstep,model,state) obj(model, states, SimulatorSetup.schedule, states_ref, true, tstep, state);
    nms = applyFunction(@(x)x.name, parameters);
    scaledGradient = cell(numel(nms), 1);
    
    switch opt.Gradient
        case 'none'
            if nargout > 2
                [varargout{2:3}] = deal(wellSols, states);
            end
            return
        case 'AdjointAD'
            if SimulatorSetup.model.toleranceCNV > 1e-6
                warning(['The accuracy in the gradient depend on the',...
                         ' acuracy on the CNV tolerance.',...
                         ' For good accuracy set  model.toleranceCNV <= 1e-6']);
            end
            gradient = computeSensitivitiesAdjointAD(SimulatorSetup, states, parameters, objh,...
                                                        'LinearSolver',opt.AdjointLinearSolver);            
            % do scaling of gradient
            for k = 1:numel(nms)
               scaledGradient{k} = parameters{k}.scaleGradient( gradient.(nms{k}), p{k});
            end
        case 'PerturbationADNUM'
            % do manual pertubuation of the defiend control variabels
            eps_scale = opt.PerturbationSize;            
            val=nan(size(p));
            try 
                parfor i=1:numel(p_org)
                    val(i) = evaluateMatch(perturb(p_org,i,eps_scale),...
                         obj,state0_org,model_org,schedule_org,objScaling,parameters, states_ref,...
                        'Gradient', 'none',...
                        'NonlinearSolver',opt.NonlinearSolver );
                end
            catch
                for i=1:numel(p_org)
                    val(i) = evaluateMatch(perturb(p_org,i,eps_scale),...
                         obj,state0_org,model_org,schedule_org,objScaling,parameters, states_ref,...
                        'Gradient', 'none',...
                        'NonlinearSolver',opt.NonlinearSolver );
                end
            end 
            gradient= (val-misfitVal)./eps_scale;
            
            gradient = mat2cell(gradient, nparam, 1);
            % do scaling of gradient
            for k = 1:numel(nms)
               scaledGradient{k} = parameters{k}.scaleGradient(gradient{k}, p{k});
            end
        otherwise
            error('Greadient method %s is not implemented',opt.Gradient)
    end

        varargout{1} = vertcat(scaledGradient{:})/objScaling;
end
if nargout > 2
    [varargout{2:3}] = deal(wellSols, states);
end
end

function p_pert = perturb(p_org,i,eps_pert) 
    p_pert = p_org;
    p_pert(i) = p_pert(i) + eps_pert;
end


    