function [misfitVal,varargout] = evaluateMatch_simple(p, obj, state0_org ,model_org,schedule_org,objScaling,parameters, states_ref)

nparam = cellfun(@(x)x.n, parameters);
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
[wellSols,states] = simulateScheduleAD(SimulatorSetup.state0, SimulatorSetup.model, SimulatorSetup.schedule);

misfitVals = obj(SimulatorSetup.model, states, SimulatorSetup.schedule, states_ref, false, [],[]);
misfitVal  = - sum(vertcat(misfitVals{:}))/objScaling ;

if nargout > 1
    objh = @(tstep,model,state) obj(model, states, SimulatorSetup.schedule, states_ref, true, tstep, state);
    gradient = computeSensitivitiesAdjointAD_V2(SimulatorSetup, states, parameters, objh);
    % do scaling of gradient
    nms = applyFunction(@(x)x.name, parameters);
    scaledGradient = cell(numel(nms), 1);
    for k = 1:numel(nms)
       scaledGradient{k} = parameters{k}.scaleGradient( gradient.(nms{k}), p{k});
    end
    varargout{1} = vertcat(scaledGradient{:})/objScaling;
end
if nargout > 2
    [varargout{2:3}] = deal(wellSols, states);
end
end
    