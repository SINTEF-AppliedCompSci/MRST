function [val, der, wellSols, states] = evalObjective(u, obj, state0, model, schedule, scaling)
% Objective (and gradient) evaluation function based on input control vector u
minu = min(u);
maxu = max(u);
if or(minu < -eps , maxu > 1+eps)
    warning('Controls are expected to lie in [0 1]\n')
end

boxLims = scaling.boxLims;
if isfield(scaling, 'obj')
    objScaling = scaling.obj;
else
    objScaling = 1;
end

% update schedule:
schedule = control2schedule(u, schedule, scaling);

% run simulation:
[wellSols, states] = simulateScheduleAD(state0, model, schedule);

% compute objective:
vals = obj(wellSols, schedule);
val  = sum(cell2mat(vals))/objScaling;

% run adjoint:
if nargout > 1
    objh = @(tstep)obj(wellSols, schedule, 'ComputePartials', true, 'tStep', tstep);
    g    = computeGradientAdjointAD(state0, states, model, schedule, objh);
    % scale gradient:
    der = scaleGradient(g, schedule, boxLims, objScaling);
    der = vertcat(der{:});
end
end

function grd = scaleGradient(grd, schedule, boxLims, objScaling)
dBox   = boxLims(:,2) - boxLims(:,1);
for k = 1:numel(schedule.control)
    grd{k} = (dBox/objScaling).*grd{k};
end
end
    
