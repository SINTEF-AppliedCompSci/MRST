function [val, der, wellSols, states] = evalObjective(u, obj, state0, model, schedule, scaling)
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

% function schedule = control2schedule(u, schedule, boxLims)
% nc = numel(schedule.control);
% nw = numel(schedule.control(1).W);
% [umin, umax] = deal(boxLims(:,1), boxLims(:,2));
% c  = 0;
% for cs = 1:nc
%     for w = 1:nw
%         c = c+1;
%         schedule.control(cs).W(w).val = u(c) *(umax(w)-umin(w))+umin(w);
%     end
% end
% end

function grd = scaleGradient(grd, schedule, boxLims, objScaling)
%dtc = accumarray(schedule.step.control, schedule.step.val);
%tScale = sum(dtc)./dtc;
dBox   = boxLims(:,2) - boxLims(:,1);
for k = 1:numel(schedule.control)
%    grd{k} = (tScale(k)*dBox/objScaling).*grd{k};
    grd{k} = (dBox/objScaling).*grd{k};
end
end
    
