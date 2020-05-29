function [val, der, wellSols, states] = evalObjective(u, obj, state0, model, schedule, scaling)
% Objective (and gradient) evaluation function based on input control vector u
%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
minu = min(u);
maxu = max(u);
if or(minu < -eps , maxu > 1+eps)
    warning('Controls are expected to lie in [0 1]')
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
    
