function [val, der, wellSols, states] = evalObjective(u, obj, state0, model, schedule_org, scaling, varargin)
% Objective (and gradient) evaluation function based on input control vector u
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
opt=struct('Gradient','adjoint','pertub',1e-3);
opt = merge_options(opt,varargin{:});
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
schedule = control2schedule(u, schedule_org, scaling);

% run simulation:
[wellSols, states,report] = simulateScheduleAD(state0, model, schedule);
timesteps = getReportMinisteps(report);
assert(all(timesteps == schedule.step.val));
%% should we fix the time stepping ??
% schedule.step = schedule_new.step

% compute objective:
vals = obj(model, states, schedule);
val  = sum(cell2mat(vals))/objScaling;

% run adjoint:
if nargout > 1
    switch opt.Gradient
        case 'adjoint'
            objh = @(tstep,model, state) obj(model, states, schedule, 'ComputePartials', true, 'tStep', tstep,'state', state);
            g    = computeGradientAdjointAD(state0, states, model, schedule, objh);
            % scale gradient:
            der = scaleGradient(g, schedule, boxLims, objScaling);
            der = vertcat(der{:});
        case 'numerical'
            u_org=u;
            val_pert = nan(size(u));
            dp=nan(size(u));
            for i=1:numel(u)
                u_pert=u_org;                
                if(abs(u_pert(i))>0)
                    dp(i) = u_pert(i)*opt.pertub;
                else
                    dp(i) = opt.pertub;
                end
                u_pert(i) = u_pert(i) + dp(i);
                val_pert(i) = evalObjective(u_pert, obj, state0, model, schedule_org, scaling)
                der = (val_pert-val)./dp;
            end
            
        otherwise
            error('Gradient method %s not implemented',opt.Gradient);
    end
end
end

function grd = scaleGradient(grd, schedule, boxLims, objScaling)
dBox   = boxLims(:,2) - boxLims(:,1);
for k = 1:numel(schedule.control)
    grd{k} = (dBox/objScaling).*grd{k};
end
end
    
