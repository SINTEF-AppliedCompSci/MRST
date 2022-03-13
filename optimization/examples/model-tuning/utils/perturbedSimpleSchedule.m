function schedule = perturbedSimpleSchedule(dt, varargin)
%Undocumented Utitlity Function

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

opt = struct('pressureFac', 0.05, ...
             'rateFac',     0.05, ...
             'perturbStep', []);
[opt, extra] = merge_options(opt, varargin{:});
schedule = simpleSchedule(dt, extra{:});

pfac = opt.pressureFac;
rfac = opt.rateFac;
p   = opt.perturbStep;
if isempty(p)
    p = (1:numel(dt))';
end
np = max(p);
assert(numel(p) == numel(dt))

schedule.step.control = p;
schedule.control = repmat(schedule.control, [np, 1]);
W  = schedule.control(1).W; 
nw = numel(W);
isPressure = arrayfun(@(w)strcmp(w.type, 'bhp'), W);

for kw = 1:nw
    if isPressure(kw)
        fac = pfac;
    else
        fac = rfac;
    end
    for kp = 1:np
        schedule.control(kp).W(kw).val = W(kw).val*(1 + fac*(rand-.5));
    end
end
end
