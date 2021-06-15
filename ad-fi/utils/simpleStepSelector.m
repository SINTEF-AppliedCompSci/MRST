function [dt, dt_history] = simpleStepSelector(dt_history, dt_current, its, varargin)
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

    opt = struct('targetIts', 5, ...
        'dt_min', 0.001*day, ...
        'dt_max', inf, ...
        'stepModifier', 1.5);

    opt = merge_options(opt, varargin{:});


    if its == opt.targetIts || its == 0
        modifier = 1;
    elseif its > opt.targetIts
        modifier = 1/opt.stepModifier;
    else
        modifier = opt.stepModifier;
    end

    dt = modifier*dt_current;

    dt_history = vertcat(dt_history, [dt, its]);
    dt = max(dt, opt.dt_min);
    dt = min(dt, opt.dt_max);
end
