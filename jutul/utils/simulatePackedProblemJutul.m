function [ws, states] = simulatePackedProblemJutul(problem, varargin)
% Variant of simulatePackedProblem that uses Jutul as the backend.
%
% SEE ALSO:
%   simulateScheduleJutul

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

    s = problem.SimulatorSetup;
    state0 = s.state0;
    schedule = s.schedule;
    model = s.model;
    name = [problem.BaseName, '_', problem.Name];
    name = name(~isspace(name));
    [ws, states] = simulateScheduleJutul(state0, model, schedule, 'name', name, varargin{:});
end
