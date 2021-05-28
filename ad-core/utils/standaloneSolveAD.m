function [state, report] = standaloneSolveAD(state0, model, dt, varargin)
% Solve a single time-step with AD solvers for given forces
%
% SYNOPSIS:
%   [state, report] = standaloneSolveAD(state0, model, dt, 'W', W, 'maxIterations', 25)
%
%
% DESCRIPTION:
%   Stand-alone solver for AD. Useful for simple problems where a full
%   schedule is not required. Calls simulateScheduleAD internally.
%
% REQUIRED PARAMETERS:
%
%   state0    - Initial state.
%
%   model     - PhysicalModel instance for simulation.
%
%   dt        - Timestep length.
%
% OPTIONAL PARAMETERS:
%
%   'Various'         - Additional inputs are given to the following
%                       functions in order of priority: First, as possible
%                       driving forces, then as inputs to
%                       simulateScheduleAD.
%                       
% RETURNS:
%   state         - Updated state.
%
%   report        - Reports for simulator.
%
%
% SEE ALSO:
%   `simulateScheduleAD`

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

    schedule = simpleSchedule(dt);
    [schedule.control, extra] = merge_options(model.getValidDrivingForces(), varargin{:});
    
    [~, states, report] = ...
      simulateScheduleAD(state0, model, schedule, extra{:});
  
    state = states{1};
end