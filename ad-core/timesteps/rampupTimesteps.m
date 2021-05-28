function dT = rampupTimesteps(time, dt, n)
% Create timesteps that ramp up geometrically
%
% SYNOPSIS:
%   dT = rampupTimesteps(1*year, 30*day)
%   dT = rampupTimesteps(1*year, 30*day, 5)
%
% DESCRIPTION:
%   This function generates a timestep sequence for a given total time
%   interval that increases geometrically until it reaches some target
%   timestep. The rest of the interval is then divided into a number of
%   target timesteps.
%
% REQUIRED PARAMETERS:
%   time   - The total simulation time so that sum(dt) = time
%
%   dt     - Target timestep after initial ramp-up
%
%   n      - (OPTIONAL) Number of rampup steps. Defaults to 8.
%
% RETURNS:
%   dt     - Array of timesteps.
%
% NOTE:
%   The final timestep may be shorter than dt in order to exactly reach T.
%

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

    if nargin < 3
        n = 8;
    end
    if time == 0
        dT = [];
        return
    end
    % Initial geometric series
    dt_init = (dt./2.^[n n:-1:1])';
    cs_time = cumsum(dt_init);
    if any(cs_time > time)
        dt_init = dt_init(cs_time < time);
    end
    
    % Remaining time that must be discretized
    dt_left = time - sum(dt_init);
    % Even steps
    dt_rem = repmat(dt, floor(dt_left/dt), 1);
    % Final ministep if present
    dt_final = time - sum(dt_init) - sum(dt_rem);
    % Less than to account for rounding errors leading to a very small
    % negative time-step.
    if dt_final <= 0
        dt_final = [];
    end
    % Combined timesteps
    dT = [dt_init; dt_rem; dt_final];
end

