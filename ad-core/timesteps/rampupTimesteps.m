function dT = rampupTimesteps(time, dt, n)
    % time: total time
    % dt  : desired timesteps
    % n   : number of rampup steps
    if nargin < 3
        n = 8;
    end
    rampup = 2*dt;
    assert(time > rampup, 'Rampup time is larger than total time!');
    % Initial geometric series
    dt_init = (dt./2.^(n:-1:1))';
    % Remaining time that must be discretized
    dt_left = time - sum(dt_init);
    % Even steps
    dt_rem = repmat(dt, floor(dt_left/dt), 1);
    % Final ministep if present
    dt_final = time - sum(dt_init) - sum(dt_rem);
    if dt_final == 0
        dt_final = [];
    end
    % Combined timesteps
    dT = [dt_init; dt_rem; dt_final];
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
