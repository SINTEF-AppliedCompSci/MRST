function schedule = simpleSchedule(dt, varargin)
%Make a schedule with varying timesteps and fixed wells/bc/src terms
%
% SYNOPSIS:
%   schedule = simpleSchedule(timesteps);
%   schedule = simpleSchedule(timesteps, 'W', W, 'src', src, 'bc', bc);
%
% PARAMETERS:
%   
%   dt    - Vector (column/row) of desired timesteps.
%
% OPTIONAL PARAMETERS:
%
%   W -  Wells to be used in the schedule. The wells will be active in
%        all timesteps.
%   
%   BC - Boundary conditions to be used in the schedule. The boundary
%        conditions will be active in all timesteps. 
%
%   src - Source terms to be used in the schedule. The sourceterms will be
%        active in all timesteps.
%
% RETURNS:
%   schedule - struct suitable for further modification, or for input to
%              `simulateScheduleAD`.
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

    opt = struct('W', [], ...
                 'src',   [], ...
                 'bc',    []);
    opt = merge_options(opt, varargin{:});
    
    sz = size(dt);
    assert(sz(1) == 1 || sz(2) == 1, ...
        'Timesteps must be given as a row or column vector');
    
    dt = reshape(dt, [], 1);
    nstep = numel(dt);
    
    schedule = struct();
    schedule.step.val = dt;
    schedule.step.control = ones(nstep, 1);
    
    % Make a single control out of the input driving forces
    schedule.control = opt;
end
