function times = getSimulationTime(states, report)
%Get the global time for a set of states produced by simulateScheduleAD
%
% SYNOPSIS:
%   times = getSimulationTime(states, report)
%
% DESCRIPTION:
%   Get the time for each state output by simulateScheduleAD.
%
% REQUIRED PARAMETERS:
%   states     - Cell array of states as given by simulateScheduleAD. Can
%                be either per control step or per ministep.
%
%   report     - Report given by simulateScheduleAD.
%
% RETURNS:
%   times      - Array with same dimensions as states, giving the values of
%                the time in the simulation model for each step. The
%                initial state passed to simulateScheduleAD is assumed to
%                be at time zero.
%
% SEE ALSO:
%   simulateScheduleAD

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

    if numel(states) == numel(report)
        timesteps = cumsum(report.SimulationTime);
    else
        timesteps = [];
        for i = 1:numel(report.ControlstepReports)
            cr = report.ControlstepReports{i};
            steprep = cr.StepReports(cellfun(@(x) x.Converged, cr.StepReports));
            timesteps = [timesteps; cellfun(@(x) x.Timestep, steprep)]; %#ok
        end
    end
    times = cumsum(timesteps);
end
