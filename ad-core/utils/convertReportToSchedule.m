function [schedule, timesteps] = convertReportToSchedule(report, schedule)
% Create a new schedule based on actual ministeps from a simulation report
%
% SYNOPSIS:
%   schedule = convertReportToSchedule(report, schedule)
%
% DESCRIPTION:
%   Running a simulation schedule with a given set of control steps may
%   lead to several extra timesteps due to time step cutting/adjustments
%   done by the nonlinear solver. This utility converts the report output
%   from 'simulateScheduleAD' along with the schedule used into a new
%   schedule that accounts for ministeps actually taken.
%
%   Note that 'simulateScheduleAD' MUST be called with the option
%   'OutputMinisteps' set to true for this to do anything, otherwise it
%   will just output the report steps.
%
% REQUIRED PARAMETERS:
%   report   - Report from 'simulateScheduleAD' with 'OutputMinisteps' set
%              to true.
%
%   schedule - The schedule used to produce the report.
%
% RETURNS:
%   schedule - New schedule modified so that the ministeps are accounted
%              for.
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

    [timesteps, controls] = deal([]);
    if isstruct(report)
        crs = report.ControlstepReports;
    else
        crs = report;
    end
    if iscell(crs)
        n = numel(crs);
    else
        n = crs.numelData();
    end
    for i = 1:n
        cr = crs{i};
        steprep = cr.StepReports(cellfun(@(x) x.Converged, cr.StepReports));
        timesteps = [timesteps; cellfun(@(x) x.Timestep, steprep)]; %#ok
        controls = [controls;...
            repmat(schedule.step.control(i), numel(steprep), 1)]; %#ok
    end
    schedule.step.val = timesteps;
    schedule.step.control = controls;
end
