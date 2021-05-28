function timesteps = getReportMinisteps(report)
% Get the timesteps used for the ministeps of a report
%
% SYNOPSIS:
%   timesteps = getReportMinisteps(report)
%
% DESCRIPTION:
%   Get the actual ministeps used by simulateScheduleAD.
%
% REQUIRED PARAMETERS:
%   report - Report from simulateScheduleAD. 
%
% RETURNS:
%   timesteps - The timesteps used to solve the problem. These differ from
%               the control steps (schedule.step.val) in that they are the
%               actual timesteps taken by the solver (due to timestep
%               cutting, adaptive timestepping and so on)
%
% SEE ALSO:
%   SimulateScheduleAD

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
    timesteps = [];
    if iscell(report)
        report = struct('ControlstepReports', {report});
    end
    for i = 1:numel(report.ControlstepReports)
        cr = report.ControlstepReports{i};
        steprep = cr.StepReports(cellfun(@(x) x.Converged, cr.StepReports));
        timesteps = [timesteps; cellfun(@(x) x.Timestep, steprep)]; %#ok
    end
end


