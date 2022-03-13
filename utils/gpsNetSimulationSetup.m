function setup = gpsNetSimulationSetup(gpsNet,iSchedule)
% Sets up a complete simulation problem for a GPSNEt model
%
% SYNOPSIS:
%     setup = gpsNetSimulationSetup(gpsNet, schedule);
%
% DESCRIPTION:
%
% The function takes GPSNet model and fine-scale simulation schedule as
% input, converts this schedule so that it can be used on the GPSNet, and
% generates a simulation setup structure for the GPSNet that holds the
% simulation model, the converted schedule, and initial data. This setup
% can be passed on to the evaluateMatch function from the optimization
% module.
%
% REQUIRED PARAMETERS:
%   gpsNet   - An instance of the GPSNet class that provides the
%              description of the GPSNet model
%
%   schedule - A simulation schedule defined on the fine grid whose
%              behavior the GPSNet tries to mimic.
%
% RETURNS:
%   setup    - A simulation setup for the GPSNet model corresponding to the
%              given simulation schedule.
% 
% SEE ALSO:
%   GPSNetModel, evaluateMatch

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

schedule = iSchedule;
for n=1:numel(schedule.control)
    for i=1:numel(gpsNet.W)
        Wi = gpsNet.W(i);
        Wi.type   = schedule.control(n).W(i).type;
        Wi.val    = schedule.control(n).W(i).val;
        Wi.WI     = sum(schedule.control(n).W(i).WI);
        Wi.status = schedule.control(n).W(i).status;
        for fn = fieldnames(Wi)'
            schedule.control(n).W(i).(fn{1}) = Wi.(fn{1});
        end
    end
end
setup = struct('model', gpsNet.model, ...
               'schedule', schedule, ...
               'state0', gpsNet.state0);
end
