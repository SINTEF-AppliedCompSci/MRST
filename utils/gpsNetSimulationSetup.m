function setup = gpsNetSimulationSetup(gpsNet,iSchedule)
% Sets up a complete simulation problem for a GPSNEt model
%
% SYNOPSIS:
%     setup = gpsNetSimulationSetup(gpsNet, schedule);
%
% DESCRIPTION:

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

schedule = iSchedule;
for n=1:numel(schedule.control)
    for i=1:numel(gpsNet.W)
        Wi = gpsNet.W(i);
        Wi.type   = schedule.control(n).W(i).type;
        Wi.val    = schedule.control(n).W(i).val;
        Wi.WI     = sum(schedule.control(n).W(i).WI);
        Wi.status = schedule.control(n).W(i).status;
        schedule.control(n).W(i) = Wi;
    end
end
setup = struct('model', gpsNet.model, ...
               'schedule', schedule, ...
               'state0', gpsNet.state0);
end

