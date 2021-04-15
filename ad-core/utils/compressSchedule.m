function scheduleCompressed = compressSchedule(schedule)
% Compress schedule to take the longest possible timesteps while honoring controls
%
% SYNOPSIS:
%   schedule = compressSchedule(schedule)
%
% DESCRIPTION:
%   Take a already defined schedule and combine all successive timesteps
%   with the same controls. This is useful to make a schedule suitable for
%   dynamic timestepping, as the control steps correspond only to the hard
%   limits set by changing well/bc controls.
%
% REQUIRED PARAMETERS:
%
%   scheduleDeck - A deck struct, typically from convertDeckScheduleToMRST
%                  or manually created. 
%
% RETURNS:
%   scheduleCompressed - Schedule ready for simulation in 'simulateScheduleAD'.

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

    finished = false;
    i = 1;
    stepno = 1;
    
    nstep = numel(schedule.step.val);
    newsteps = nan(nstep, 1);
    newctrls = nan(nstep, 1);
    oldctrl = schedule.step.control;
    oldvals = schedule.step.val;
    while ~finished
        % Find the next point in the schedule where the control changes
        ctrl = oldctrl(i);
        nextctrl = find(oldctrl ~= ctrl & ~isnan(oldctrl), 1, 'first');
        if isempty(nextctrl)
            nextctrl = nstep + 1;
            finished = true;
        end
        % Indices of concurrent controls that are the same at the current
        % position
        current = i:(nextctrl-1);
        
        % Combine into a large step
        newsteps(stepno) = sum(oldvals(current));
        newctrls(stepno) = ctrl;
        
        % Mask out old entries
        oldctrl(current) = nan;
        % Jump to next changed control, since we are done up to i-1.
        i = nextctrl;
        stepno = stepno + 1;
    end
    scheduleCompressed = schedule;
    strip = @(x) x(~isnan(x));
    scheduleCompressed.step.val     = strip(newsteps);
    scheduleCompressed.step.control = strip(newctrls);
end