function schedule = combineSchedules(schedule, varargin)
% Combine multiple schedules to form a schedule with multiple controls

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
    opt = struct('makeConsistent', true);
    firstString = find(cellfun(@ischar, varargin), 1, 'first');
    if ~isempty(firstString)
        opt = merge_options(opt, varargin{firstString:end});
        varargin = varargin(1:firstString-1);
    end
    for i = 1:numel(varargin)
        schedule = combine(schedule, varargin{i});
    end
    if opt.makeConsistent
        schedule = makeScheduleConsistent(schedule);
    end
end

function s = combine(s, s2)
    offset = numel(s.control);
    s.control = [s.control; s2.control];
    s.step.val = [s.step.val; s2.step.val];
    s.step.control = [s.step.control; s2.step.control + offset];
end
