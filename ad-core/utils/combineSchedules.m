function schedule = combineSchedules(schedule, varargin)
% Combine multiple schedules to form a schedule with multiple controls
    for i = 1:numel(varargin)
        schedule = combine(schedule, varargin{i});
    end
end

function s = combine(s, s2)
    offset = numel(s.control);
    s.control = [s.control; s2.control];
    s.step.val = [s.step.val; s2.step.val];
    s.step.control = [s.step.control; s2.step.control + offset];
end
