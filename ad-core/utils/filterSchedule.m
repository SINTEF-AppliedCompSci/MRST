function s = filterSchedule(s)
    % Filter unused controls from a schedule
    [subs, ~, renum] = unique(s.step.control);
    s.step.control = renum;
    s.control = s.control(subs);
end
