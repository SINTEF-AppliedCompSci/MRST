function [subschedule, mappings] = getSubSchedule(schedule, mappings)
    subschedule = schedule;
    control     = subschedule.control;
    for i = 1:numel(control)
        [control(i), mappings] = getSubForces(control(i), mappings);
    end
    subschedule.control = control;
end