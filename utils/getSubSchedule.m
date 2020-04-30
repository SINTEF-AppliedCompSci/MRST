function [subschedule, mappings] = getSubSchedule(schedule, mappings)
    % Get subschedule by extracting subforces for a subset of a full model
    subschedule = schedule;
    control     = subschedule.control;
    for i = 1:numel(control)
        [control(i), mappings] = getSubForces(control(i), mappings);
    end
    subschedule.control = control;
end