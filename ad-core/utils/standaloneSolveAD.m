function [state, report] = standaloneSolveAD(state0, model, dt, varargin)    
    schedule = simpleSchedule(dt);
    [schedule.control, extra] = merge_options(model.getValidDrivingForces(), varargin{:});
    
    [~, states, report] = ...
      simulateScheduleAD(state0, model, schedule, extra{:});
  
    state = states{1};
end