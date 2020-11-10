function problem = getSubSchedule(problem, steps)
% Example function for optional QoI class property `processProblemFn`. This
% function extracts a subset of the schedule after setting the sample.
% Useful for e.g. when the schedule read from a deck.
    schedule = problem.SimulatorSetup.schedule;
    schedule.step.val     = schedule.step.val(steps);
    schedule.step.control = schedule.step.control(steps);
    problem.SimulatorSetup.schedule = schedule;
end