function schedule = convertReportToSchedule(report, schedule)
    [timesteps, controls] = deal([]);
    for i = 1:numel(report.ControlstepReports)
        cr = report.ControlstepReports{i};
        steprep = cr.StepReports(cellfun(@(x) x.Converged, cr.StepReports));
        timesteps = [timesteps; cellfun(@(x) x.Timestep, steprep)]; %#ok
        controls = [controls;...
            repmat(schedule.step.control(i), numel(cr.StepReports), 1)]; %#ok
    end
    schedule.step.val = timesteps;
    schedule.step.control = controls;
end