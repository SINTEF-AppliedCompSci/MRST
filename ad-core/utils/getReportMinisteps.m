function timesteps = getReportMinisteps(report)
    timesteps = [];
    for i = 1:numel(report.ControlstepReports)
        cr = report.ControlstepReports{i};
        steprep = cr.StepReports(cellfun(@(x) x.Converged, cr.StepReports));
        timesteps = [timesteps; cellfun(@(x) x.Timestep, steprep)]; %#ok
    end
end