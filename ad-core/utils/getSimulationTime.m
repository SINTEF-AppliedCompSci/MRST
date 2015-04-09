function times = getSimulationTime(states, report)
    if numel(states) == numel(report)
        timesteps = cumsum(report.SimulationTime);
    else
        timesteps = [];
        for i = 1:numel(report.ControlstepReports)
            cr = report.ControlstepReports{i};
            steprep = cr.StepReports(cellfun(@(x) x.Converged, cr.StepReports));
            timesteps = [timesteps; cellfun(@(x) x.Timestep, steprep)]; %#ok
        end
    end
    times = cumsum(timesteps);
end
