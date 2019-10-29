function stats = getPressureTransportIterations(seqReport)
    % Get pressure, transport and outer loop iterations for sequential or
    % fully-implicit model
    if isa(seqReport, 'ResultHandler')
        cr = seqReport;
        nstep = cr.numelData();
    elseif iscell(seqReport)
        cr = seqReport;
        nstep = numel(cr);
    else
        cr = seqReport.ControlstepReports;
        nstep = numel(cr);
    end
    cuts = zeros(nstep, 3);
    its = zeros(nstep, 3);
    steps = zeros(nstep, 2);
    isSeq = false;
    for j = 1:nstep
        % Loop over all control time-steps
        local_report = cr{j};
        srep = local_report.StepReports;
        for sr = 1:numel(srep)
            % Loop over each mini step report
            nouter = numel(srep{sr}.NonlinearReport);
            cuts(j, 3) = cuts(j, 3) + ~srep{sr}.Converged;
            for nl = 1:nouter
                % Loop over outer iterations
                NL = srep{sr}.NonlinearReport{nl};
                its(j, 3) = its(j, 3) + double(~NL.Converged);
                isSeq = isSeq | isfield(NL, 'PressureSolver');
                if isfield(NL, 'PressureSolver')
                    pr = srep{sr}.NonlinearReport{nl}.PressureSolver;
                    tr = srep{sr}.NonlinearReport{nl}.TransportSolver;

                    its(j, 1) = its(j, 1) + getIts(pr.StepReports);
                    if ~isempty(tr)
                        its(j, 2) = its(j, 2) + getIts(tr.StepReports);
                    end
                    steps(j, 1) = steps(j, 1) + numel(pr.StepReports);
                    if ~isempty(tr)
                        steps(j, 2) = steps(j, 2) + numel(tr.StepReports);
                    end
                end
            end
            if isSeq
                % Sequential models do "one" outer iteration extra as the
                % step which is fully converged is considered just as a
                % residual check.
                its(j, 3) = its(j, 3) + 1;
            end
        end
    end
    stats = struct('pressure',  its(:, 1),...
                   'transport', its(:, 2), ...
                   'outer',     its(:, 3), ...
                   'steps',     steps, ...
                   'cuts',      cuts);
end


function its = getIts(steprep)
    % Sum up non-converged iterations
    its = sum(cellfun(@(x) x.Iterations*x.Converged, steprep));
end