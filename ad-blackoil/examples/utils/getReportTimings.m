function timings = getReportTimings(report)
    if isstruct(report)
        report = report.ControlstepReports;
    end
    nstep = numel(report);
    timings = repmat(makeTiming(nan, nan, nan, nan, nan), nstep, 1);
    for stepNo = 1:nstep
        rep = report{stepNo};
        assembly = 0;
        lsolve = 0;
        lsolve_prep = 0;
        total = rep.WallTime;
        for i = 1:numel(rep.StepReports)
            for j = 1:numel(rep.StepReports{i}.NonlinearReport)
                r = rep.StepReports{i}.NonlinearReport{j};
                assembly = assembly + r.AssemblyTime;
                % Linear solver
                ls = r.LinearSolver;
                if not(r.Converged)
                    prep = ls.preparationTime + ls.postprocessTime;
                    lsolve = lsolve + ls.LinearSolutionTime;
                    lsolve_prep = lsolve_prep + prep;
                end
            end
        end
        timings(stepNo) = makeTiming(total, assembly, lsolve, lsolve_prep, rep.Iterations);
    end
end

function s = makeTiming(total, assembly, linear, linear_prep, its)
    s = struct('Total', total, ...
               'Assembly', assembly, ...
               'LinearSolve', linear, ...
               'LinearSolvePrep', linear_prep, ...
               'Iterations', its);
end