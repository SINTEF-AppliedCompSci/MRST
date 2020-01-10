function timings = getReportTimings(report)
%Undocumented Utility Function

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
                    prep = 0;
                    if isfield(ls, 'PreparationTime')
                        % Fine grained statistics is available
                        prep = prep + ls.PreparationTime;
                        prep = prep + ls.PostProcessTime;
                        lsolve = lsolve + ls.LinearSolutionTime;
                    elseif isfield(ls, 'SolverTime')
                        lsolve = lsolve + ls.SolverTime;
                    end
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
