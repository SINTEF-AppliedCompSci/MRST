function timings = getReportTimings(report, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
    opt = struct('skipConverged', false, 'total', false);
    opt = merge_options(opt, varargin{:});
    if isstruct(report)
        report = report.ControlstepReports;
    end
    nstep = numel(report);
    timings = repmat(makeTiming(nan, nan, nan, nan, nan, nan), nstep, 1);
    for stepNo = 1:nstep
        rep = report{stepNo};
        assembly = 0;
        lsolve = 0;
        lsolve_prep = 0;
        total = rep.WallTime;
        n_assemblies = 0;
        for i = 1:numel(rep.StepReports)
            nasm = numel(rep.StepReports{i}.NonlinearReport);
            n_assemblies = n_assemblies + nasm;
            for j = 1:nasm
                r = rep.StepReports{i}.NonlinearReport{j};
                conv = r.Converged;
                if opt.skipConverged && conv
                    % We are only measuring the time in each iteration,
                    % skipping the cost of each assembly that just
                    % confirmed convergence.
                    continue;
                end
                assembly = assembly + r.AssemblyTime;
                % Linear solver
                ls = r.LinearSolver;
                if ~conv
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
        timings(stepNo) = makeTiming(total, assembly, lsolve, lsolve_prep, rep.Iterations, n_assemblies);
    end
    if opt.total
        % Sum up over all time-steps.
        fn = fieldnames(timings);
        tmp = timings(1);
        for i = 1:numel(fn)
            f = fn{i};
            tmp.(f) = sum(vertcat(timings.(f)));
        end
        timings = tmp;
    end
end

function s = makeTiming(total, assembly, linear, linear_prep, its, asm)
    s = struct('Total', total, ...
               'Assembly', assembly, ...
               'LinearSolve', linear, ...
               'LinearSolvePrep', linear_prep, ...
               'NumberOfAssemblies',  asm, ...
               'Iterations', its);
end
