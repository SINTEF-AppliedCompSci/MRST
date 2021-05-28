function [all_ws, all_states, all_reports, names, T, grids, act] = getMultiplePackedSimulatorOutputs(problems, varargin)
%Short description
%
% SYNOPSIS:
%   [all_ws, all_states, all_reports, names, T, grids, act] = ...
%   getMultiplePackedSimulatorOutputs(problems)
%
% REQUIRED PARAMETERS:
%   problems - Cell-array of packed simulation problems
%
% OPTIONAL PARAMETERS:
%   minSteps - Minimum number of steps that must be complete in order for
%              outputs to be read for a individual problem.
%
% RETURNS:
%
%   all_ws      - Cell array of wellSols for each problem with more than
%                 minSteps completed report steps.
%
%   all_states  - Cell array of states (see all_ws)
%
%   all_reports - Cell array of reports (see all_ws)
%
%   names       - Cell array of names for the simulations with more than
%                 minSteps completed report steps.
%
%   T           - Cell array of the times after each time-step for each
%                 problem.
%
%   grids       - Cell array of simulation grids.
%
%   act         - Boolean indicators for each problem, with true indicating
%                 that the output from the problem is present in the above
%                 cell arrays.
% EXAMPLE:
%   demoPackedProblems
%
% SEE ALSO:
%   getPackedSimulatorOutput, packSimulationProblem, simulatePackedProblem

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

    opt = struct('minSteps', 0);
    [opt, arg] = merge_options(opt, varargin{:});
    np = numel(problems);
    counts = zeros(np, 1);
    [T, all_states, all_ws, all_reports, names, grids] = deal(cell(np, 1));
    for i = 1:np
        problem = problems{i};
        ministeps = problem.SimulatorSetup.OutputMinisteps;
        [all_ws{i}, all_states{i}, all_reports{i}] = getPackedSimulatorOutput(problem, arg{:});
        names{i} = problem.Name;
        counts(i) = numelData(all_ws{i});
        schedule_time = cumsum(problem.SimulatorSetup.schedule.step.val);
        if ministeps
            [~, dt] = convertReportToSchedule(all_reports{i}, problem.SimulatorSetup.schedule);
            time = cumsum(dt);
        else
            time = schedule_time;
        end
        T{i} = time(1:counts(i));
        grids{i} = problem.SimulatorSetup.model.G;
        % Then, modify counts since they are interpreted as control steps
        if ministeps
            counts(i) = find(time <= T{i}(end), 1, 'last');
        end
    end
    act = counts >= opt.minSteps;
    all_ws = all_ws(act);
    all_states = all_states(act);
    all_reports = all_reports(act);
    names = names(act);
    T = T(act);
    grids = grids(act);
end