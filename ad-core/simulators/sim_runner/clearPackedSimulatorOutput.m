function clearPackedSimulatorOutput(problems, varargin)
%Remove stored data for one or more packed simulation problem
%
% SYNOPSIS:
%   clearPackedSimulatorOutput(problems)
%   clearPackedSimulatorOutput(problem)
%   clearPackedSimulatorOutput(problems, 'Prompt', false);
%
% REQUIRED PARAMETERS:
%   problems - Either a single packed problem or multiple as a cell array.
%
% OPTIONAL PARAMETERS:
%
%   Prompt  - Should we prompt before deleting results? Default: true.
%
%   Default - Default choice for prompt. Either 'y' or 'n'. Default: 'n'.
%
%   start   - Delete all entries from step 'start' onwards. Default: 1.
%
%   stop    - Delete all entries up to and including 'stop. Default is the
%             all data present above start.
% RETURNS:
%   Nothing.
%
% EXAMPLE:
%   demoPackedProblems
%
% SEE ALSO:
%   packSimulationProblem

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

    opt = struct('Prompt', true,...
                 'Default', 'n', ...
                 'start', 1, ...
                 'stop', inf);
    opt = merge_options(opt, varargin{:});
    if isstruct(problems)
        problems = {problems};
    end
    for i = 1:numel(problems)
        problem = problems{i};
        s = problem.OutputHandlers.states;
        r = problem.OutputHandlers.reports;
        w = problem.OutputHandlers.wellSols;
        ns = s.numelData();
        nr = r.numelData();
        nw = w.numelData();
        if ns == 0 && nr == 0 && nw == 0
            continue
        end
        range = opt.start:min(opt.stop, ns);
        if opt.Prompt
            prompt = sprintf(['Do you want to delete %d states and ',...
                              'reports \nfor %s [%s]? y/n [%s]: '], ...
                              numel(range), problem.BaseName, problem.Name, opt.Default);
            str = input(prompt,'s');
            if ~strcmpi(str, 'y') && ~(strcmpi(opt.Default, 'y') && isempty(str))
                fprintf('Ok, will not remove files.\n');
                continue
            end
        end
        doPrint = mrstVerbose || opt.Prompt;
        dispif(doPrint, 'Removing files...');
        if opt.start == 1 && numel(range) == ns
            s.resetData();
            r.resetData();
            w.resetData();
        else
            s.resetData(range);
            r.resetData(range(range <= nr));
            w.resetData(range(range <= nw));
        end
        dispif(doPrint, ' Files removed.\n');
    end
end