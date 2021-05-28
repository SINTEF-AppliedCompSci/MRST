function problem = initEclipsePackedProblemAD(deck, varargin)
% Set up a packed problem based on Eclipse input with reasonable defaults
%
% SYNOPSIS:
%   problem = initEclipsePackedProblemAD('/path/to/case.data')
%   problem = initEclipsePackedProblemAD(deck)
%
% REQUIRED PARAMETERS:
%   deck   - Either a deck struct or the valid path to a .DATA file.
%
% OPTIONAL PARAMETERS:
%   Name   - The name of the test case. Will use the name in the input file
%            if not provided.
%
%   Other  - Additional input arguments are passed onto
%            initEclipseProblemAD.
% RETURNS:
%   problem - Self-contained packed simulation problem which can be passed
%             to a number of routines, e.g:
%                 - simulatePackedProblem: Simulate and store output.
%                 - getPackedSimulatorOutput: Retrieve results.
%                 - simulatePackedProblemBackground: Simulate in seperate
%                 thread.
% EXAMPLE:
%   demoPackedProblems, simulateSPE1, simulateSPE9
%
% SEE ALSO:
%   getPackedSimulatorOutput, simulatePackedProblem, packSimulationProblem

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

    opt = struct('BaseName', [], 'Name', []);
    [opt, extra] = merge_options(opt, varargin{:});
    
    [state0, model, schedule, nls] = initEclipseProblemAD(deck, extra{:});
    name = opt.BaseName;
    if isempty(name)
        if isfield(model.inputdata.RUNSPEC, 'TITLE')
            name = model.inputdata.RUNSPEC.TITLE;
        else
            s = randi(26, 1, 20);
            letters = 'a':'z';
            name = letters(s);
        end
    end
    problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls, 'name', opt.Name);
end