% SIM_RUNNER
%
% Files
%   clearPackedSimulatorOutput        - Remove stored data for one or more packed simulation problem
%   copyPackedProblem                 - Make a copy of a packed problem with a new name and updated
%   getMultiplePackedSimulatorOutputs - Short description
%   getPackedSimulatorOutput          - Get output from a packed simulation problem
%   initEclipsePackedProblemAD        - Set up a packed problem based on Eclipse input with reasonable defaults
%   monitorBackgroundSimulations      - Monitor simulations running in the background or another session
%   PackedProblemManager              - Class for managing a set of packed simulation problems. See
%   packSimulationProblem             - Pack simulation inputs into a single atomic representation of problem
%   plotPackedProblem                 - Plot simulation results from a packed problem
%   simulatePackedProblem             - Simulate one or more packed simulation problems
%   simulatePackedProblemBackground   - Simulate a packed simulation problem as a seperate Matlab thread
%   simulatePackedProblemStandalone   - Stand-alone solver for running packed problems programmatically

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
