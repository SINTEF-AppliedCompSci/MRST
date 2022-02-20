% OPTIMIZATION
%
% Files
%   argmaxCubic                      - find max of cubic polynom through p1, p2
%   computeDefaultBasis              - Undocumented Utility Function
%   control2well                     - Update val-fields of W from controls u
%   evalObjectiveDiagnostics         - Undocumented Utility Function
%   getObjectiveDiagnostics          - Get an objective function for the diagnostics optimization routines
%   lineSearch                       - lineSearch -- helper function which performs line search based on
%   linsolveWithTimings              - Undocumented Utility Function
%   optimizeDiagnosticsBFGS          - Optimize well rates/bhps to minimize/maximize objective 
%   optimizeTOF                      - Optimize well rates based on diagnostics-based objective function
%   optimizeWellPlacementDiagnostics - Optimize the placement of wells using flow diagnostics
%   plotWellRates                    - Simple utility for plotting well rates. Used by optimizeTOF.
%   plotWellsPrint                   - Plots wells as simple colored circles. Helper for diagnostics examples.
%   solveStationaryPressure          - Solve incompressible, stationary pressure without gravity with optional TOF output
%   SolveTOFEqsADI                   - Solve the time of flight equations for solveStationaryPressure.
%   tofRobustFix                     - Presently does nothing. Can be overriden if needed.
%   unitBoxBFGS                      - Undocumented Utility Function
%   well2control                     - Produce control vector u from target-wells

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
