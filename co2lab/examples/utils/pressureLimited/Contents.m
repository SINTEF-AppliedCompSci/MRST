% PRESSURELIMITED
%
% Files
%   cellSubtract                     - we assume that a and b should be single-indexed
%   computeOverburdenPressure        - Compute overburden pressure acting on top of a formation, which is due to
%   findMaxPercentagePlimitReached   - Determine maximum percentage of pressure limit (or overburden pressure)
%   leak_penalizer_at_infinity_Rerun - Purpose: to compute the obj value at each time-step:
%   leakAtInfinity                   - computes mass (in Gt) leaked by time infinity as long as penalty is set
%   leakPenalizerAtInfinity          - new objective function to penalize leakage 'felt' at infinity, using
%   makeBjarmelandModel              - Simple Bjarmeland model for testing pressure-limited injection
%   maxPressureViolation             - states.pressure is a cell array.
%   postProcessExample               - Post-processing of Bjarmeland pressure-limited example
%   pressureAtCells                  - states.pressure is a cell array.
%   pressurePenalizer                - Penalize pressure of specific or all cells.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
