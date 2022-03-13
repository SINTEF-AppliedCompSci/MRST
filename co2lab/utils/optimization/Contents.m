% OPTIMIZATION
%
% Files
%   computeToptraps    - Compute subscale trapping potential for all cells in a top surface grid,
%   optimizeControls   - Compute an optimal set of well controls ('rate' or 'bhp') for a proposed
%   optimizeFormation  - Optimize a injection scenario for a formation from the CO2-atlas.  Well
%   optimizeRatesIPOPT - Compute an optimal set of injection rates for a proposed
%   setSchedule        - Construct a schedule object that is convenient for use with

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
