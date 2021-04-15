% INITIALIZATION
%
% Files
%   getDeckOPMData                        - Undocumented Utility Function
%   getInitializationRegionsBase          - Undocumented Utility Function
%   getInitializationRegionsBlackOil      - Undocumented Utility Function
%   getInitializationRegionsCompositional - Undocumented Utility Function
%   getInitializationRegionsDeck          - Undocumented Utility Function
%   getMinimumPhaseSaturations            - Undocumented Utility Function
%   getMinMaxPhaseSaturations             - Undocumented Utility Function
%   getMinMaxPhaseSaturationsFromRelPerm  - Undocumented Utility Function
%   initEclipseProblemAD                  - Set up all inputs to simulateScheduleAD from deck-file
%   initializeEquilibriumPressures        - Undocumented Utility Function
%   initializeEquilibriumSaturations      - Undocumented Utility Function
%   initStateBlackOilAD                   - Undocumented Utility Function
%   initStateDeck                         - Undocumented Utility Function
%   plotRegionContacts                    - Undocumented Utility Function

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
