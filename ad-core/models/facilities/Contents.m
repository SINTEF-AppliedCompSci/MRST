% FACILITIES
%
% Files
%   combineMSwithRegularWells           - Combine regular and MS wells, accounting for missing fields
%   computeWellContributionsSingleWell  - Main internal function for computing well equations and source terms
%   GenericFacilityModel                - 
%   FacilityModel                       - A model coupling facilities and wells to the reservoir
%   MultisegmentWell                    - Derived class implementing multisegment wells
%   nozzleValve                         - Nozzle valve model
%   packPerforationProperties           - Extract variables corresponding to a specific well
%   selectFacilityFromDeck              - Pick FacilityModel from input deck
%   setupMSWellEquationSingleWell       - Setup well residual equations for multi-segmented wells - and only those.
%   setupWellControlEquationsSingleWell - Setup well controll (residual) equations 
%   SimpleWell                          - Base class implementing a single, instantaneous equilibrium well model
%   unpackPerforationProperties         - Unpack the properties extracted by packPerforationProperties. Internal function.
%   wellBoreFriction                    - Empricial model for well-bore friction
%   UniformFacilityModel                - Simplified facility model which is sometimes faster
%   reorderWellPerforationsByDepth      - Undocumented Utility Function

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

