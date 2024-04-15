% CALIBRATESLEIPNER
%
% Files
%   addHeightData                     - Undocumented Utility Function
%   calibrateSleipnerSetup            - Setup the IEAGHG Sleipner model, including grid, entry rates, etc.
%   evaluatePlumeMatch                - Undocumented Utility Function
%   getrhoGmultfun                    - Undocumented Utility Function
%   getSleipnerOriginalInjectionRates - The following data corresponds to 32 years of CO2 entry into Sleipner
%   getSleipnerPlumeHeights           - Undocumented Utility Function
%   matchToCo2Surface                 - "surface.h" is permitted to contain nan's (i.e., missing data), in which
%   matchToPlumeData                  - Undocumented Utility Function
%   plotObsAndSim                     - Plot simulated and observed plume thicknesses
%   setIsotropicPermeabilityFun       - this is the 'setter' function for the ModelParameter handling
%   setrhoGmultfun                    - Undocumented Utility Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
