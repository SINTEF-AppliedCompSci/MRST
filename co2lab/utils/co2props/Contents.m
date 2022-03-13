% CO2PROPS
%
% Files
%   addSampledFluidProperties - Add density, viscosity or enthalpy properties to a fluid object.  The
%   boCO2                     - CO2 formation volume factor function
%   CO2CriticalPoint          - Values from paper by Span % Wagner
%   CO2props                  - Generate a set of CO2 property functions based on sampled data.
%   CO2VaporPressure          - Compute the CO2 vapor pressure for a given temperature
%   generatePropsTable        - Generates and saves a sampled table of fluid properties, using 'coolprops'.
%   propFilename              - Standardized filename generator for a sampled table of a given property of
%   SampledProp2D             - Create a structure with functions to interpolate a 2D sampled property and its

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
