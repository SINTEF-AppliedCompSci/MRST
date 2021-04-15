% Files
%   blackOilPVT.m            - Evaluate blackoil pvt functions and return phase properties.
%   eclipsePhaseProperties.m - Construct PVT evaluators for each active phase of an ECLIPSE/FrontSim deck
%   eclipseRelperm.m         - Construct relperm and capillary pressure evaluators from tabulated data.
%   processBlackOilDec.m     - Construct Black Oil PVT model evaluation function.
%   pvcdo.m                  - Evaluate dead oil properties for oil with constant compressibility.
%   pvdx.m                   - Interpolate dead oil/dead gas PVT properties.
%   pvtg.m                   - Interpolate live gas PVT properties.
%   pvto.m                   - Interpolate live oil PVT properties.
%   pvtw.m                   - Evaluate water PVT properties.

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
