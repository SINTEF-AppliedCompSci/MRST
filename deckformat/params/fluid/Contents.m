% Models for three-phase compressible fluids.
%
% Files
%   CompressibleFluidBase                 - Base class for compressible fluids
%   ConstantCompressibilityFluids         - Class representing compressible fluids with constant compressibility.
%   Fluid                                 - Class representing a fluid composition as a collection of fluids.
%   IdealGases                            - Simple model of an ideal gas.
%   initBlackoilFluid                     - Return function that evaluate Black Oil PVT model
%   initBlackOilPVT                       - Construct Black Oil PVT model evaluation function.
%   initBlackOilRelPerm                   - Construct Black Oil relperm evaluation function.
%   initCompressibleFluid                 - Construct PVT and relperm functions, constant compr., single phase.
%   initSimpleThreephaseCompressibleFluid - Construct PVT and relperm functions, constant compressibility three-phase
%   initTrangensteinBellFluid             - Construct PVT and relperm functions, constant compressibility three-phase

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
