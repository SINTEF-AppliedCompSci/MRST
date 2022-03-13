% UPSCALING
%
% Files
%   getCapVisDist               - Upscale two-way oil-water distribution
%   pcVsUpscaledSwGravityBinary - Create fine scale capillary pressure pc as a function of upscaled sW.
%   upAbsPerm                   - Undocumented Utility Function
%   upAbsPermAvg                - Power averaging of of the absolute permeability
%   upAbsPermPres               - Undocumented Utility Function
%   upFracFlowOW                - Upscale fractional flow curves.
%   upPcOW                      - Upscale capillary pressure curves.
%   upPolyAds                   - Upscale polymer adsorption isotherm using a simple average.
%   upPolyRk                    - Undocumented Utility Function
%   upPoro                      - Upscale porosity of the block by pore volume averaging
%   upRelPerm                   - Upscaling of relative permeability
%   upRelPermEPS                - Upscaling of relperm curves by end-point scaling (EPS)
%   upRelPermPV                 - Undocumented Utility Function
%   upscaleTransNew             - Calculate upscaled transmissibilities for a coarse model
%   upWells                     - Upscale wells

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
