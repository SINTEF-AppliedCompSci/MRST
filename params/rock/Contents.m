% ROCK
%   Routines related to rock data (permeability and porosity).
%
% Files
%   gaussianField         - Compute realization of lognormal, isotropic permeability field.
%   grdecl2Rock           - Extract rock properties from input deck (grdecl structure)
%   logNormLayers         - Compute realization of lognormal, isotropic permeability field.
%   makeRock              - Create rock structure from given permeabilty and porosity values
%   permeabilityConverter - Add tensor permeability to GRDECL struct
%   permTensor            - Expand permeability tensor to full format.
%   poreVolume            - Compute pore volumes of individual cells in grid.

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
