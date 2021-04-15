% SPE10
%
% Files
%   getSPE10_model_1_fluid   - Construct ADI Fluid Object for Model 1 of Tenth SPE CSP
%   getSPE10_model_1_relperm - Define Oil/Gas Relative Permeability Properties for Model 1 of tenth SPE CSP
%   getSPE10_model_1_rock    - Define rock properties for Model 1 of tenth SPE CSP
%   getSPE10rock             - Define rock properties for Model 2 of tenth SPE CSP
%   getSPE10setup            - Initialise properties for Model 2 of tenth SPE Comparative Solution Project
%   grdeclBox                - Make a GRDECL structure for simple corner-point grid, possibly faulted.
%   make_spe10_data          - Create on-disk (MAT file) representation of SPE 10 'rock' data.
%   setupSPE10_AD            - Undocumented Utility Function
%   SPE10_rock               - Define rock properties for Model 2 of tenth SPE CSP
%   SPE10_setup              - Initialise properties for Model 2 of tenth SPE Comparative Solution Project

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
