% UTILS
%
% Files
%   grdeclSloping            - Make a GRDECL structure for simple corner-point grid, possibly faulted.
%   makeJohansenVEgrid       - Make a VE model based upon a data set of the Johansen formation
%   makeSleipnerVEmodel      - Make a VE model based upon the Sleipner data set from ieaghg.org
%   makeSlopingAquifer       - Make an example model of a sloping aquifer with heterogeneous rock props
%   makeSlopingAquiferBig    - Make an VE model based upon a data set obtained from the IGEMS project
%   readIGEMSIRAP            - Read and process IGEMS .irap grids.
%   readIrapClassicAsciiSurf - Read an .irap classic ASCII surface
%   sinusDeck                - Make a GRDECL structure for simple sloping sinus-formed reservoir.
%   sinusDeckAdiVE           - Make a GRDECL structure for simple sloping sinus-formed reservoir.

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
