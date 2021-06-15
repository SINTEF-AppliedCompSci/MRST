% Files
%   assignDENSITY     - dens of size ntpvtx3
%   assignGRAVITY     - Define Fluid Densites from API and Specific Fluid Gravity
%   assignEHYSTR      - Assigns flags for hysteresis computation, as well as input options for
%   assignMISC        - 
%   assignMSFN        - 
%   assignOILVIST     - Compute tables (static data)
%   assignPLMIXPAR    - 
%   assignPLYADS      - 
%   assignPLYMAX      - 
%   assignPLYROCK     - 
%   assignPLYSHEAR    - Polymer shear thinning/thickening
%   assignPLYSHLOG    - 
%   assignPLYVISC     - 
%   assignPMISC       - 
%   assignPVCDO       - 
%   assignPVDG        - 
%   assignPVDO        - 
%   assignPVDS        - 
%   assignPVTG        - 
%   assignPVTO        - 
%   assignPVTW        - 
%   assignROCK        - 
%   assignROCKTAB     - Undocumented Utility Function
%   assignRSCONSTT    - 
%   assignRelPerm     - 
%   assignRelPermScal - 
%   assignSDENSITY    - dens of size ntpvtx3
%   assignSGFN        - 
%   assignSGOF        - 
%   assignSGWFN       - Undocumented Utility Function
%   assignSHRATE      - 
%   assignSLGOF       - 
%   assignSOF2        - 
%   assignSOF3        - 
%   assignSPECHEAT    - Compute tables (static data)
%   assignSPECROCK    - Compute tables (static data)
%   assignSSFN        - 
%   assignSURFADS     - 
%   assignSURFCAPD    - Polymer shear thinning/thickening
%   assignSURFROCK    - 
%   assignSURFST      - 
%   assignSURFVISC    - 
%   assignSWFN        - 
%   assignSWOF        - 
%   assignTLMIXPAR    - 
%   assignVISCREF     - dens of size ntpvtx3
%   assignWATERVIST   - Compute tables (static data)
%   initDeckADIFluid  - Initialize AD-solver fluid from ECLIPSE-style input deck

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
