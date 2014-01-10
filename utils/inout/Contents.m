% Low-level input routines for ECLIPSE deck.
%
% Files
%   addToTimeStruct           - Append single time step report structure to cumulative report structure
%   cacheDir                  - Define platform-dependent directory for temporary files.
%   convertUnitsOfReport      - Convert report object to standard metric units
%   eclipseNumberingToGrid    - convert eclipse output to rsol object
%   eclipseOutput2Sol         - Convert ECLIPSE output to MRST rSol object.
%   face2CellFace             - calculate cell faces of given faces
%   facesOfFault              - [faces] = facesOfFault(field,grdecl,G)
%   findEclipseNumberOfWells  - Enumerate ECLIPSE summary wells.
%   getEclipseKeyword         - Extract the next keyword from an ECLIPSE input deck.
%   grdecl2FaultFaces         -
%   grdecl2Rock               - Extract rock properties from input deck (grdecl structure)
%   grdecl2Saturation         - convert grdecl structure with water cut to saturation
%   permeabilityConverter     - add tensor permeability to the grdecl.permtensor
%   plotFaults                - plot faults of the grid
%   readCache                 - Read cached variables associated with given file into caller's workspace.
%   readDefaultedKW           - Read data, possibly containing default designators, for a single keyword.
%   readDefaultedRecord       - Read data, possibly containing default designators, for a single record.
%   readEQUIL                 - Read EQUIL keyword
%   readEclipseFile           - Input ECLIPSE or FrontSim Deck specification
%   readEclipseIncludeFile    - Read an ECLIPSE INCLUDE file.
%   readEclipseOutput         - Read single formatted ECLIPSE output/result file.
%   readEclipseOutputField    - Read next field of a formatted ECLIPSE output file
%   readEclipseTimeData       - Read formatted ECLIPSE summary file.
%   readEclipseResults        - Read formatted output/results from ECLIPSE run.
%   readFaults                - Read FAULTS keyword
%   readGRDECL                - Read Eclipse file assuming corner-point grid format.
%   readMultfelt              - Read MULTFELT keyword
%   readOperator              - Read MULTFELT keyword
%   readOperatorModField      - Read operator KEYWORDS
%   readPARAMS                - Read Eclipse or FrontSim data file.
%   readpvt                   - Read Black Oil PVT data from ECLIPSE parameter file into table.
%   readRelPermTable          - Input ECLIPSE/FrontSim relperm table (SWOF/SGOF &c).
%   readSCHEDULE              - Read (simplified version of) SCHEDULE section of ECLIPSE input deck.
%   readSelectedEclipseOutput - Read selected fields from single formatted ECLIPSE output/result file.
%   readVector                - Input vector of floating point numbers from GRDECL file.
%   readWellKW                - Read well definitions from an ECLIPSE Deck specification.
%   wellCalculateProduction   - Convert MRST well solution data to ECLIPSE keyword representation.
%   writeCache                - Write callers workspace to matfile in directory ./.cache/
%   writeGRDECL               - Write a GRDECL structure out to permanent file on disk.

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
