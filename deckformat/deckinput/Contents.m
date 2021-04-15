% DECKINPUT
%
% Files
%   applyOperator          - Apply ECLIPSE/FrontSim operator to input array.
%   boxKeyword             - Uniform handling of BOX keyword
%   convertDeckUnits       - Convert ECLIPSE/FrontSim input deck units to MRST conventions
%   endboxKeyword          - Uniform handling of ENDBOX keyword
%   initializeDeck         - Initialise simple structure of ECLIPSE format type
%   matchWells             - Match wells to keyword records
%   readDefaultedKW        - Read data, possibly containing default designators, for a single keyword.
%   readDefaultedRecord    - Read data, possibly containing default designators, for a single record.
%   readEclipseDeck        - Read Simplified ECLIPSE input deck.
%   readEDIT               - Read edit
%   readGRID               - deck = readGRID(fid, dirname, deck)
%   readGridBoxArray       - Input grid array (global cell based) corresponding to current input box
%   readImmisciblePVTTable - Input immiscible ECLIPSE/FrontSim PVT table (e.g., PVDO, PVDG).
%   readMisciblePVTTable   - Input miscible ECLIPSE/FrontSim PVT table (e.g., PVTO, PVTG).
%   readPROPS              - deck = readPROPS(fid, dirname, deck)
%   readREGIONS            - Read regions
%   readRelPermTable       - Input ECLIPSE/FrontSim saturation function table (SWOF/SGOF &c).
%   readRUNSPEC            - Read runspec
%   readSCHEDULE           - Read schedule
%   readSOLUTION           - Read solution
%   readSUMMARY            - Read summary
%   readVFPINJ             - Read VFP tables for injector
%   readVFPPROD            - Read VFP tables for producer
%   readWellKW             - Read well definitions from an ECLIPSE Deck specification.
%   replaceKeywords        - Replace keywords. Internal function

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
