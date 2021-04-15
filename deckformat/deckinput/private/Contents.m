% Files
%   boxIndices.m             - Construct linear (cell) indices corrsponding to current input box
%   completePVTO.m           - Fill in missing u-sat data from scaled copies of the FVF and viscosity.
%   defaultBox.m             - Default input box accessor function.
%   getEclipseKeyword.m      - Extract the next keyword from an ECLIPSE input deck.
%   gridBox.m                - Input box accessor function.
%   readEclipseIncludeFile.m - Read an ECLIPSE INCLUDE file.
%   readRecordString.m       - Read single record string data (unconverted).
%   readVector.m             - Input vector of floating point numbers from GRDECL file.
%   splitString.m            - Split string on multiple whitespace into a cell array of strings.
%   stringReplace.m          - Replace string within another

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
