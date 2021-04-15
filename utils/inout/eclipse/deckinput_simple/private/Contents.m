% PRIVATE
%   Low-level support routines for simplified ECLIPSE data input
%
% Files
%   readEQUIL            - Read EQUIL keyword
%   readFaults           - Read FAULTS keyword
%   readMultfelt         - Read MULTFELT keyword
%   readOperator         - Read MULTFELT keyword
%   readOperatorModField - Read operator KEYWORDS
%   readRecordString     - Read single record string data (unconverted).
%   readSCHEDULE_orig    - Read (simplified version of) SCHEDULE section of ECLIPSE input deck.
%   readVector           - Input vector of floating point numbers from GRDECL file.
%   splitString          - Split string on multiple whitespace into a cell array of strings.

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
