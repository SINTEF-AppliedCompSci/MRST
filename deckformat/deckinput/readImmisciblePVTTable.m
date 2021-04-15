function T = readImmisciblePVTTable(fid, ntab, ncol)
%Input immiscible ECLIPSE/FrontSim PVT table (e.g., PVDO, PVDG).
%
% SYNOPSIS:
%   Table = readImmisciblePVTTable(fid, ntab, ncol)
%
% PARAMETERS:
%   fid  - Valid file identifier (as obtained using FOPEN) to file
%          containing the PVT table.  FTELL(fid) is assumed to be at the
%          start of the PVT table (i.e., after the keyword, if any, which
%          prompted the reading of this table).
%
%   ntab - Number of PVT tables to input.  Typically corresponds to the
%          ECLIPSE/FrontSim parameter 'NTPVT' (item two of 'TABDIMS').
%
%   ncol - Number of columns in each of the 'ntab' tables.
%
% RETURNS:
%   Table - An ntab-by-1 cell array of immiscible PVT tables.  An
%           immiscible PVT table is an m-by-ncol array of type DOUBLE.  It
%           is the caller's responsibility to correctly interpret the
%           individual columns of a specific table.
%
% NOTE:
%   Immiscible PVT tables in ECLIPSE or FrontSim contain an arbitrary
%   number of records, each record separated by whitespace (typically
%   newline, '\n').  A table is terminated by a slash character ('/').
%
% SEE ALSO:
%   `readRelPermTable`, `readMisciblePVTTable`.

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


   T = cell([ntab, 1]);

   for i = 1 : ntab,
      % An immiscible PVT table is textually represented as a sequence of
      % rows, each row containing 'ncol' numbers.  The sequence is
      % terminated by a single slash ('/') character.
      %
      % Read data into a character string, split on whitespace, convert to
      % DOUBLE.
      %
      data = readRecordString(fid);    % Treat table data as single record.
      data = splitString(strtrim(data));

      T{i} = cellfun(@(s) sscanf(s, '%f'), reshape(data, ncol, []) .');
   end
end
