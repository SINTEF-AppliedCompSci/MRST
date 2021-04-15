function output = readEclipseOutputFileFmt(fname)
%Read formatted (ASCII) ECLIPSE output/result file
%
% SYNOPSIS:
%   output = readEclipseOutputFileFmt(filename)
%
% PARAMETERS:
%   filename - Name (string) of text file containing formatted ECLIPSE
%              results (typically restart data, grid specification or
%              metadata on a particular result set).
%
%              The file represented by 'filename' will be opened using
%              function 'fopen' in mode 'rt'.
%
%              Use function 'readEclipseOutputFileUnFmt' to read
%              unformatted (binary) ECLIPSE-compatible results.
%
% RETURNS:
%   output   - Data structure containing all field data represented in the
%              input stream 'filename'.  One structure field for each data
%              field in the file.
%
% SEE ALSO:
%   `readEclipseOutputFileUnFmt`, `fopen`.

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


   [fid, msg] = fopen(fname, 'rt');
   if fid < 0, error([fname, ': ', msg]); end

   output = readEclipseOutputStream(fid, @readFieldFmt);

   fclose(fid);
end
