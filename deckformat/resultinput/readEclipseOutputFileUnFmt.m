function varargout = readEclipseOutputFileUnFmt(fname)
%Read unformatted (binary) ECLIPSE output/result file
%
% SYNOPSIS:
%   output = readEclipseOutputFileUnFmt(filename)
%
% PARAMETERS:
%   filename - Name (string) of file containing unformatted ECLIPSE results
%              (typically restart data, grid specification or metadata on a
%              particular result set).
%
%              The file represented by 'filename' will be opened using
%              function 'fopen' in mode ('r','ieee-be').
%
%              Use function 'readEclipseOutputFileFmt' to read formatted
%              (ASCII) ECLIPSE-compatible results.
%
% RETURNS:
%   output   - Data structure containing all field data represented in the
%              input stream 'filename'.  One structure field for each data
%              field in the file.
%
% SEE ALSO:
%   `readEclipseOutputFileFmt`, `fopen`.

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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


   [fid, msg] = fopen(fname, 'rb', 'ieee-be');
   if fid < 0, error([fname, ': ', msg]); end
   
   % skip first 4  
   fseek(fid, 4, 'cof');
   
   [varargout{1:nargout}] = readEclipseOutputStream(fid, @readFieldUnFmt);

   fclose(fid);
end
