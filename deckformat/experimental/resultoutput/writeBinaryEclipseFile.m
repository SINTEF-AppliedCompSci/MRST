function writeBinaryEclipseFile(filename, eclipse_struct)
%Write Eclipse data to binary files.
%
% SYNOPSIS:
%   writeBinaryEclipseFile(filename, eclipse_struct)
%
% PARAMETERS:
%   filename       - Filename of output file.
%   eclipse_struct - Data structure containing all field data of an Eclipse
%                    data file. One structure field for each data field.
%
% SEE ALSO:
%   `writeBinaryEclipseField`

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


fid = fopen(filename, 'wb');
fldnames = fieldnames(eclipse_struct);
% Write all fields
   for i = 1 : numel(fldnames),
      header.name  = fldnames{i};
      ecl_field = eclipse_struct.(header.name);
      header.data_type = ecl_field.type;
      header.number = numel(ecl_field.values);
      fid = writeBinaryEclipseField(ecl_field, header, fid);
   end
end
