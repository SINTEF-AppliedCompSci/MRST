function writeBinaryEclipseRestartData(prefix, restart_data)
%Write unformatted Eclipse restart data to binary files. Multiple file output.
%
% SYNOPSIS:
%   writeBinaryEclipseRestartData(prefix, restart_data)
%
% PARAMETERS:
%   prefix  - Path-name prefix from which to construct restart file names.
%             Filenames = 'prefix'.Xnnnn, where first nnnn is 0000 and is
%             incremented by one for every new file.
%
%   restart_data - Struct array. One struct for each report step.
%
%
% SEE ALSO:
%   `writeBinaryEclipseFile`, `writeBinaryEclipseField`

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


% Write one file for every report step
  for i = 1 : numel(restart_data),
    filename = sprintf('%s.X%04d', prefix, i-1)
    writeBinaryEclipseFile(filename, restart_data{i});
  end
end
