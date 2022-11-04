function stepinfo = readInfoStep(fname)
%Read OPM Flow's INFOSTEP File
%
% SYNOPSIS:
%   stepinfo  = readInfoStep(fname)
%
% PARAMETERS:
%   fname - File name of INFOSTEP file.  Full or partial path.  Passed to
%           function FOPEN using mode 'rt'.  File columns assumed to be
%           separated by whitespace, meaning this function will break if
%           the unit strings are separated from the column names too.
%
% RETURNS:
%   stepinfo - Structure with fields inferred from the column headers.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
   if fid < 0
      error('Open:Failure', ...
            'Unable to Open INFOSTEP File ''%s'': %s', fname, msg);
   end

   clean = onCleanup(@() fclose(fid));

   vars = textscan(fgetl(fid), '%s');
   vars = regexprep(vars{1}, '\([^)]+\)', '');

   vals = textscan(fid, '%f');
   vals = num2cell(reshape(vals{1}, numel(vars), []) .', 1);

   stepinfo = cell2struct(vals, vars, 2);
end
