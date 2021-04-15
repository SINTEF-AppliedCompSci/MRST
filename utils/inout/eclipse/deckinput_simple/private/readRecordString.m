function data = readRecordString(fid)
%Read single record string data (unconverted).
%
% SYNOPSIS:
%   rstring = readRecordString(fid)
%
% PARAMETERS:
%   fid     - Valid file identifier as obtained from FOPEN.
%
% RETURNS:
%   rstring - Record string.  Unconverted, raw character string.
%             Corresponds to single record data, not including record
%             terminator (i.e., '/' character).
%
% SEE ALSO:
%   `fopen`, `fgetl`.

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


   done = false;
   data = [];

   comment = '(--|#).*$';

   while ~done,
      lin   = fgetl(fid);
      valid = ischar(lin);

      if ~valid,
         if feof(fid),
            done = true;
         else
            [msg, errnum] = ferror(fid);

            assert (errnum ~= 0, ...
                    'Internal error in ''%s'' while reading ''%s''.\n', ...
                    mfilename, fopen(fid));

            error(['Input failed in ''%s'' while reading ''%s''.\n', ...
                   'System reports: %s.'], ...
                  mfilename, fopen(fid), msg);
         end
      end

      if ~done,
         lin  = regexprep(lin, comment, '');  % Remove comments.
         done = any(lin == '/');              % Record complete?
         data = [data, ' ', lin];  %#ok       % Willfully ignore MLINT.
      end
   end

   p = find(data == '/') - 1;
   if ~isempty(p),
      data = data(1 : p);                  % Exclude terminator character.
   end
   data(data == '''') = '';                % Exclude any quote characters.
end
