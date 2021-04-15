function data = removeQuotes(data)
%Remove explicit quote characters from string data
%
% SYNOPSIS:
%   data = removeQuotes(data)
%
% PARAMETERS:
%   data - String or cell array of strings that possibly contains quote (')
%          characters.
%
% RETURNS:
%   data - String or cell array of strings from which explicit quote
%          charaters are removed.  If the input is a cell array of strings,
%          then the individual strings are also subject to white-space
%          trimming using function 'strtrim'.
%
% SEE ALSO:
%   `strtrim`, `private/readRecordString`.

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

   if ischar(data),

      data = dequote(data);

   elseif iscellstr(data),

      data = cellfun(@dequote, data, 'UniformOutput', false);

   else

      error('Type:Error', 'Type ''%s'' Not Supported in ''%''', ...
            class(data), mfilename);

   end
end

%--------------------------------------------------------------------------

function s = dequote(s)
   s(s == '''') = '';

   s = strtrim(s);
end
