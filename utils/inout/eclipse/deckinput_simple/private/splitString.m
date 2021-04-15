function a = splitString(s)
%Split string on multiple whitespace into a cell array of strings.
%
% SYNOPSIS:
%   a = splitString(s)
%
% PARAMETERS:
%   s - A string consisting of words separated by whitespace.
%
% RETURNS:
%   a - A cell array of strings, each element being a separate word from
%       the string 's'.
%
% SEE ALSO:
%   `regexp`.

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


% This function exists only as an Octave compatibility shim.  MATLAB's
% REGEXP is already sufficiently capable to handle this task.

   patt = '\s+';
   try
      a = regexp(s, patt, 'split');
   catch %#ok
      % Octave compatibility shim.
      %
      % STRSPLIT exists only in Octave.
      % It is an error to get here in a MATLAB run.
      %
      s = regexprep(s, patt, ' ');  % Multiple -> single whitespace char.
      a = strsplit (s, ' ', true);  % Split on single whitespace char.
   end
end
