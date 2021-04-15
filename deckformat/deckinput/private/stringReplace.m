function s = stringReplace(s, s1, s2)
%Replace string within another
%
% SYNOPSIS:
%   s = stringReplace(s1, s2, s3)
%
% DESCRIPTION:
%   Replaces all occurrences of the string 's2' in string 's1' with the
%   string 's3'.  The new string is returned.
%
% PARAMETERS:
%   s1 - Original input string.  May be a cell array of strings.
%
%   s2 - The string to replace.
%
%   s3 - Replacement string.
%
% RETURNS:
%   s - New string.  Is a cell array of strings if 's1' is a cell array of
%       strings.
%
% SEE ALSO:
%   `strrep`.

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
% STRREP is already sufficiently capable.

   try
      s = strrep(s, s1, s2);
   catch  %#ok
      % Octave compatibility.
      % Octave's STRREP does not accept a cell array of strings.
      %
      assert (iscell(s));
      s = cellfun(@(s0) strrep(s0, s1, s2), s, ...
                  'UniformOutput', false);
   end
end
