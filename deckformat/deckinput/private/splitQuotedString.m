function S = splitQuotedString(s)
%Separate a string with quoted substrings into constituent parts
%
% SYNOPSIS:
%   S = splitQuotedString(s)
%
% PARAMETERS:
%   s - Character array, possibly containing quoted substrings.
%
% RETURNS:
%   S - Scalar structure with the following fields
%         quoted   - Cell array of strings, each element of which
%                    corresponds to one quoted substring from 's'.
%                    Elements are ordered from left to right in order of
%                    appearance in 's'.
%
%         unquoted - Cell array of strings, each element of which
%                    corresponds to one unquoted substring from 's'.
%                    Elements are ordered from left to right in order of
%                    appearance in 's'.
%
% EXAMPLE:
%   s = '1.0 1* ''H--'' Z -1.234e5 ''*PROD 1'' ''MULTX-''';
%   S = splitQuotedString(s);
%
%   assert (numel(S.unquoted) == 4)
%   assert (isequal(S.unquoted(1 : (end - 1)), ...
%                   { '1.0 1* ', ' Z -1.234e5 ', ' '}))
%
%   assert (isequal(S.quoted, { '''H--''', '''*PROD 1''', '''MULTX-'''}))
%
% NOTE:
%   If the input string 's' does not contain any quoted substrings, then
%   the entire string is returned in a single element of '.unquoted' and
%   the 'quoted' cell array is empty.
%
%   Use function private/assembleString to reassemble the input string 's'
%   from the return value 'S'.
%
% SEE ALSO:
%   `private/assembleString`.

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

   s   = strtrim(s);
   pos = [ 0, strfind(s, ''''), numel(s) + 1 ];

   if numel(pos) > 2
      assert (mod(numel(pos), 2) == 0, ...
              'Non-terminated quoted string in input');

      % Separate input string into quoted and unquoted parts.  Order of
      % appearance within input string is
      %
      %    u q u q u ... u q u
      %
      % First and/or last unquoted string may be empty.

      % 1) Extract quoted substrings.
      qs = pos(2 : end - 1);
      q  = arrayfun(@(b, e) s(b : e), ... % Include quote char (')
                    qs(1 : 2 : end), ...
                    qs(2 : 2 : end), ...
                    'UniformOutput', false);

      % 2) Extract unquoted substrings.
      u = arrayfun(@(b, e) s((b + 1) : (e - 1)), ... % Exclude quote char
                   pos(1 : 2 : end), ...
                   pos(2 : 2 : end), 'UniformOutput', false);
   else
      % No quoted substring in 's'.  Treat entire string as unquoted.

      u = { s };
      q = { '' };

   end

   S = struct('quoted', { q }, 'unquoted', { u });
end
