function tokens = tokenizeRecord(data)
%Split record into individual tokens while preserving quoted strings
%
% SYNOPSIS:
%   tokens = tokenizeRecord(record)
%
% PARAMETERS:
%   record - Character string, possibly containing quoted substrings.
%
% RETURNS:
%   tokens - Cell array of strings containing non-empty tokens from the
%            input record.  Quoted substrings are preserved as individual
%            tokens.
%
% EXAMPLE:
%   record = '1.0 1* ''H--'' Z -1.234e5 ''*PROD 1'' ''MULTX-''';
%   tokens = tokenizeRecord(record);
%   expect = { '1.0', '1*', '''H--''', 'Z', '-1.234e5', ...
%              '''*PROD 1''', '''MULTX-''' };
%
%   assert (isequal(tokens, expect))
%
% SEE ALSO:
%   `private/splitString`, `private/removeQuotes`.

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

   assert (ischar(data), 'Input must be character string');

   S = splitQuotedString(data);

   if numel(S.unquoted) > 1,

      % Input data contains quoted strings that must be protected from
      % automatic splitting on white-space.  Split only the unquoted
      % strings and explode the resulting tokens into the output data.

      u  = cellfun(@splitString, S.unquoted, 'UniformOutput', false);

      tokens = [ reshape([ u(1:end-1) ; S.quoted ], 1, []), u{end} ];
      tokens = [ tokens{:} ];

      if ischar(tokens), tokens = { tokens }; end

      tokens = tokens(~cellfun(@isempty, tokens));

   else

      tokens = splitString(strtrim(data));

   end
end
