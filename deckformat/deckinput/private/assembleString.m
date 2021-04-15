function s = assembleString(tokens)
%Assemble list of string tokens into single string
%
% SYNOPSIS:
%   s = assembleString(tokens)
%
% PARAMETERS:
%   tokens - Cell array of strings or return value from function
%            splitQuotedString.
%
% RETURNS:
%   s - Single character array derived from append each string token
%       sequentially from front to back.
%
% NOTE:
%   If the input is a cell array of strings, then the output is created by
%   separating each token by a single space character.  Otherwise, the
%   output is created by altering the unquouted tokens with the quoted
%   tokens.  No separator is inserted into the output string.
%
% SEE ALSO:
%   `private/splitQuotedString`.

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

   assert (iscellstr(tokens) || ...
           (isstruct(tokens) && ...
            all(isfield(tokens, { 'quoted', 'unquoted' }))), ...
          ['Input must be cell array of strings or return value ', ...
           'from function ''splitQuotedString''']);

   row     = @(a) reshape(a, 1, []);
   blank   = @(n) repmat({ ' ' }, [ 1, n ]);
   explode = @(s) [ s{:} ];

   if iscellstr(tokens),
      % Insert single blank between each token and at end of string.

      s = [ row(tokens) ; blank(numel(tokens)) ];

   elseif numel(tokens.unquoted) == 1,

      s = row(tokens.unquoted);

   else

      s = [ row(tokens.unquoted)           ; ...
            [ row(tokens.quoted), { '' } ] ];
   end

   s = explode(row(s));
end
