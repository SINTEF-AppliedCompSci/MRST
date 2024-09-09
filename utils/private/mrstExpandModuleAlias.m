function names = mrstExpandModuleAlias(names, lst)
%Expand Module Names Aliases to Canonical Names
%
% SYNOPSIS:
%   names = mrstExpandModuleAlias(names, lst)
%
% PARAMETERS:
%   names - Module names, possibly containing alises.  Cell array of
%           character vectors.
%
%   lst   - List of canonical module names.  Cell array of character
%           vectors.
%
% RETURNS:
%   names - Expanded list of names.  Aliases in 'names' replaced by their
%           corresponding canonical names from 'lst'.
%
% SEE ALSO:
%   mrstPath, mrstModule.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

   names = reshape(names, 1, []);
   aliases = get_module_aliases();

   [ia, ja] = identify_aliases(aliases, names);

   if ~ isempty(ia)
      % There is at least one alias name in 'names'.  Replace each alias
      % with their corresponding list of canonical names from 'lst'.
      names = expand_aliases(names, lst, aliases, ia, ja);
   end
end

%--------------------------------------------------------------------------

function aliases = get_module_aliases()
   % Aliases is an M-by-2 cell array of character vectors
   %
   % First column is the module name alias (simple string)
   %
   % Second column is a pattern matching the modules corresponding to that
   % alias.  This should normally be an anchored regular expression.
   aliases = { ...
      'co2lab', '^co2lab-.*$'; ...
   };
end

%--------------------------------------------------------------------------

function [ia, ja] = identify_aliases(aliases, names)
   % Aliases vary across rows, module names vary across columns.

   [ia, ja] = ...
      find(strcmp(repmat(reshape(aliases(:,1), [], 1), ...
                         [1, numel(names)]), ...
                  repmat(names, [size(aliases, 1), 1])));
end

%--------------------------------------------------------------------------

function names = expand_aliases(names, lst, aliases, ia, ja)
   % Canonical names vary across rows.  Alias pattern vary across columns.
   alias_match = regexp(repmat(lst(:,1), [1, numel(ia)]), ...
                        repmat(reshape(aliases(ia, 2), 1, []), ...
                               [numel(lst), 1]));

   [imatch, jmatch] = find(~ cellfun('isempty', alias_match));

   pos = accumarray(reshape(ja(jmatch) + 1, [], 1), 1, ...
                    [numel(names) + 1, 1], [], 1);
   pos = cumsum(pos);

   % Replace alias name with its expanded list of canonical names.
   names = rldecode(names, diff(pos), 2);
   names(mcolon(pos(ja), pos(ja + 1) - 1)) = reshape(lst(imatch), 1, []);
end
