function [input, pos, rec] = matchWells(lookup, source)
%Match wells to keyword records
%
% SYNOPSIS:
%   [input, pos, rec] = matchWells(lookup, source)
%
% PARAMETERS:
%   lookup - Cell array of strings defining what wells to look for.  The
%            strings must be either well names or well templates.  Well
%            lists are currently not supported.
%
%   source - List of wells to match by name.
%
% RETURNS:
%   input - Array of numeric indicies into 'lookup' that specify which of
%           the 'lookup' records where successfully matched.  The unmatched
%           lookups can be determined by code like
%
%               matched        = false([numel(lookup), 1]);
%               matched(input) = true;
%               unmatched      = find(~ matched);
%
%   pos   - Indirection array into source record list.  Size equal to
%           (NUMEL(input) + 1)-by-1.  Those 'source' records that are
%           matched by lookup record 'input(i)' are found in positions
%
%               pos(i) : (pos(i + 1) - 1)
%
%           of the record list.
%
%   rec   - Source record list.  Array of numeric indices into 'source'
%           that are matched by records in 'lookup'.
%
% NOTE:
%   Well names and well template names are matched case insensitively.

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

   assert (iscellstr(lookup) && iscellstr(source), ...
          ['Inputs ''lookup'' and ''source'' must both be ', ...
           'cell arrays of strings']);

   [input, pos, rec] = deal([], 1, []);

   for i = 1 : numel(lookup),
      if is_well_list(lookup{i}),
         dispif(mrstVerbose, ...
                'Well List (''%s'') unsupported. Ignored.\n', lookup{i});
         continue
      end

      r = match(lookup{i}, source);

      if ~ isempty(r),
         [input, pos, rec] = append(input, pos, rec, i, r);
      end
   end
end

%--------------------------------------------------------------------------

function r = match(lookup, source)
% Identify 'source' records matching 'lookup' record

% Note: Names matched case insensitively.

   if is_template(lookup),
      search = [lookup(1:end-1), '.*'];

      i = ~ cellfun(@isempty, regexpi(source, search));
   else
      i = strcmpi(lookup, source);
   end

   r = find(i);
end

%--------------------------------------------------------------------------

function [input, pos, rec] = append(input, pos, rec, i, r)
% Work around Code Analyzer (MLINT) messages of the form
%
%   VAR changes size on every iteration.
%   Consider preallocation for speed.

   input = [ input ; i ];
   pos   = [ pos   ; pos(end) + numel(r) ];
   rec   = [ rec   ; reshape(r, [], 1) ];
end

%--------------------------------------------------------------------------

function b = is_well_list(lookup)
% Well LISTs identified by initial asterisk character

   b = strncmpi('*', lookup, 1);
end

%--------------------------------------------------------------------------

function b = is_template(lookup)
% Templates (well or well list) identified by final asterisk character.

   b = strncmpi('*', fliplr(lookup), 1);
end
