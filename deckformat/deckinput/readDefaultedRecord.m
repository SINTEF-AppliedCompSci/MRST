function rec = readDefaultedRecord(fid, rec)
%Read data, possibly containing default designators, for a single record.
%
% SYNOPSIS:
%   rec = readDefaultedRecord(fid, template)
%
% DESCRIPTIONS:
%   Reads a single record, containing a specific number of fields, whilst
%   replacing input default value designators of the form 'n*' (where 'n'
%   is an integer) with 'n' copies of a user-supplied default value.
%
%   Function 'readDefaultedRecord' supports early record termination (i.e.,
%   the terminator character '/' occurring before all fields have been
%   input) in which case any trailing fields will be assigned the default
%   value supplied by the caller.
%
% PARAMETERS:
%   fid      - Valid file identifier as obtained from FOPEN.
%
%   template - An n-element CELL array constituting a template for the next
%              record.  Assumed to contain 'n' copies of a default value.
%              The typical default is the string 'NaN', but any value may
%              be used as 'readDefaultedRecord' does not inspect this value
%              in any way.
%
%              To conveniently generate an n-element CELL array, each
%              element of which contains the string 'NaN', issue a
%              statement of the form:
%
%                  template(1 : n) = { 'NaN' }
%
% RETURNS:
%   rec - Fully, or partially, filled 'template' CELL array.  Any
%         non-defaulted input data will be entered, in *STRING* form, into
%         the corresponding 'rec' element whilst leaving defaulted data in
%         its input 'template' form.
%
% NOTE:
%   It is the caller's responsibility to supply default values in the input
%   'template' which can be reliably distinguished from all (expected) data
%   for any given keyword record.  Moreover, an all-default return value
%   (i.e., rec = template) typically constitutes the end of the keyword
%   data.
%
% SEE ALSO:
%   `readDefaultedKW`, `readWellKW`.

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

   data = readRecordString(fid);

   % An all-default (i.e., empty) record terminates the keyword.
   %
   % Only enter this part if there actually *is* any data.
   %
   if ~isempty(data) && ~all(isspace(data))
      % Tokenize data.
      data = removeQuotes(tokenizeRecord(strtrim(data)));

      % Detect "default" strings of the form n*.
      is_def = regexp(data, '^(\d+)\*$', 'tokens', 'once');
      i      = ~cellfun(@isempty, is_def);

      if ~all(i)
         % Compute index offsets for non-defaulted values.
         add    = zeros([numel(rec), 1]);
         add(i) = cellfun(@(c) sscanf(c{1}, '%f'), is_def(i));
         add    = cumsum(add);

         % Assign non-defaulted values from input data.
         ix      = (1 : sum(~i)).' + add(~i);
         rec(ix) = data(~i);
      end
   end
end
