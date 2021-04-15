function data = readFixedNumRecords(fid, template, nrec)
%Read fixed number of data records for a single keyword
%
% SYNOPSIS:
%   data = readFixedNumRecords(fid, template, nrec)
%
% DESCRIPTION:
%   Reads a fixed number of records of constant size-number of items-from
%   an input stream.  The typical use case is keywords like PVTW which have
%   an implied copy behaviour.
%
% PARAMETERS:
%   fid      - Valid file identifier as obtained from FOPEN.
%
%   template - An n-element CELL array constituting a template for the next
%              set of records.  Assumed to contain 'n' default values.  The
%              typical default is the string 'NaN', but any value may be
%              used as 'readFixedNumRecords' does not inspect this value in
%              any way.
%
%              To conveniently generate an n-element CELL array, each
%              element of which contains the string 'NaN', use a statement
%              of the form:
%
%                  template = repmat({ 'NaN' }, [1, n])
%
%  nrec      - Number of records to input.  Finite, positive integer.
%              Typically derived from a maximum or exact dimension in the
%              RUNSPEC section.
%
% RETURNS:
%   data - Keyword data.  An nrec-by-n CELL array, the rows of which
%          corresponds to the individual records while the columns
%          represent the individual fields.  Default field values are
%          derived from the input 'template' or copied from the previous
%          record if applicable.
%
% NOTE:
%   It is the caller's responsibility to supply default values in the input
%   'template' which can be reliably distinguished from all (expected) data
%   for any given keyword.
%
%   Data reading terminates when the input reader detects an empty record.
%
% SEE ALSO:
%   `readDefaultedRecord`.

%{
Copyright 2020-2021 SINTEF Digital, Mathematics & Cybernetics.

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

   check_num_records(nrec);

   data = repmat(template, [nrec, 1]);
   for rec = 1 : nrec
      data(rec, :) = readDefaultedRecord(fid, template);

      % Use previous record as template for next.  Simplifies copying.
      template  = data(rec, :);
   end
end

%--------------------------------------------------------------------------

function check_num_records(nrec)
   if ~ (isnumeric(nrec) && isscalar(nrec))
      error('NumRec:NonNumeric', ...
           ['%s: Number of records must be a scalar of numeric ', ...
            'type.  Got ''%s'' of %d elements.'], ...
            mfilename(), class(nrec), numel(nrec));
   end

   if ~isfinite(nrec)
      error('NumRec:NonFinite', ...
            '%s: Number of records must be finite.  Got %g.', nrec);
   end

   if (nrec < 1) || (fix(nrec) ~= nrec)
      error('NumRec:NonPositive', ...
           ['%s: Number of records must be a strictly ', ...
            'positive integer.  Got %g.'], nrec);
   end
end
