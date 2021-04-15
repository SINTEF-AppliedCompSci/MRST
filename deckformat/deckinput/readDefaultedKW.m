function data = readDefaultedKW(fid, template, varargin)
%Read data, possibly containing default designators, for a single keyword.
%
% SYNOPSIS:
%   data = readDefaultedKW(fid, template)
%   data = readDefaultedKW(fid, template, 'pn1', pv1, ...)
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
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters.
%              The supported options are:
%                 - NRec -- Maximum number of records to input.  Integer.
%                           Default value: Inf whence input will terminate
%                           upon detecting an empty input record (usually
%                           consisting only of a slash character).
%
% RETURNS:
%   data - Keyword data.  An m-by-n CELL array, the rows of which
%          corresponds to the individual records while the columns
%          represent the individual fields.  Default field values are
%          derived from the input 'template'.
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


   opt = struct('NRec', inf);
   opt = merge_options(opt, varargin{:});

   data = [];
   rec  = readDefaultedRecord(fid, template);
   nrec = 1;
   while ~isequal(rec, template) && nrec < opt.NRec
      % Accumulate data.  We don't know how many records there are, so we
      % can unfortunately not easily pre-allocate storage for this data.
      %
      data = [data; rec];  %#ok    % Willfully disregard MLINT warnings.
      rec  = readDefaultedRecord(fid, template);
      nrec = nrec + 1;
   end

   if ~isequal(rec, template)
      % Terminated due to maximum number of records read.  Put final record
      % into return data.
      %
      data = [data; rec];
   end
end
