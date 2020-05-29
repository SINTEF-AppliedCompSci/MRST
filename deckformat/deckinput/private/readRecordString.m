function data = readRecordString(fid)
%Read single record string data (unconverted).
%
% SYNOPSIS:
%   rstring = readRecordString(fid)
%
% PARAMETERS:
%   fid     - Valid file identifier as obtained from FOPEN.
%
% RETURNS:
%   rstring - Record string.  Unconverted, raw character string.
%             Corresponds to single record data, not including record
%             terminator (i.e., '/' character).
%
% SEE ALSO:
%   `fopen`, `fgetl`.

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

   done = false;
   data = [];

   while ~done
      [lin, done] = get_line(fid, done);

      if ~done
         [data, done] = append_line(data, lin);
      end
   end

   p = strfind(data, '/') - 1;
   if ~isempty(p)
      data = data(1 : p(end));              % Exclude terminator character.
   end
end

%--------------------------------------------------------------------------

function [lin, done] = get_line(fid, done)
   lin   = fgetl(fid);
   valid = ischar(lin);

   if ~valid
      % Unable to input character data.  Try to determine cause of failure.

      if feof(fid)
         % Record terminated by end of file.  Mark as "done".
         done = true;

      else
         % Some other input error--disk on fire?
         [msg, errnum] = ferror(fid);

         assert (errnum ~= 0, ...
                 'Internal error in ''%s'' while reading ''%s''.\n', ...
                 mfilename, fopen(fid));

         error(['Input failed in ''%s'' while reading ''%s''.\n', ...
                'System reports: %s.'], ...
               mfilename, fopen(fid), msg);
      end
   end
end

%--------------------------------------------------------------------------

function [data, done] = append_line(data, lin)
   comment = '(--|#).*$';

   if regexp(lin, ['^\s*', comment])
      % Line is a comment (first non-blank is a comment marker).  Skip.
      done = false;
      return
   end

   S = splitQuotedString(lin);
   i = find(~ cellfun(@isempty, regexp(S.unquoted, comment)), 1, 'first');

   if ~ isempty(i)
      S.quoted   = S.quoted  (1 : (i - 1));
      S.unquoted = S.unquoted(1 : i);

      % Remove comments.
      S.unquoted(end) = regexprep(S.unquoted(end), comment, '');
   end

   done = is_complete(S.unquoted{end});                  % Record complete?

   if done
      % Discard everything after FIRST '/' character in LAST UNquoted
      % component of 'lin'.  This is to handle cases like
      %
      %   EQUALS
      %      'MULTZ'   0.05  14  25  26  30  10  10  /  C-segm mid/B-2H
      %
      % which features a '/' character in otherwise informational text.
      p = strfind(S.unquoted{end}, '/');

      assert (~isempty(p), 'Internal Logic Error');

      S.unquoted{end} = S.unquoted{end}(1 : p(1));
   end

   lin  = assembleString(S);
   data = [data, ' ', lin];
end

%--------------------------------------------------------------------------

function tf = is_complete(lin)
   persistent HAVE_CONTAINS

   if isempty(HAVE_CONTAINS)
      HAVE_CONTAINS = exist('contains', 'builtin');
   end

   if HAVE_CONTAINS
      % Contains was introduced in MATLAB 9.1.0 (R2016b)
      tf = contains(lin, '/');
   else
      tf = ~isempty(strfind(lin, '/'));
   end
end
