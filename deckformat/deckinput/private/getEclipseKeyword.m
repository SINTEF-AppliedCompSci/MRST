function kw = getEclipseKeyword(fid)
%Extract the next keyword from an ECLIPSE input deck.
%
% SYNOPSIS:
%   kw = getEclipseKeyword(fid)
%
% PARAMETERS:
%   fid - Valid file identifier as obtained from FOPEN.
%         Assumed to point to a seekable input stream (i.e., not a POSIX
%         pipe(2)).
%
% RETURNS:
%   kw  - Valid ECLIPSE keyword string (stripped of any terminating
%         character), or -1 if file input failed.  When kw==-1, the
%         functions FEOF and FERROR may be used to determine the cause of
%         the input failure.
%
% NOTE:
%   Function 'getEclipseKeyword' reads the input file one line at a time.
%   This is detrimental to performance when reading large amounts of data.
%   Consequently, function 'getEclipseKeyword' should only be used to input
%   the small portion of file data following any data from a previous
%   keyword and up to and including the next keyword string.
%
%   In particular, function 'getEclipseKeyword' should generally not be
%   used as a replacement for the POSIX 'grep' utility when determining
%   which keywords are present in any given input deck.
%
% SEE ALSO:
%   `fopen`, `ftell`, `fseek`, `feof`, `ferror`.

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

   % An ECLIPSE keyword consists of at least one, and at most eight,
   % upper-case English (ASCII) characters at the beginning of a line.  Any
   % file data from column 9 up to (and including) the newline are
   % considered comments and ignored.
   %
   % The keyword must be followed by at least one white space character
   % (possibly '\n') or a keyword-terminating slash character.  Moreover,
   % any subsequent keyword data (e.g., porosity values) must start on a
   % new line.  IOW, the data cannot follow the keyword on the same line.
   %
   % Consequently, reading input using the newline-stripping FGETL function
   % will correctly position the file pointer for the purpose of any
   % further data reading once a keyword is returned to the caller.
   %
   matches_kw = @(s) ~isempty(regexp(s, '^[A-Z][A-Z0-9]{0,7}(|/)', 'once'));

   kw = fgetl(fid);
   while ischar(kw) && ~matches_kw(kw)
      kw = fgetl(fid);
   end

   if ischar(kw) && kw(end) == '/'
      % Found valid keyword, but the kw string ends in '/'.  Remove slash.
      kw = kw(1:end-1);
   end

   if ischar(kw)
      % The keyword consists solely of upper case English characters and
      % digits.  The following is usually, but not always, an expensive
      % no-op.
      kw = regexprep(kw, '([A-Z][A-Z0-9]*).*', '$1');
   end
end
