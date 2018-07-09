function varargout = readEclipseIncludeFile(fun, fid, dirname, varargin)
%Read an ECLIPSE INCLUDE file.
%
% SYNOPSIS:
%   [ret{1:nret}] = readEclipseIncludeFile(fun, fid, dirname, ...)
%
% PARAMETERS:
%   fun     - Callback function handle.  Assumed to support the syntax
%
%                [ret{1:nret}] = fun(fid, dirname, ...)
%
%             Function 'fun' will be called with 'fid' set to the FOPEN
%             return value of the INCLUDE file name, while 'dirname' will
%             be the complete directory name of the INCLUDE file name.
%
%   fid     - Valid file identifier as obtained by FOPEN.  The file pointer
%             FTELL(fid) is assumed to be placed directly after the
%             'INCLUDE' keyword and strictly before the name of the file
%             which will be INCLUDEd.
%
%             Function 'readEclipseIncludeFile' will, upon successful
%             return, place FTELL(fid) directly after the keyword-closing
%             slash character.
%
%   dirname - Complete directory name of file from which the input file
%             identifier 'fid' was derived through FOPEN.
%
%   ...     - Additional function parameters.  These parameters will be
%             passed unchanged on to function 'fun'.
%
% RETURNS:
%   ret     - Any and all return values from function 'fun'.  It is the
%             responsibility of the caller of 'readEclipseIncludeFile' to
%             supply sufficient number of output arrays (i.e., 'nret') to
%             store these return values.

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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


   lin = '';
   while ischar(lin) && isempty(strtrim(lin)),
      lin = fgetl(fid);
      lin = regexprep(lin, '--.*$', '');
   end

   if ischar(lin),
      lin = strtrim(lin);

      quotes = strfind(lin, '''');

      if ~isempty(quotes)
         % There is a quoted substring somewhere in 'lin'.  Extract that
         % substring if possible (i.e., if correctly delimited.)
         if numel(quotes) ~= 2
            fclose(fid);
            error('INCLUDE argument (%s) not correctly delimited.', lin);
         end

         % Extract pathname portion of INCLUDE argument.
         inc_fn = lin((quotes(1) + 1) : (quotes(2) - 1));
      else
         % Extract first (hopefully only) non-blank portion of 'lin'.
         inc_fn = sscanf(lin, '%s');
      end

      p = strfind(lin, '/');
      if isempty(p)
         terminated = false;
      else
         p = p(end);
         terminated = (p == numel(lin)) || isspace(lin(p + 1));
      end
   end

   if strcmp(inc_fn(end), '/'), inc_fn = inc_fn(1 : end - 1); end

   % Gobble up keyword-closing '/' character if not already read.
   if ~terminated,
      p     = ftell(fid);
      fn    = fopen(fid);

      slash = fscanf(fid, '%s', 1);  % Possibly too weak.

      if ~strcmp(slash, '/')
         fclose(fid);
         error(msgid('Include:WrongfulTermination'), ...
              ['INCLUDE keyword not correctly terminated at ', ...
               'position %lu in file ''%s'''], p, fn);
      end
   end

   if ~isempty(regexp(inc_fn, regexptranslate('escape', '\'), 'once'))
      % The filename has Windows-style directory name separators.
      % Guarantee forward slashes only.
      inc_fn = regexprep(inc_fn, regexptranslate('escape', '\'), '/');
   end

   % Replace forward slashes with native directory name separators (no
   % change on Unix/Linux/MacOS X).
   inc_fn = regexprep(inc_fn, '/', filesep);

   if ~strcmp(inc_fn(1), filesep)
      % Translate relative pathname to absolute pathname.
      inc_fn = fullfile(dirname, inc_fn);
   end

   [inc_fid, msg] = fopen(inc_fn, 'rt');
   if inc_fid < 0, error([inc_fn, ': ', msg]); end

   try
      % Call back to our (likely) caller with the new file.
      [varargout{1:nargout}] = fun(inc_fid, ...
                                   dirname, ...  % or fileparts(inc_fn)
                                   varargin{:});
   catch %#ok
      err = lasterror;  %#ok
      try  %#ok
         % Don't leak fids, but don't gripe about a child already closing
         % the fid.
         fclose(inc_fid);
      end
      rethrow(err);
   end
   fclose(inc_fid);
end
