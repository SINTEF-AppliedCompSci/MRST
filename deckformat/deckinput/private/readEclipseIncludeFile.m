function varargout = readEclipseIncludeFile(fun, fid, dirname, rspec, varargin)
%Read an ECLIPSE INCLUDE file.
%
% SYNOPSIS:
%   [ret{1:nret}] = readEclipseIncludeFile(fun, fid, dirname, rspec, ...)
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
%   rspec   - RUNSPEC section of current simulation case.  Needed to handle
%             pathname aliases entered in the PATHS keyword.
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

   inc_rec = read_include_record(fid);

   if ischar(inc_rec)
      inc_fn = extract_filename(inc_rec, fid);
   end

   inc_fid = open_include_file(inc_fn, dirname, rspec);

   try
      % Call back to our (likely) caller with the new file.
      [varargout{1:nargout}] = fun(inc_fid, dirname, varargin{:});

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

%--------------------------------------------------------------------------

function rec = read_include_record(fid)
   tmpl = {'|+|*File+++Name*|+|'};

   % Function readDefaultedRecord knows about quoted substrings that may
   % contain slash ('/') charaters and informational comments following the
   % record termination character.  Leverage that support here.
   rec = readDefaultedRecord(fid, tmpl);

   % Extract actual INCLUDE keyword data.  Trailing termination character
   % already discarded by readDefaultedRecord.
   rec = rec{1};
end

%--------------------------------------------------------------------------

function inc_fn = extract_filename(lin, fid)
   lin = strtrim(lin);

   quotes = strfind(lin, '''');

   if ~isempty(quotes)
      % There is a quoted substring somewhere in 'lin'.  Extract that
      % substring if possible (i.e., if correctly delimited.)
      if numel(quotes) ~= 2
         fclose(fid);

         error(msgid('Include:MissingDelimiter'), ...
              ['INCLUDE argument (%s) not correctly delimited at ', ...
               'position %lu in file ''%s''.'], ...
               lin, ftell(fid), fopen(fid));
      end

      % Extract pathname portion of INCLUDE argument.
      inc_fn = lin((quotes(1) + 1) : (quotes(2) - 1));
   else
      % Extract first (hopefully only) non-blank portion of 'lin'.
      inc_fn = sscanf(lin, '%s');
   end
end

%--------------------------------------------------------------------------

function inc_fid = open_include_file(inc_fn, dirname, rspec)
   if strcmp(inc_fn(end), '/')
      inc_fn = inc_fn(1 : end - 1);
   end

   inc_fn = normalise_filename(inc_fn, dirname, rspec);

   [inc_fid, msg] = fopen(inc_fn, 'rt');
   if inc_fid < 0
      if isunix() && ~ismac()
         % Case-sensitive system (i.e., Linux/Unix).  Simulation model MAY
         % have been created on a case insensitive system (Windows) or case
         % preserving (macOS) system.  TRY to open 'inc_fn' while ignoring
         % filename casing.
         inc_fid = open_include_file_case_insensitively(inc_fn, dirname);
      else
         error('Open:Failed', ...
               'Failed to Open INCLUDE file ''%s'': %s', inc_fn, msg);
      end
   end
end

%--------------------------------------------------------------------------

function inc_fn = normalise_filename(inc_fn, dirname, rspec)
   if ~isempty(regexp(inc_fn, regexptranslate('escape', '\'), 'once'))
      % The filename has Windows-style directory name separators.
      % Guarantee forward slashes only.
      inc_fn = regexprep(inc_fn, regexptranslate('escape', '\'), '/');
   end

   if ~isempty(regexp(inc_fn, dollar(), 'once'))
      if isfield(rspec, 'PATHS')
         inc_fn = substitute_path_aliases(inc_fn, rspec.PATHS);
      else
         error('PathAlias:Missing', ...
              ['Include File Name ''%s'' References a Path Alias ', ...
               'But Simulation Case Does Not Define PATHS'], inc_fn);
      end
   end

   % Replace forward slashes with native directory name separators (no
   % change on Unix/Linux/MacOS X).
   inc_fn = regexprep(inc_fn, '/', filesep);

   if ~strcmp(inc_fn(1), filesep)
      % Translate relative pathname to absolute pathname.
      inc_fn = fullfile(dirname, inc_fn);
   end
end

%--------------------------------------------------------------------------

function inc_fid = open_include_file_case_insensitively(inc_fn, dirname)
   assert (isunix() && ~ismac(), 'Internal Logic Error');

   [root, pth, search_fn] = split_filename_path(inc_fn, dirname);

   fname = search_filename_case_insensitively(root, pth, search_fn);

   inc_fid = open_case_insensitive_filename(fname, inc_fn);
end

%--------------------------------------------------------------------------

function [root, pth, search_fn] = split_filename_path(inc_fn, dirname)
   if ~strcmp(dirname(end), '/'), dirname = [ dirname, '/' ]; end

   [pth, search_fn, ext] = fileparts(regexprep(inc_fn, dirname, ''));

   if isempty(pth) || ~strcmp(pth(1), '/')
      root = dirname;
   else
      root = '/';  pth = pth(2:end);
   end

   search_fn = [search_fn, ext];
end

%--------------------------------------------------------------------------

function fname = ...
      search_filename_case_insensitively(fname, pth, search_fn)

   for e = [ regexp(pth, '/', 'split'), { search_fn } ]
      if isempty(e) || isempty(e{1}) || all(isspace(e{1})), continue; end

      d     = dir(fname);
      elems = { d.name };
      ix    = find(strcmpi(elems, e{1}));

      if numel(ix) == 1
         fname = fullfile(fname, elems{ix});
      else
         error('CaseMatch:Fail', ...
              ['Case Insensitive Filename Mathcing Failed Trying ', ...
               'to Match Filename Component ''', e{1}, '''']);
      end
   end
end

%--------------------------------------------------------------------------

function inc_fid = open_case_insensitive_filename(fname, inc_fn)
   if exist(fname, 'file')
      [inc_fid, msg] = fopen(fname, 'rt');

      if inc_fid < 0
         error('CaseMatch:OpenFailure', ...
               'Failed to Open Case Insensitive Filename ''%s'': %s', ...
               fname, msg);

      else
         dispif(mrstVerbose(), ...
               ['Case sensistive INCLUDE filename ''%s'' did not ', ...
                'name existing file on disk.\n  -> Using case ', ...
                'insensitive filename ''%s'' instead.\n\n'], ...
                inc_fn, fname);
      end
   else
      error('CaseMatch:None', ...
           ['Case Insensitive Filename ''%s'' Does Not Name an ', ...
            'Existing Filename on Disk'], fname);
   end
end

%--------------------------------------------------------------------------

function inc_fn = substitute_path_aliases(inc_fn, paths)
   for alias = regexp(inc_fn, [dollar(), '(\w{1,8})'], 'tokens')
      inc_fn = substitute_single_alias(inc_fn, alias{1}{1}, paths);
   end
end

%--------------------------------------------------------------------------

function patt = dollar()
   patt = regexptranslate('escape', '$');
end

%--------------------------------------------------------------------------

function inc_fn = substitute_single_alias(inc_fn, alias, paths)
   i = strcmp(alias, paths(:,1));

   nmatch = sum(i);

   if nmatch == 1
      inc_fn = strrep(inc_fn, ['$', alias], paths{i, 2});
   else
      substitution_failure(nmatch, alias);
   end
end

%--------------------------------------------------------------------------

function substitution_failure(nmatch, alias)
   msg = ['Unable to Substitute Path Alias ''$', alias, ''': '];

   if nmatch < 1
      msg = [ msg, 'No Matching Alias in PATHS Keyword' ];
   else
      msg = [ msg, 'More Than One Matching Alias in PATHS Keyword' ];
   end

   error('PathAlias:SubstituteFailure', msg);
end
