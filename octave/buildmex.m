function buildmex(varargin)
%Wrapper around MEX which abstracts away details of pathname generation.
% 
% This version of 'buildmex' is a patch intended for OCTAVE users.
%
% SYNOPSIS:
%   buildmex [options ...] file [files ...]
%
% DESCRIPTION:
%   `buildmex` accepts the same file and option parameters as `mex`, with
%   the additional provision that filenames which do not start with
%   `filesep` are interpreted as pathnames relative to the directory
%   containing the M file which calls `buildmex`. Absolute pathnames (those
%   starting with `filesep`) are left untouched.
%
%   This function has no return values, but a compiled mex-function with
%   the same name as the caller is produced in the directory containing the
%   M file which calls `buildmex`.
%
% PARAMETERS:
%   Various - Same parameters as `mex`. See description for more details.
%
% RETURNS:
%   Nothing - No return values.
%
% NOTES:
%   - All parameters must be character strings.
%   - Function `buildmex` requests `-largeArrayDims` from `mex`.
%   - Function `buildmex` must be called from an M file only; it cannot be
%     invoked from the base workspace.
%
% SEE ALSO:
%   `filesep`, `mex`.

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

   % Check input.
   assert (iscellstr(varargin), 'All inputs must be strings.');

   % Compute caller function and containing M file.
   st = dbstack(1, '-completenames');
   try
      caller =           st(1).name ;
      pth    = fileparts(st(1).file);
   catch %#ok
      error(msgid('CalledFrom:BASE'), ...
            'Function BUILDMEX cannot be called from BASE workspace');
   end

   args = canonicalise_filenames(pth, varargin{:});
   args = canonicalise_options  (args);

   % Support MEX option files in 'pth', but only if our caller did not
   % already specify '-f'.
   optfile = select_options_file(pth, args);
   if ~isempty(optfile) && exist(optfile, 'file') == 2
      % Option file really does exist.  Use it.
      args = [{'-f', optfile}, args];
   end

   % Determine CXX flags
   
   [args, new_cxxflags] = extract_cxx_flags_from_args(args);

   cxxflags_orig = mkoctfile('-p', 'CXXFLAGS');
   %cxxflags = [cxxflags_orig(1:end-1), ' ', new_cxxflags];
   cxxflags = new_cxxflags; % ignore the original flags - we shouldn't need them
   
   setenv('CXXFLAGS', cxxflags);
   setenv('DL_LDFLAGS', "-shared"); % get rid of the -Wl,Bsymbolic linker option

   % if regexpi(mkoctfile('-p', 'CXX'), 'g\+\+')
   %    % ensure we do not use a more recent compiler, which creates linking problems
   %    setenv('CXX', 'g++-7');
   %    setenv('LD_CXX', 'g++-7');
   % end
   
   fsep = filesep;
   if mrstVerbose()
      %mkoctfile('--mex', '-v', '-o', [pth, fsep, caller], ['-I', pth], ['-L', pth], args{:});
      mkoctfile('--mex', '-v', '-o', [pth, fsep, caller], ['-I', pth], args{:});
      %mkoctfile('--mex', '-v', '-o', [pth, fsep, caller], ['-I', pth], '-c', args{:});
      %keyboard;
   else
      %mkoctfile('--mex', '-o', [pth, fsep, caller], ['-I', pth], ['-L', pth], args{:});
      mkoctfile('--mex', '-o', [pth, fsep, caller], ['-I', pth], args{:});
   end
      
   setenv('CXXFLAGS', cxxflags_orig(1:end-1));
   
   rehash
end

%--------------------------------------------------------------------------
function [args, new_cxxflags] = extract_cxx_flags_from_args(args)

   new_cxxflags = ''; % default return value is an empty string
   
   % determine arguments that gives CXXFLAGS
   cxxflag_ind = find(cellfun(@(x) numel(x) >= 8 && strcmpi('cxxflags', x(1:8)), args));
   
   if isempty(cxxflag_ind)
      return;
   end

   assert(numel(cxxflag_ind) == 1);
   
   str = args{cxxflag_ind};
   args(cxxflag_ind) = [];
   if any(str=='-')
      new_cxxflags = str(strfind(str, '-')(1):end);
   end
end

%--------------------------------------------------------------------------

function args = canonicalise_filenames(pth, varargin)
   args = varargin;
   fsep = filesep;
   if fsep == '\'
     % when used in a regular expression, '\' must be preceded by '\'
     fsep = '\\';  
   end
   patt = [ strcat('^\s*', switchChar()), ...
            { '=', ['^\s*', fsep], '^\s*\w:\\' } ];

   % Prepend absolute pathnames to parameters not matching either pattern
   % and that actually represent file-names in the caller's output
   % directory.
   i = false([numel(args), numel(patt)]);
   for p = 1 : numel(patt)
      i(:,p) = cellfun(@isempty, regexp(args, patt{p}));
   end

   for k = reshape(find(all(i, 2)), 1, [])  % all(i,2) true for "barewords"
      fname = strcat(pth, filesep, args{k});

      if exist(fname, 'file') == 2
         args{k} = fname;
      end
   end
end

%--------------------------------------------------------------------------

function args = canonicalise_options(args)
   if ~ispc
      % Nothing to do if we're not on MS Windows (MSVC).
      return
   end

   % MSVC uses leading '/' to identify command line options (e.g. '/MD' for
   % a multithreaded DLL linking to MSVCRT.DLL).  This conflicts with MEX's
   % internal filename handling logic so replace the '/' character with the
   % character '-'. This is also supported by MSVC for compatibility.
   args = regexprep(args, '^\s*/', '-');
end

%--------------------------------------------------------------------------

function optfile = select_options_file(pth, args)
   optfile = [];
   entries = dir([pth, filesep, 'mexopts*']);

   if ~any(strcmpi(args, '-f')) && ~isempty(entries)
      % At least one 'mexopts' file exists in 'pth'.  Construct canonical,
      % architecture dependent mexopts filename, then check if this file is
      % among the 'entries'.

      optfile = fullfile(pth, 'mexopts');
      if ispc
         optfile = [optfile, '.bat'];  % Windows
      else
         optfile = [optfile, '.sh' ];  % Unix or Mac.
      end
   end
end

%--------------------------------------------------------------------------

function sc = switchChar()
   sc = { '-' };  % MEX and/or Unix-like compilers
   if ispc
      sc = [ sc, { '/' } ]; % MSVC switch character
   end
end
