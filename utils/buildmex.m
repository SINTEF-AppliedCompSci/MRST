function buildmex(varargin)
%Wrapper around MEX which abstracts away details of pathname generation.
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

   % Check input.
   assert (iscellstr(varargin), 'All inputs must be strings.');

   isOctave = mrstPlatform('octave');
   [pth, caller] = identify_caller(dbstack(1, '-completenames'), isOctave);

   args = canonicalise_filenames(isOctave, pth, varargin{:});
   args = canonicalise_options  (args);
   if isOctave
      buildmex_octave(pth, args, caller);
   else
      buildmex_matlab(pth, args, caller);
   end

   if ispc || isOctave
      % Unconditionally expose output (= caller) to MATLAB environment.
      % This is needed on Windows to avoid infinite loops when BUILDMEX is
      % used to affect automatic compilation of MEX-based accelerators
      % through a gateway routine of the form
      %
      %     function varargout = accel_mex(varargin)
      %        buildmex(...)
      %
      %        [varargout{1:nargout}] = accel_mex(varargin{:})
      %     end
      %
      % Without this sledgehammer (i.e., unconditionally invoking REHASH)
      % MATLAB does not detect that there is a new (MEXed) version of a
      % particular M-language gateway function and the recursive call ends
      % up in the gateway routine rather than the accelerator.
      %
      rehash();
   end
end

%--------------------------------------------------------------------------

function [pth, caller] = identify_caller(st, isOctave)
   if isOctave
      % some versions of octave include 'buildmex' in the stack trace
      % returned by 'dbstack'
      if strcmpi(st(1).name, 'buildmex')
         st = st(2:end);
      end
   end

   try
      caller =           st(1).name ;
      pth    = fileparts(st(1).file);
   catch
      error(msgid('CalledFrom:BASE'), ...
            'Function BUILDMEX cannot be called from BASE workspace');
   end
end

%--------------------------------------------------------------------------

function buildmex_matlab(pth, args, caller)
% Support MEX option files in 'pth', but only if our caller did not
% already specify '-f'.

   optfile = select_options_file(pth, args);
   if ~isempty(optfile) && (exist(optfile, 'file') == 2)
      % Option file really does exist.  Use it.
      args = [{'-f', optfile}, args];
   end

   verbose = '';
   if mrstVerbose()
      verbose = '-v';
   end

   % Setup complete, now let MEX do its thing...
   mex(verbose, '-outdir', pth, '-output', caller, ...
       ['-I', pth], ['-L', pth], '-largeArrayDims', args{:});
end

%--------------------------------------------------------------------------

function buildmex_octave(pth, args, caller)
% Support MEX option files in 'pth', but only if our caller did not
% already specify '-f'.

   optfile = select_options_file(pth, args);
   if ~isempty(optfile) && exist(optfile, 'file') == 2
      % Option file really does exist.  Use it.
      args = [{'-f', optfile}, args];
   end

   % Determine CXX flags
   [args, new_cxxflags] = extract_cflags_from_args('cxxflags', args);
   [args, new_cflags] = extract_cflags_from_args('cflags', args);

   cxxflags_orig = mkoctfile('-p', 'CXXFLAGS');
   cflags_orig = mkoctfile('-p', 'CFLAGS');

   %cxxflags = [cxxflags_orig(1:end-1), ' ', new_cxxflags];
   cxxflags = new_cxxflags; % ignore the original flags - we shouldn't need them
   cflags = [new_cflags, ' -std=c99']; % add c99 flag to avoid problems with mex.h

   % containing c++-style comments
   setenv('CXXFLAGS', cxxflags);
   setenv('CFLAGS', cflags);
   % setenv('DL_LDFLAGS', "-shared"); % get rid of the -Wl,Bsymbolic linker option

   % if regexpi(mkoctfile('-p', 'CXX'), 'g\+\+')
   %    % ensure we do not use a more recent compiler, which creates linking problems
   %    setenv('CXX', 'g++-7');
   %    setenv('LD_CXX', 'g++-7');
   % end

   fsep = filesep;
   fp = [pth, fsep, caller];
   if exist([fp, '.cpp'], 'file')
      cpath = [fp, '.cpp'];
   else
      cpath = [fp, '.c'];
   end

   if is_octfile(cpath)
      fprintf('Building OCT-file\n');
      arg = {'-DMRST_OCTEXT'};
   else
      fprintf('Building MEX-file\n');
      arg = {'--mex'};
   end

   if mrstVerbose()
      arg = [arg, {'-v'}];
   end

   [out, status] = mkoctfile(arg{:}, '-o', fp, ['-I', pth], args{:});
   if status ~= 0
      error('Unable to build: %s: mkoctfile returned message: "%s"', fp, out);
   end

   setenv('CXXFLAGS', cxxflags_orig(1:end-1));
   setenv('CFLAGS', cflags_orig(1:end-1));
end

%--------------------------------------------------------------------------

function args = canonicalise_filenames(isOctave, pth, varargin)
   args = varargin;

   fsep = filesep();
   if isOctave
      fsep = regexptranslate('escape', fsep);
   end

   non_bareword_patt = [ strcat('^\s*', switchChar()), ...
                         { '=', ['^\s*', fsep], '^\s*\w:\\' } ];

   % Prepend absolute pathnames to parameters not matching either pattern
   % and that actually represent file-names in the caller's output
   % directory.
   i = false([numel(args), numel(non_bareword_patt)]);
   for p = 1 : numel(non_bareword_patt)
      i(:,p) = cellfun(@isempty, regexp(args, non_bareword_patt{p}));
   end

   for k = reshape(find(all(i, 2)), 1, [])  % all(i,2) true for "barewords"
      fname = fullfile(pth, args{k});

      if exist(fname, 'file') == 2
         args{k} = fname;
      end
   end
end

%--------------------------------------------------------------------------

function ok = is_octfile(filename)
   if exist(filename, 'file')
      % If the file exists, we read the entire thing and look for the
      % magic signature (and hope noone got too create with formatting or
      % preprocessor)
      tmp = fileread(filename);
      ok = ~isempty(strfind(tmp, 'DEFUN_DLD')); %#ok
   else
      ok = false;
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
   entries = dir(fullfile(pth, 'mexopts*'));

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

%--------------------------------------------------------------------------

function [args, new_flags] = extract_cflags_from_args(flag, args)
   new_flags = ''; % default return value is an empty string

   % determine arguments that gives 'flag' (CXXFLAGS or CFLAGS)
   flag_ind = find(cellfun(@(x) numel(x) >= 8 && strcmpi(flag, x(1:numel(flag))), args));

   if isempty(flag_ind)
      return;
   end

   assert(numel(flag_ind) == 1);

   str = args{flag_ind};

   % check if multiple entries need to be concatenated
   if any(str == '"')
      % this entry is just the beginning of a longer set of options that
      % should be concatenated

      % search for end of flag list
      end_ix = find(arrayfun(@(n) any(n{:} == '"'), args(flag_ind:end)));
      if numel(end_ix) > 1
         end_ix = end_ix(2); % otherwise, remain at end_ix(1), which equals flag_ind
      end

      flagargs = args(flag_ind:end_ix);
      args(flag_ind:end_ix) = [];

      tmp = arrayfun(@(n) [n{:}, ' '], flagargs, 'UniformOutput', false);
      new_flags = horzcat(tmp{:});

      % remove quotes and what becomes before
      tmp = find(new_flags == '"');
      new_flags = new_flags(tmp(1)+1:tmp(2)-1);

      % remove possible reference to $CFLAGS or $CXXFLAGS
      new_flags = erase(erase(new_flags, '$CFLAGS'), '$CXXFLAGS');

   else
      % the specification of flags is limited to this single entry
      args(flag_ind) = [];
      if any(str == '-')
         tmp = strfind(str, '-');
         new_flags = str(tmp(1) : end);
      end
   end
end
