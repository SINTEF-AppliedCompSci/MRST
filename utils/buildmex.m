function buildmex(varargin)
%Wrapper around MEX which abstracts away details of pathname generation.
%
% SYNOPSIS:
%   buildmex [options ...] file [files ...]
%
% PARAMETERS:
%   BUILDMEX accepts the same file and option parameters as MEX, with the
%   additional provision that filenames which do not start with FILESEP are
%   interpreted as pathnames relative to the directory containing the M
%   file which calls BUILDMEX.  Absolute pathnames (those starting with
%   FILESEP) are left untouched.
%
%   Note: All parameters must be character strings.
%   Note: Function BUILDMEX requests -largeArrayDims from MEX.
%
% RETURNS:
%   Nothing, but a compiled mex-function with the same name as the caller
%   is produced in the directory containing the M file which calls BUILDMEX.
%
% CAVEAT:
%   Function BUILDMEX must be called from an M file only; it cannot be
%   invoked from the base workspace.
%
% SEE ALSO:
%   filesep, mex.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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

   args = varargin;
   patt = {'^\s*-', '=', ['^\s*', filesep], '^\s*\w:\\'};

   % Prepend absolute pathnames to parameters not matching either pattern.
   i = false([numel(args), numel(patt)]);
   for p = 1 : numel(patt),
      i(:,p) = cellfun(@isempty, regexp(args, patt{p}));
   end
   i = all(i, 2);
   args(i) = strcat(pth, filesep, args(i));

   % Support MEX option files in 'pth', but only if our caller did not
   % already specify '-f'.
   entries = dir([pth, filesep, 'mexopts*']);
   if ~any(strcmpi(args, '-f')) && ~isempty(entries),
      % At least one 'mexopts' file exists in 'pth'.  Construct canonical,
      % architecture dependent mexopts filename, then check if this file is
      % among the 'entries'.

      optfile = fullfile(pth, 'mexopts');
      if ispc,
         optfile = [optfile, '.bat'];  % Windows
      else
         optfile = [optfile, '.sh' ];  % Unix or Mac.
      end

      if exist(optfile, 'file') == 2,
         % Option file really does exist.  Use it.
         args = [{'-f', optfile}, args];
      end
   end

   % Setup complete, now let MEX do its thing...
   mex('-outdir', pth, '-output', caller, ...
       ['-I', pth], ['-L', pth], '-largeArrayDims', args{:});
end
