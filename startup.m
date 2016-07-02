function startup
%Amend MATLAB PATH to handle MRST implementation.

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

% $Date: 2014-01-10 10:30:11 +0100 (Fri, 10 Jan 2014) $
% $Revision: 12139 $

   d = fileparts(mfilename('fullpath'));
   m = fullfile(d, 'modules');

   p = split_path(genpath(d));

   i =     strncmp(m, p, length(m));
   i = i | ~cellfun(@isempty, regexp(p, '\.(git|hg|svn)'));
   i = i | ~cellfun(@isempty, regexp(p, '3rdparty'));
   i = i | cellfun(@isempty, p);

   addpath(p{~i});

   % Add modules as module root directory
   mrstPath('addroot', m);

   % Register known third-party modules
   mod_3rdparty = { 'matlab_bgl' };
   thirdparty   = @(m) fullfile(d, 'utils', '3rdparty', m);
   for mod = reshape(mod_3rdparty, 1, []),
      mrstPath('add', mod{1}, thirdparty(mod{1}));
   end

   % If there exists a startup_user.m file in the root dir of MRST, we
   % execute this file.
   local = fullfile(ROOTDIR, 'startup_user.m');
   if exist(local, 'file') == 2
       run_local(local);
   end

   % Automatically load selected modules for backwards compatibility.
   autoload = { 'incomp' };
   p = mrstPath('search', autoload{:});

   if isempty(p),
      autoload = {};
   elseif iscellstr(p),
      autoload = autoload(cellfun(@isempty, p));
   end

   if ~isempty(autoload),
      pl = 's'; if numel(autoload) == 1, pl = ''; end

      fprintf(['Note: Automatically loading selected ', ...
               'module%s for backwards compatibility:\n'], pl);

      fprintf('  * %s\n', autoload{:});

      mrstModule('add', autoload{:})
   end
   % Display welcome message
   mrstStartupMessage();
end

%--------------------------------------------------------------------------

function p = split_path(p)
   try
      p = regexp(p, pathsep, 'split');
   catch  %#ok
      % Octave compatibility.  It is an error to get here in an M run.
      p = strsplit(p, pathsep);
   end
end

%--------------------------------------------------------------------------

function run_local(local)
   run(local)
end
