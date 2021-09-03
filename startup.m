function startup
%Amend MATLAB PATH to handle MRST implementation.

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

   build_mrst_path_tree();

   % Register known third-party modules
   mod_3rdparty = { 'matlab_bgl' };
   activate_3rdparty_modules(mod_3rdparty);

   % If there exists a startup_user.m file in the root dir of MRST, we
   % execute this file.
   run_local();

   % Run platform (Octave/Matlab) specific startup routines
   run_platform_specific();
   
   % Automatically load selected modules for backwards compatibility.
   autoload = {};
   load_compat_modules(autoload);

   % Display welcome message
   mrstStartupMessage();

   % Attempt to load settings - and perform first-time setup if required.
   settings = mrstSettings(); %#ok
end

%--------------------------------------------------------------------------

function d = rootdir
   d = fileparts(mfilename('fullpath'));
end

%--------------------------------------------------------------------------

function build_mrst_path_tree
   d = rootdir();

   m = fullfile(d, 'modules');
   p = split_path(genpath(d));

   % Remove octave_only folders from path
   oct_only = ~cellfun(@isempty, regexp(p, 'octave_only'));

   i =     strncmp(m, p, length(m));
   i = i | ~cellfun(@isempty, regexp(p, '\.(git|hg|svn)'));
   i = i | ~cellfun(@isempty, regexp(p, '3rdparty'));
   i = i | ~cellfun(@isempty, regexp(p, 'output'));
   i = i | cellfun(@isempty, p);
   i = i | oct_only;

   addpath(p{~i});

   if mrstPlatform('octave')
       % If we are on Octave, we add in the octave_only paths at the end to
       % ensure that they can overwrite any of the core files where needed.
       % We first turn off the shadowing warning since we may have
       % compatibility fixes included in this folder.
       os = 'Octave:shadowed-function';
       warnstate = warning('query', os);
       warning('off', os);
       addpath(p{oct_only});
       warning(warnstate);
   end

   % Add modules as module root directory
   mrstPath('addroot', m);
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

function activate_3rdparty_modules(mod_3rdparty)
   d          = rootdir();
   thirdparty = @(m) fullfile(d, 'utils', '3rdparty', m);

   for mod = reshape(mod_3rdparty, 1, [])
      mrstPath('add', mod{1}, thirdparty(mod{1}));
   end
end

%--------------------------------------------------------------------------

function load_compat_modules(mlist)
   if isempty(mlist), return, end

   p = mrstPath('search', mlist{:});

   if isempty(p)
      mlist = {};
   elseif iscellstr(p)
      mlist = mlist(~ cellfun(@isempty, p));
   end

   if ~isempty(mlist)
      pl = 's'; if numel(mlist) == 1, pl = ''; end

      fprintf(['Note: Automatically loading selected ', ...
               'module%s for backwards compatibility:\n'], pl);

      fprintf('  * %s\n', mlist{:});

      mrstModule('add', mlist{:})
   end
end

%--------------------------------------------------------------------------

function run_local()
   do_run_local(fullfile(rootdir(), 'startup_user.m'));
end

%--------------------------------------------------------------------------

function run_platform_specific()
   if mrstPlatform('octave')
      do_run_local(fullfile(rootdir(), 'utils', ...
                            'octave_only', 'startup_octave.m'));
   end
end

%--------------------------------------------------------------------------

function do_run_local(local)
   if exist(local, 'file') == 2
      run(local);
   end
end
