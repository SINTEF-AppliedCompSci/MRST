function varargout = mrstModule(varargin)
%Query or modify list of activated add-on MRST modules
%
% SYNOPSIS:
%   Either of the modes
%      1) mrstModule <command> [module list]
%      2) modules = mrstModule
%
% PARAMETERS:
%     command     - *Mode 1 only*: One of the explicit verbs 'add', 'clear',
%                   'list' or 'reset'. The semantics of the command verbs
%                   are as follows:
%
%                       * add - Activate specified modules from the
%                         [module list].  Modules already activated are
%                         moved to the beginning of MATLAB's search path
%                         and remain active.
%
%                       * clear - Deactivate all modules.  An explicit
%                         module list, if present, is ignored.
%
%                       * list - Display list of currently active modules
%                         in command window.  An explicit module list, if
%                         present, is ignored.
%
%                       * reset - Convenience verb.  Equivalent to the
%                         verb sequence::
%                           mrstModule clear
%                           mrstModule add [module list]
%
%                       * gui - Launch user interface for loading and
%                         unloading modules. Equivialent to calling
%                         `moduleGUI()` directly.
%
%     module_list - *Mode 1 only*: A sequence of strings naming individual
%                   add-on modules for MRST.  The mapping of module names
%                   to system paths is performed by function `mrstPath`.
%
% RETURNS:
%     modules - *Mode 2 only*: List, represented as a cell array of
%               strings, of the currently active add-on modules.
%
% EXAMPLES:
%    mrstModule add deckformat ad-core spe10
%
%    mrstModule list
%    mrstModule clear
%
% SEE ALSO:
%   `mrstPath`.

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

   persistent MODLIST
   persistent BEHAVIOUR

   % Enforce cell array type.
   if isempty(MODLIST), MODLIST = {}; end

   if isempty(BEHAVIOUR), BEHAVIOUR = 'experimental'; end

   if mrstPlatform('matlab')
      purge_octave_only();
   end

   if nargin > 0
      assert (iscellstr(varargin), 'All parameters must be strings.');

      cmd  = varargin{1};
      mods = varargin(2 : end);

      switch lower(cmd)
         case 'add'
            MODLIST = prune_modules(MODLIST);

            mods = deduplicate(cmd, mods);

            if ~ (isempty(mods) || all(cellfun(@isempty, mods)))
               MODLIST = add_modules(MODLIST, mods, BEHAVIOUR);
               mlock
            else
               munlock
            end

         case 'clear'
            munlock
            MODLIST = clear_modules(MODLIST, BEHAVIOUR);

         case 'list'
            MODLIST = prune_modules(MODLIST);

            if ~ isempty(MODLIST)
               fprintf('Currently active MRST modules\n');
               fprintf('  * %s\n', MODLIST{:});
            else
               if mislocked, munlock, end
               fprintf('No active MRST modules\n');
            end

         case 'reset'
            munlock
            MODLIST = clear_modules(MODLIST, BEHAVIOUR);

            if ~ (isempty(mods) || all(cellfun(@isempty, mods)))
               mods = deduplicate(cmd, mods);

               mlock
               MODLIST = add_modules(MODLIST, mods, BEHAVIOUR);
            end

         case 'gui'
            moduleGUI();

         case {'behavior', 'behaviour'}
            if nargin == 1
               % Print to screen the current behaviour
               fprintf('Current mrstModule behaviour is "%s".\n', BEHAVIOUR);
            else
               % We have been given a new behaviour
               nextmode = lower(mods{1});
               if ~strcmp(nextmode, BEHAVIOUR)
                  % The behaviour has changed.
                  munlock

                  assert (any(strcmp(nextmode, {'experimental', 'release'})), ...
                         ['Supported behaviours are ''experimental'' ', ...
                          'and ''release''']);
                  mods = MODLIST;

                  % Clear modules with previous behavior mode
                  MODLIST = clear_modules(MODLIST, BEHAVIOUR);

                  % Set new behaviour mode
                  BEHAVIOUR = nextmode;

                  % Add modules, this time changing the included/excluded
                  % folders according to the current directive.
                  mods = deduplicate(cmd, mods);

                  MODLIST = add_modules(MODLIST, mods, BEHAVIOUR);

                  mlock
               end
            end

         otherwise
            error(msgid('Command:Unsupported'), ...
                 ['Command word ''%s'' unsupported. Must be one of ', ...
                  '''add'', ''clear'', ''list'', ''gui'', or ''reset''.'], cmd);
      end

   elseif nargout == 0

      mrstModule list

   elseif nargout == 1

      MODLIST = prune_modules(MODLIST);
      varargout{1} = MODLIST;

   else

      error(msgid('Syntax:Error'), ...
           ['Call syntax is\n\t', ...
            mfilename, ' <command> [module list]  or\n\t', ...
            'mods = ', mfilename]);
   end
end

%--------------------------------------------------------------------------

function mods = deduplicate(cmd, mods)
   [umod, iu, iu] = unique(mods);                               %#ok<ASGLU>

   if numel(umod) == numel(mods)
      % No duplicates in module list (common case).  Return input list
      % unaltered.
      return
   end

   % If we get here, the list ist of input modules contains duplicates.
   % Prune those (moderately expensive) and warn that we've altered the
   % list (to preserve the *last* of each duplicated entry).

   dup = umod(accumarray(iu, 1) > 1);

   assert (~isempty(dup), 'Internal Logic Error');
   pl = ''; if numel(dup) ~= 1, pl = 's'; end

   warnmsg = sprintf(['Encountered Duplicated Module%s in ''%s'' ', ...
                      'Command Verb ''%s'' (Using Last Duplicated ', ...
                      'Entry):\n'], pl, mfilename, cmd);
   warnmsg = [ warnmsg, sprintf('  * %s\n', dup{:}) ];

   warning('ModuleNames:Duplicated', warnmsg);

   [i, i] = sort(accumarray(reshape(iu, [], 1), (1 : numel(iu)).', ...
                            [numel(umod), 1], @max));           %#ok<ASGLU>
   mods   = umod(i);
end

%--------------------------------------------------------------------------

function lst2 = prune_modules(lst)
   lst2 = {};

   if ~isempty(lst)
      pth  = split_path(path);
      mdir = path_search(lst);

      [i, j] = blockDiagIndex(numel(pth), numel(mdir));
      match  = strcmp(reshape(pth (i), [], 1), ...
                      reshape(mdir(j), [], 1));

      lst2 = lst(unique(j(match)));
   end
end

%--------------------------------------------------------------------------

function lst = add_modules(lst, mods, behaviour)
   pth = path_search(mods);
   fnd = ~ cellfun(@isempty, pth);

   if any(~ fnd)
      fprintf('No module mapping found for\n');
      fprintf('  * %s\n', mods{~fnd});
   end

   for k = reshape(find(fnd), 1, [])
      mroot = pth{k};

      try
         mload   = fullfile(mroot, 'private', 'modload.m');
         mloadfb = fullfile(mroot, 'private', 'modload_fallback.m');

         if exist(mload, 'file') == 2
            % If module provides a 'modload', then run it
            run(mload);

         elseif exist(mloadfb, 'file') == 2
            % Otherwise, if module provides a 'modload_fallback', then run
            % that.  This is an escape clause that supports code generation
            % approaches to constructing 'modload' without putting the
            % generated code into a VCS (e.g., Git).
            %
            run(mloadfb);
         end

         dirs = filter_module_dirs(mroot, behaviour);
         addpath(dirs{:});

      catch ME
         % 'modload' failed, proceed to next module.
         [m, m] = fileparts(mroot);                             %#ok<ASGLU>

         fprintf('Failed to load ''%s'':\n%s\n', m, ME.message);
         fprintf('Module ''%s'' IGNORED\n', m);

         fnd(k) = false;

         continue;
      end
   end

   if ~ any(fnd)
      % No pathname known for any of the requested modules (or they all
      % failed to load).  Return without attempting to modify the list of
      % active modules.
      return
   end

   mods   = mods(fnd);
   [i, j] = blockDiagIndex(numel(lst), numel(mods));

   match = strcmp(reshape(lst (i), [], 1), ...
                  reshape(mods(j), [], 1));

   lst(unique(i(match))) = [];
   lst = [ fliplr(reshape(mods, 1, [])), lst ];
end

%--------------------------------------------------------------------------

function lst = clear_modules(lst, behaviour)
   if ~isempty(lst)
      p = path_search(lst);
      p = p(~ cellfun(@isempty, p));

      for r = reshape(p, 1, [])
         dirs = filter_module_dirs(r{1}, behaviour);
         rmpath(dirs{:});
      end

      lst = {};
   end
end

%--------------------------------------------------------------------------

function pth = path_search(mods)
   assert (iscellstr(mods), 'Internal error in %s', mfilename);

   pth = mrstPath('search', mods{:});

   if ischar(pth)
      assert (numel(mods) == 1, 'Internal error in %s', mfilename);
      pth = { pth };
   end

   if isunix
      % Implement (partial) tilde expansion akin to FOPEN to avoid false
      % negatives in 'prune_modules'.  Paths of the form '~/<something>' is
      % by far the most common, so that's what we support.  We assume that
      % every user has a valid 'HOME' value.

      pth = regexprep(pth, '^~/', [getenv('HOME'), filesep]);
   end
end

%--------------------------------------------------------------------------

function dirs = filter_module_dirs(root, behaviour)
   if ~isempty(root)
      dirs  = split_path(genpath(root));
      vcdir = [regexptranslate('escape', filesep), '\.(git|hg|svn)'];

      is_vcdir = ~cellfun(@isempty, regexp(dirs, vcdir));
      exclude  = is_vcdir | cellfun(@isempty, dirs);

      if strcmpi(behaviour, 'release')
          exdir = [regexptranslate('escape', filesep), 'experimental'];
          is_exdir = ~cellfun(@isempty, regexp(dirs, exdir));
          exclude = exclude | is_exdir;
      end

      dirs = dirs(~exclude);
   else
      dirs = {};
   end
end

%--------------------------------------------------------------------------

function pth = split_path(pth)
   try
      pth = regexp(pth, pathsep, 'split');
   catch  %#ok
      % Octave compatiblity.  It is an error to get here in an M run.
      pth = strsplit(pth, pathsep);
   end
end

%--------------------------------------------------------------------------

function purge_octave_only()
   octave = fullfile('utils', 'octave_only');
   pth    = regexp(path, regexptranslate('escape', pathsep), 'split');
   match  = regexp(pth , regexptranslate('escape', octave));
   purge  = pth(~cellfun('isempty', match));

   if ~isempty(purge)
      rmpath(purge{:});
   end

   if mrstVerbose && ~isempty(purge)
      dstring = sprintf(' * %s\n', purge{:});

      pl = 'ies were';
      if sum(purge) == 1, pl = 'y was'; end

      warning('ExlcudeOctave:InMATLAB', ...
             ['The following Octave compatibility director%s ', ...
              'removed from the search path due to running ', ...
              'in MATLAB:\n%s'], pl, dstring);
   end
end
