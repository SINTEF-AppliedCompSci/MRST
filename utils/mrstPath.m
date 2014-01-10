function varargout = mrstPath(varargin)
%Establish and maintain mapping from module names to system directory paths
%
% SYNOPSIS:
%   Either of the modes
%      1) mrstPath register list
%      2) mrstPath <command> [module list]
%      3) mrstPath search module list
%         paths = mrstPath('search', module list)
%
% PARAMETERS:
%   Mode 1)
%     Register (insert) new mappings of module names to system paths
%     (directories) into the current list.  In this case, the 'list' must
%     contain an even number of strings, the odd numbered ones being
%     interpreted as module names and the even numbered interpreted as
%     directories.  The directories must exist (i.e., function ISDIR must
%     return TRUE) to register a mapping.
%
%     The pathname resolution algorithm is as follows (note: DIRARG is one
%     of the directory arguments in the above list while MODNAME is the
%     corresponding module name):
%
%        if ISDIR(fullfile(pwd, DIRARG)),
%
%           register MODNAME fullfile(pwd, DIRARG)
%
%        elseif ISDIR(fullfile(ROOTDIR, DIRARG)),
%
%           register MODNAME fullfile(ROOTDIR, DIRARG)
%
%        elseif DIRARG is a subdirectory of a any entry in PATH
%
%           register MODNAME fullfile(PATH entry, DIRARG)
%
%        elseif isdir(DIRARG)
%           % Assume DIRARG is an absolute pathname
%           register MODNAME DIRARG
%
%        else
%           % DIRARG is not a directory, skip it
%        end
%
%     EXAMPLE:
%        mrstPath register AGMG       modules/mex/AGMG       ...
%                          matlab_bgl modules/mex/matlab_bgl ...
%                          triangle   modules/mex/triangle
%
%   Mode 2)
%     <Command> -
%        One of the explicit verbs 'addroot', 'clear', 'list', 'remove'
%        'reregister' or 'reset'.  The semantics of the command verbs are
%        as follows:
%
%           o) addroot -- Register a module root directory.  List
%                         interpreted as a list of directories, subject, the
%                         subdirectories of which will be
%
%           o) clear   -- Deactivate all modules.  An explicit module list,
%                         if present, is ignored.
%
%           o) list    -- Display list of currently registered modules in
%                         command window.  An explicit module list, if
%                         present, is ignored.
%
%           o) remove  -- Deregister selected modules.  List interpreted as
%                         module names.  No action if empty.  Unknown
%                         modules ignored.
%
%           o) reregister --
%                         Convenience verb to reestablish certain module
%                         mappings.  Equivalent to the verb sequence:
%
%                            mrstPath('remove', list{1:2:end})
%                            mrstPath register list
%
%           o) reset   -- Convenience verb.  Equivalent to the verb
%                         sequence:
%
%                            mrstPath clear
%                            mrstPath register [list]
%
%   Mode 3)
%     None.
%
% RETURNS:
%   Modes 1) and 2)
%     Nothing.
%
%   Mode 2)
%     modules - List, represented as a cell array of strings, of the
%               currently active add-on modules.
%
% SEE ALSO:
%   ROOTDIR, isdir, filesep, fullfile.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   persistent CACHE

   if isempty(CACHE), CACHE = cell([0, 2]); end

   if nargin > 0,
      assert (iscellstr(varargin), 'All parameters must be strings.');

      cmd  = varargin{1};
      mods = varargin(2 : end);

      switch lower(cmd),
         case 'addroot',
            if ~ (isempty(mods) || all(cellfun(@isempty, mods))),
               mlock
               CACHE = register_root(CACHE, mods);
            else
               munlock;
            end

         case { 'add', 'register' },
            % Register a new set of modules.  Assume that 'mods' is a list
            % of key/value pairs.

            if ~ (isempty(mods) || all(cellfun(@isempty, mods))),
               mlock
               CACHE = register_modules(CACHE, mods);
            else
               munlock
            end

         case 'clear',
            % Clear cache, thus removing all currently registered modules
            % leaving only the modules that are known by default.

            munlock
            CACHE = clear_modules;

         case 'list',

            list_modules(CACHE);

         case 'remove',
            % Remove a module or set of modules.

            CACHE = remove_modules(CACHE, mods);

            if is_empty(CACHE),
               munlock
            end

         case 'reregister',
            % Convenience verb: Remove selected modules, then register new
            % mappings for those modules.

            mrstPath('remove'  , mods{1 : 2 : end});
            mrstPath('register', mods{:});

         case 'reset',
            % Convenience verb: Clear existing list, then register new
            % modules.

            mrstPath clear
            mrstPath('register', mods{:});

         case { 'query', 'search' },
            % Look for a particular module or set of modules.  The 'mods'
            % list is interpreted as a list of modules to search for.

            [pth, notfound] = search_modules(CACHE, mods);

            if nargout > 0,
               lst = repmat({ '' }, [1, numel(mods)]);
               [lst{~notfound}] = pth{:};

               if (nargout == 1) && (numel(lst) > 1),

                  varargout{1} = lst;

               elseif nargout == numel(mods),

                  [varargout{1:numel(mods)}] = lst{:};

               end

            else

               map = [ reshape(mods(~notfound), [], 1) , ...
                       reshape(pth,             [], 1) ; ...
                       reshape(mods(notfound),  [], 1) , ...
                       repmat({ 'Not Found' }, [sum(notfound), 1]) ];

               fprintf('Module mapping results:\n');
               print_list(map);
            end

         otherwise

            error(msgid('Verb:Unknown'), ...
                 ['Unknown command verb ''%s''.  Must be one of ', ...
                  '''register'', ''clear'', ''list'', ''remove'', ', ...
                  '''reset'', or ''query''.'], cmd);

      end

   elseif nargout == 0,

      mrstPath list

   else

      error(msgid('Syntax:Error'), ...
           ['Call syntax is\n\t', ...
            mfilename, ' <command> [module list]  or\n\t', ...
            'map = ', mfilename]);

   end
end

%--------------------------------------------------------------------------

function cache = register_root(cache, mods)
   mods = mods(cellfun(@isdir, mods));
   m2d  = @(r, m) cellfun(@(x) fullfile(r, x), m, 'UniformOutput', false);

   for mroot = reshape(mods, 1, []),
      m = dir(mroot{1});
      m = { m([ m.isdir ]).name };
      m = m(~ strncmp('.', m, 1));

      rlist = reshape([m ; m2d(mroot{1}, m)], 1, []);
      cache = register_modules(cache, rlist);
   end
end

%--------------------------------------------------------------------------

function cache = register_modules(cache, mods)
   assert (mod(numel(mods), 2) == 0, ...
           'Register: Must be list of key/value pairs');

   mods = reshape(mods, 2, []) .';

   dsrc = { pwd };
   root = ROOTDIR;

   % The FULLFILE(x, '.') construct is to transparently handle presence or
   % absence of a terminating FILESEP on either of the directory strings.
   %
   cmp         = match_algorithm;
   canonical   = @(d) fileparts(fullfile(d, '.'));
   dir_is_same = @(d1, d2) cmp(canonical(d1), canonical(d2));

   if ~ dir_is_same(dsrc{1}, root),
      dsrc = [ dsrc , { root } ];
   end
   dsrc = [ dsrc , { '' } ];

   [I, J] = blockDiagIndex(size(mods, 1), numel(dsrc));

   t = cellfun(@(r, m) fullfile(r, m), ...
               reshape(dsrc(J), [],  1), ...
               reshape(mods(I, 2), [], 1), 'UniformOutput', false);
   t = reshape(t, size(mods, 1), []);

   i = cellfun(@isdir, t);
   e = any(i, 2);

   if ~ all(e),
      if mrstVerbose,
         nonext = mods(~ e, :);

         if size(nonext, 1) > 1,
            [pl1, pl2, pl3] = deal('s', 'ies', ''  );
         else
            [pl1, pl2, pl3] = deal('' , 'y'  , 'es');
         end

         fprintf('The following module%s could not be registered.\n', pl1);
         fprintf('Reason: The director%s do%s not exist.\n', pl2, pl3);

         print_list(nonext);
      end

      mods = mods(e, :);
      t    = t(e, :);
      i    = i(e, :);
   end

   % 'i' is a LOGICAL matrix.  Consequently, MAX(i, [], 2) is the *first*
   % TRUE entry on each row.  Its location (within the row) is exactly the
   % column from which we must extract the pathname expansion of
   % 'mods(:,1)'.
   %
   [col, col] = max(i, [], 2);   %#ok
   pathix     = sub2ind(size(t), (1 : size(t,1)) .', reshape(col, [], 1));
   mods(:, 2) = cellfun(canonical, t(pathix), 'UniformOutput', false);

   map = cache_to_map(cache);
   [ii, jj] = look_for(map(:,1), mods(:,1));
   if ~isempty(ii),
      map(ii, 2) = mods( jj, 2);
      mods       = mods(~jj, :);
   end

   cache = map_to_cache([ map ; mods ]);
end

%--------------------------------------------------------------------------

function cache = clear_modules
   cache = cell([0, 2]);
end

%--------------------------------------------------------------------------

function tf = is_empty(cache)
   tf = isempty(cache_to_map(cache));
end

%--------------------------------------------------------------------------

function list_modules(cache)
   map = cache_to_map(cache);

   if isempty(map),
      fprintf('No modules registered\n');
   else
      fprintf('Currently registered modules\n');
      print_list(map);
   end
end

%--------------------------------------------------------------------------

function cache = remove_modules(cache, mods)
   map = cache_to_map(cache);
   i   = look_for(map(:,1), mods);

   pick    = false([size(map,1), 1]);
   pick(i) = true;

   cache = map_to_cache(map(~pick, :));
end

%--------------------------------------------------------------------------

function [pth, notfound] = search_modules(cache, mods)
   map    = cache_to_map(cache);
   [i, j] = look_for(map(:,1), mods);

   pth = map(i, 2);

   notfound     = false([numel(mods), 1]);
   notfound(~j) = true;
end

%--------------------------------------------------------------------------

function map = cache_to_map(map)
   % Identity, no-op
end

%--------------------------------------------------------------------------

function cache = map_to_cache(cache)
   % Identity, no-op
end

%--------------------------------------------------------------------------

function [match_i, match_j] = look_for(reg, mods)
   [i, j] = blockDiagIndex(numel(reg), numel(mods));

   matches = match(reshape(reg(i), [], 1), reshape(mods(j), [], 1));

   match_i = i(matches);
   match_j = false([numel(mods), 1]);
   match_j(j(matches)) = true;
end

%--------------------------------------------------------------------------

function tf = match(m1, m2)
   alg = match_algorithm;
   tf  = alg(m1, m2);
end

%--------------------------------------------------------------------------

function alg = match_algorithm
   if ispc,
      % Match module names irrespective of case on PC/WIN.
      alg = @strcmpi;
   else
      % Consider case when matching module names on Unix-like systems.
      alg = @strcmp;
   end
end

%--------------------------------------------------------------------------

function print_list(map)
   nchar = max(cellfun('prodofsize', map(:,1)));

   for k = 1 : size(map, 1),
      fprintf('  * %-*s -> %s\n', nchar, map{k,1}, map{k,2});
   end
end
