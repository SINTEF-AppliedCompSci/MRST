function varargout = mrstExamples(varargin)
%Discover Example M-Files Pertaining to One or More MRST Modules
%
% SYNOPSIS:
%            mrstExamples [module list]
%            mrstExamples all
%   exList = mrstExamples(...)
%
% PARAMETERS:
%   varargin -
%           Sequence of strings that are treated as names of registered
%           MRST modules (see functions `mrstPath` and `mrstModule`).  The
%           special module name 'core', although not a module in a strict
%           sense, represents those examples that are available in the base
%           MRST package--i.e., without activating any modules at all.
%
%           Alternatively, the single string 'all' can be given to list the
%           examples in all registered modules including the base MRST
%           package.
%
% RETURNS:
%   exList - A cell array of same length as the number of input arguments,
%            where each entry itself is a list (cell array of strings) of
%            the paths to the examples of the corresponding module.  For
%            instance, if called as ::
%
%                exList = mrstExamples('ad-blackoil', 'diagnostics')
%
%            then `exList{1}` will be a cell array of the examples relating
%            to module `ad-blackoil` while `exList{2}` is a cell array of
%            the examples of module `diagnostics`.
%
%            If called using the special string 'all', then `exList`
%            contains one element for each known module and one additional
%            element, specifically `exList{1}`, corresponding to the `core`
%            (i.e., base) MRST package.
%
% NOTES:
%   If no output argument is given, the routine instead prints clickable
%   editor links of all examples to the command window.
%
%   Examples are defined as all M files found in the `examples` directory
%   of a module--excluding any M files in a subdirectory called `utils`.
%   Subdirectories of the `examples` directory other than `utils` are
%   searched recursively.
%
% SEE ALSO:
%   `mrstModule`, `mrstPath`.

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

    varargin = parse_arguments(varargin{:});

    doPrint      = nargout == 0;
    npaths       = numel(varargin);
    examplePaths = cell([npaths, 1]);

    for i = 1:npaths
        name = varargin{i};
        pth  = get_module_path(name);

        if isempty(pth) && doPrint
           fprintf('Module ''%s'' is not known to MRST (skipped)\n', name);
           continue;
        end

        examplePaths{i} = module_examples(pth);

        if doPrint
           print_examples(name, examplePaths{i});
        end
    end

    if ~doPrint
        varargout{1} = examplePaths;
    end
end

%--------------------------------------------------------------------------

function varargin = parse_arguments(varargin)
   if nargin == 0
      varargin = {'core'};
   end

   if any(strcmpi(varargin, 'all'))
      % List all modules, including core
      assert (nargin == 1, ...
             ['Module designation ''all'' must be ', ...
              'only argument if present']);

      % Put core first in the module list
      varargin = [{ 'core' }; reshape(sort(mrstPath()), [], 1)];
   end
end

%--------------------------------------------------------------------------

function pth = get_module_path(name)
   if isempty(name) || strcmpi(name, 'core')
      pth = ROOTDIR();
   else
      pth = mrstPath(name);
   end
end

%--------------------------------------------------------------------------

function allfiles = module_examples(mroot)
   dirs = {fullfile(mroot, 'examples')};
   ix   = 1;
   ndir = numel(dirs);

   allfiles = {};

   while ix <= ndir
      [d, paths] = getSub(dirs{ix});

      dirs     = [dirs    , reshape(d    , 1, [])];                    %#ok
      allfiles = [allfiles, reshape(paths, 1, [])];                    %#ok

      ix   = ix + 1;
      ndir = numel(dirs);
   end
end

%--------------------------------------------------------------------------

function print_examples(module, allfiles)
   nex = numel(allfiles);

   fprintf('Module "%s" has %d example', module, nex);
   if nex == 1
      fprintf(':\n');
   elseif nex > 1
      fprintf('s:\n');
   else
      % nex == 0 -- No examples in this module.
      fprintf('s.\n');
      return
   end

   % Split pathnames on the token sequence
   %
   %    [filesep, 'examples', filesep]
   %
   % and derive the example file names relative to the 'examples' directory
   % ('ex') for purpose of printing shortened, clickable EDIT links (URLs)
   % in the Command Window.  We assume that there is but one occurrence of
   % the above token sequence within each element of 'allfiles' so the
   % relative name is the second component of the split pathname.

   options = { 'split' };
   if ispc
      % Case insensitive pathname matching on MS Windows.
      options = [ options, { 'ignorecase' } ];
   end

   sep = regexptranslate('escape', [filesep, 'examples', filesep]);

   ex  = cellfun(@(p) p(2), regexp(allfiles, sep, options{:}));
   if mrstPlatform('richtext')
      prt = [ reshape(allfiles, 1, []) ; ...
              reshape(ex,       1, []) ];
      fprintf('    <a href="matlab: edit ''%s''">%s</a>\n', prt{:});
   else
      prt = reshape(ex, 1, []);
      fprintf('    %s\n', prt{:});
   end
end

%--------------------------------------------------------------------------

function [subdir, mpaths] = getSub(path)
    cand = dir(path);
    % Filter very short names, hidden/up directory files and Contents.m
    % files.
    filter = @(x) ~strcmpi(x.name(1), '.') && ...
                  ~strcmpi(x.name, 'contents.m') && ...
                  ~(strcmpi(x.name, 'utils') && x.isdir) &&...
                  (numel(x.name) > 2 || x.isdir);
    ok = arrayfun(filter, cand);
    cand = cand(ok);
    % Get subdirectories
    subdir = cand([cand.isdir]);
    % Get files with .m extension
    ismfile = arrayfun(@(x) strcmpi(x.name(end-1:end), '.m'), cand);
    mfiles = cand(ismfile);
    
    subdir = arrayfun(@(x) fullfile(path, x.name), subdir, 'UniformOutput', false);
    mpaths = arrayfun(@(x) fullfile(path, x.name), mfiles, 'UniformOutput', false);
end
