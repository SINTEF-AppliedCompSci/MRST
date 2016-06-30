function varargout = mrstExamples(varargin)
%PLOTGRID plots exterior grid faces to current axes.
%
% SYNOPSIS:
%       mrstExamples('ad-blackoil')
%       mrstExamples()
%       exList = mrstExamples('ad-blackoil', 'diagnostics')
%
% PARAMETERS:
%   input - Any number of strings that are the names of registered mrst
%           modules (see mrstPath / mrstModule).
%
% RETURNS:
%   ex    - A cell array of same length of the number of input arguments,
%           where each entry is a list of the paths to the examples
%           corresponding to that module. For instance, if called as 
%                ex = mrstExamples('firstmodule', 'secondmodule')
%           ex{1} will be a cell array of the examples corresponding to the
%           module "firstmodule" and so on.
%
% NOTES:
%   If no output arguments are given, the routine instead prints clickable
%   editor links of all examples to the command window.
%
%   Examples are defined as all *.m files found in the examples directory
%   of a module. Subdirectories of the examples directory are also included
%   in the search.
%
% SEE ALSO:
%   mrstModule, mrstPath

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
    if nargin == 0
        varargin = {'core'};
    end
    doPrint = nargout == 0;

    examplePaths = cell(nargin, 1);
    for i = 1:numel(varargin)
        name = varargin{i};

        if isempty(name) || strcmpi(name, 'core')
            pth = ROOTDIR();
        else
            pth = mrstPath(name);
            if isempty(pth) && doPrint
                fprintf('Module "%s" is not known to MRST\n', name);
                continue;
            end
        end
        dirs = {fullfile(pth, 'examples')};
        ix = 1;
        ndir = numel(dirs);
        allfiles = {};
        allnames = {};
        while ix <= ndir
            [d, names, paths] = getSub(dirs{ix});
            dirs = [dirs, d{:}];                                           %#ok
            allfiles = [allfiles, paths{:}];                               %#ok
            allnames = [allnames, names{:}];                               %#ok
            ix = ix + 1;
            ndir = numel(dirs);
        end
        nex = numel(allnames);
        if doPrint
            fprintf('Module "%s" has %d example', name, nex);
            if nex == 1
                fprintf(':\n');
            elseif nex > 1
                fprintf('s:\n');
            else
                fprintf('s.\n');
            end
            for j = 1:nex
                % Split by path + examples to get base folder path
                spl = strsplit(allfiles{j}, fullfile(pth, ['examples', filesep]));
                fprintf('    <a href="matlab: edit %s">%s\n</a>',...
                                                allfiles{j}, spl{end});
            end
        end
        examplePaths{i} = allfiles;
    end

    if ~doPrint
        varargout{1} = examplePaths;
    end
end

function [subdir, mnames, mpaths] = getSub(path)
    cand = dir(path);
    % Filter very short names, hidden/up directory files and Contents.m
    % files.
    filter = @(x) ~strcmpi(x.name(1), '.') && ...
                  ~strcmpi(x.name, 'contents.m') && ...
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
    mnames = arrayfun(@(x) x.name, mfiles, 'UniformOutput', false);
end