function buildMexExtensions(rebuild, names, varargin)
% (Re)build a set of mex extensions located in a specific folder

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
    assert(islogical(rebuild))
    [filenames, paths] = get_names(names, varargin{:});
    
    if rebuild
        delete_compiled_files(filenames, paths);
    end
    build_files(filenames, paths);
end

function build_files(names, paths)
    for i = 1:numel(names)
        n = names{i};
        p = paths{i};
        if exist(p, 'file')
            fprintf('%s is already compiled -> OK!\n', names{i});
        else
            try
                timer = tic();
                fprintf('Building MEX file %s...\n', n)
                eval(n);
                fprintf('Extension built in %2.2fs -> OK!\n', toc(timer));
            catch
                fprintf('Compilation FAILED for %s. Run %s() to manually recompile and see errors.\n', n, n);
            end
        end
    end
end

function [names, paths] = get_names(names, varargin)
    if ~iscell(names)
        names = {names};
    end
    isOctave = mrstPlatform('octave');
    n = numel(names);
    ext = mexext();
    if nargin == 1
        % Recieved full paths
        paths = names;
        names = cell(n, 1);
        for i = 1:n
            [~, names{i}, ~] = fileparts(paths{i});
        end
    else
        paths = cell(n, 1);
        % Recieved named files + folders
        pth = varargin{1};
        for i = 1:n
            n = names{i};
            if isOctave
                if is_octfile(sprintf('%s.%s', n, 'cpp'))
                    fe = 'oct';
                else
                    fe = ext;
                end
            else
                fe = ext;
            end
            paths{i} = fullfile(pth, sprintf('%s.%s', n, fe));
        end
    end
end

function delete_compiled_files(names, paths)
    for i = 1:numel(names)
        n = names{i};
        fp = paths{i};
        if exist(fp, 'file')
            fprintf('Removing MEX file %s ', n)
            clear(n);
            if mrstPlatform('octave')
                autoload(n, pwd(), 'remove');
            end
            delete(fp);
            if exist(fp, 'file')
                fprintf('-> FAILURE. Unable to delete. Is it loaded in another session?\n');
            else
                fprintf('-> OK!\n');
            end
        end
    end
    rehash
end

function ok = is_octfile(filename)
    tmp = fileread(filename);
    ok = ~isempty(strfind(tmp, 'DEFUN_DLD')); %#ok
end
