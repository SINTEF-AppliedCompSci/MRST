function buildMexExtensions(rebuild, names, varargin)
% (Re)build a set of mex extensions located in a specific folder
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
            paths{i} = fullfile(pth, sprintf('%s.%s', n, ext));
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