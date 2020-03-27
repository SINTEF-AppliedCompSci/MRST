function buildMexOperators(varargin)
% Build MEX operators for automatic differentiation
    if mod(nargin, 2) == 1
        rebuild = varargin{1};
        varargin = varargin(2:end);
    else
        rebuild = false;
    end
    opt = struct('names', {{}});
    opt = merge_options(opt, varargin{:});

    
    [filenames, paths] = get_names(opt);
    
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

function [names, paths] = get_names(opt)
    if isempty(opt.names)
        names = {'mexDiscreteDivergenceJac', ...
                 'mexDiscreteDivergenceVal', ...
                 'mexDiagonalSparse',        ...
                 'mexFaceAverageDiagonalJac', ...
                 'mexSinglePointUpwindVal', ...
                 'mexSinglePointUpwindDiagonalJac', ...
                 'mexTwoPointGradientDiagonalJac', ...
                 'mexTwoPointGradientVal', ...
                 'mexDiagMult', ...
                 'mexDiagProductMult'};
    else
        names = opt.names;
        if ~iscell(names)
            names = {names};
        end
    end
    n = numel(names);
    pth = fullfile(mrstPath('ad-core'), 'backends',...
                            'diagonal', 'operators', 'mex');

    ext = mexext();
    paths = cell(n, 1);
    for i = 1:n
        n = names{i};
        paths{i} = fullfile(pth, sprintf('%s.%s', n, ext));
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