function buildMexOperators(rebuild, varargin)
    opt = struct('names', {{}});
    opt = merge_options(opt, varargin{:});
    if nargin == 0
        rebuild = false;
    end
    
    filenames = get_names(opt);
    
    if rebuild
        delete_compiled_files(filenames);
    end
    build_files(filenames);
end

function build_files(names)
    for i = 1:numel(names)
        n = names{i};
        try
            timer = tic();
            fprintf('Building MEX file %s\n', n)
            eval(n);
            fprintf('%s built in %2.2fs\n', n, toc(timer));
        catch
            fprintf('Compilation FAILED for %s. Run %s() to manually recompile and see errors.\n', n, n);
        end
    end
end

function names = get_names(opt)
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
end

function delete_compiled_files(names)
    pth = fullfile(mrstPath('ad-core'), 'backends',...
                            'diagonal', 'operators', 'mex');

    ext = mexext();
    for i = 1:numel(names)
        n = names{i};
        fp = fullfile(pth, sprintf('%s.%s', n, ext));
        if exist(fp, 'file')
            fprintf('Removing MEX file %s\n', n)
            clear(n);
            delete(fp);
        end
    end
    rehash
end