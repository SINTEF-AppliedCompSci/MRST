function ok = buildMexOperators(rebuild, varargin)
    opt = struct('names', {{}});
    opt = merge_options(opt, varargin{:});
    if nargin == 0
        rebuild = false;
    end
    
    if rebuild
        delete_compiled_files(opt);
    end
    ok = build_files();
end

function ok = build_files()
    G = cartGrid([2, 2, 2]);
    G = computeGeometry(G);
    try
        res = testMexDiagonalOperators(G);
        ok = true;
    catch
        ok = false;
    end
end

function delete_compiled_files(opt)
    if isempty(opt.names)
        names = {'mexDiscreteDivergenceJac', ...
                 'mexDiscreteDivergenceVal', ...
                 'mexDiagonalSparse',        ...
                 'mexFaceAverageDiagonalJac', ...
                 'mexSinglePointUpwindVal', ...
                 'mexSinglePointUpwindDiagonalJac', ...
                 'mexTwoPointGradientDiagonalJac', ...
                 'mexTwoPointGradientVal'};
    else
        names = opt.names;
        if ~iscell(names)
            names = {names};
        end
    end
    pth = fullfile(mrstPath('ad-core'), 'backends',...
                            'diagonal', 'operators', 'mex');

    ext = mexext();
    for i = 1:numel(names)
        fn = names{i};
        fp = fullfile(pth, sprintf('%s.%s', fn, ext));
        if exist(fp, 'file')
            clear(fn);
            delete(fp);
        end
    end
    rehash
end