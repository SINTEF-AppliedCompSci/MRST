function ok = buildMexOperators(rebuild)
    if nargin == 0
        rebuild = false;
    end
    pth = fullfile(mrstPath('ad-core'), 'backends',...
                            'diagonal', 'operators', 'mex');
    ls(pth)
    ext = mexext();
    names = {'mexDiscreteDivergenceJac', ...
             'mexDiscreteDivergenceVal', ...
             'mexDiagonalSparse',        ...
             'mexFaceAverageDiagonalJac', ...
             'mexSinglePointUpwindVal', ...
             'mexSinglePointUpwindDiagonalJac', ...
             'mexTwoPointGradientDiagonalJac', ...
             'mexTwoPointGradientVal'};
    if rebuild
        for i = 1:numel(names)
            fn = names{i};
            fp = fullfile(pth, sprintf('%s.%s', fn, ext));
            if exist(fp, 'file')
                clear(fn);
                delete(fp);
            end
        end
    end
    G = cartGrid([2, 2, 2]);
    G = computeGeometry(G);
    try
        res = testMexDiagonalOperators(G);
        ok = true;
    catch
        ok = false;
    end
end
