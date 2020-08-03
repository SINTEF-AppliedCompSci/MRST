function [L, x, y, Z_L, Z_V, p, T] = getFlashTable(p_range, T_range, z, EOS)
% Utility for getting meshgrid plots of flash results (e.g. for phase diagram)
    [p, T] = meshgrid(p_range, T_range);
    np = numel(p_range);
    nT = numel(T_range);
    
    n = numel(p);
    ncomp = numel(z);
    Z = cell(1, ncomp);
    for i = 1:ncomp
        if iscell(z)
            zi = z{i};
        else
            zi = z(i);
        end
        Z{i} = repmat(zi, n, 1);
    end
    
    [L, x, y, Z_L, Z_V] = standaloneFlash(p(:), T(:), Z, EOS);
    
    fix = @(x) reshape(x, nT, np)';
    L = fix(L);
    Z_L = fix(Z_L);
    Z_V = fix(Z_V);
    for i = 1:numel(x)
        x{i} = fix(x{i});
        y{i} = fix(y{i});
    end
    T = fix(T);
    p = fix(p);
end
