function z = ensureMinimumFraction(z, tol)
    if nargin == 1
        tol = 1e-8;
    end
    sz = 0;
    for i = 1:numel(z)
        z{i} = max(z{i}, tol);
        sz = sz + z{i};
    end
    z = cellfun(@(x) x./sz, z, 'unif', false);

end