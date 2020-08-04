function alphas = computealphas(eta, nz)
    alphas = NaN(nz, 1);
    options = optimset('fzero');
    options = optimset(options, 'tolfun', 1e-16, 'tolx', 1e-16);
    for i = 1 : nz
        f = @(x) (tan(x) - 2*eta*(x + (i - 1)*pi));
        z = fzero(f, [eps, pi/2 - eps], options);
        alphas(i) = (i - 1)*pi + z;
    end
    
end
