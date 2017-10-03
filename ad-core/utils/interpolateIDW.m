function v = interpolateIDW(x, f, xq, order)
    % Number of sample points and dimension of space
    [ns, nx] = size(x);
    % Dimension of output space
    nf = size(f, 2);
    % Number of query points
    nq = size(xq, 1);
    
    assert(size(f, 1) == ns);
    assert(size(xq, 2) == nx);
    
    v = nan(nq, nf);
    for i = 1:nq
        d = bsxfun(@minus, x, xq(i, :));
        dist = max(sum(abs(d).^order, 2), 1e-20);
        wi = 1./dist;
        W = repmat(wi, 1, nf);
        for j = 1:nf
            disabled = isnan(f(:, j));
            W(disabled, j) = 0;
            f(disabled, j) = 0;
            s = sum(W(:, j));
            if s < 1e-18
                W(:, j) = repmat(1/ns, ns, 1);
            else
                W(:, j) = W(:, j)./s;
            end
        end
        v(i, :) = sum(W.*f, 1);
    end
end
