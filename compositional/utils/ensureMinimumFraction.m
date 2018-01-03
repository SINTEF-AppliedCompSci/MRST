function z = ensureMinimumFraction(z, tol)
% Set a minimum value on a composition matrix
    if nargin == 1
        tol = 1e-12;
    end
    z = max(z, tol);
    z = bsxfun(@rdivide, z, sum(z, 2));
end