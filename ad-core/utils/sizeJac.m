function sz = sizeJac(jac, dim)
    if ~isa(jac, 'DiagonalJacobian')
        sz = size(jac, dim);
        return
    end
    if nargin < 2
        dim = ':';
    end
    sz(1) = size(jac.diagonal, 1);
    sz(2) = prod(jac.dim);
    sz    = sz(dim);
end