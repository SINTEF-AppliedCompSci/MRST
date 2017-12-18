function [x, D] = diagMult(v, M, D)
    if ~any(v)
%         x = 0*M;
        sz = size(M);
%         x = ZeroJacobian(sz);
        x = sparse([],[],[],sz(1), sz(2));
    elseif nnz(M) == 0
        x = M;
    elseif isscalar(v)
        x = v*M;
    else
        if isempty(D)
            n = numel(v);
            ix = (1:n)';
            D = sparse(ix, ix, v, n, n);
        end
        x = D*M;
    end
end
