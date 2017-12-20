function [x, D] = diagMult(v, M, D)
% Internal function for diagonal multiplication in AD code
    if ~any(v)
        sz = size(M);
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
