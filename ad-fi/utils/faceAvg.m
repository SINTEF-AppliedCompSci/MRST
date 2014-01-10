function avg = faceAvg(expr, N)
persistent M
if isempty(M)
    if nargin < 2
        N = getGridNeighbors;
    end
    nc = max(max(N));
    nf = size(N,1);
    M  = sparse((1:nf)'*[1 1], N, .5*ones(nf,2), nf, nc);
end
avg = M*expr;
end
