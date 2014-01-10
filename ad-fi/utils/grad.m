function g = grad(N)
persistent C
if isempty(C)
    if nargin == 0
        N  = getGridNeighbors;
    end
    nc = max(max(N));
    nf = size(N,1);
    C = sparse((1:nf)'*[1 1], N, ones(nf,1)*[1 -1], nf, nc);
end
g = C;
end
