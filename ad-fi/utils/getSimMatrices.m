function [N, T] = getSimMatrices(G, T)
N  = double(G.faces.neighbors);
intInx = (prod(N,2)~=0);
N  = N(intInx, :);
clear getGridNeighbors
getGridNeighbors(N);
% half-trans -> trans and reduce to interior
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T  = T(intInx);



