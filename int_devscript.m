mrstModule add vem vemmech

%%

n = 2;
G = computeGeometry(cartGrid([n,n], [1,1]));
G.nodes.coords = G.nodes.coords;
G = computeVEMGeometry(G);

%%

degree = 1;
basis = dgBasis(degree, G.griddim, 'legendre');

cells = (1:G.cells.num)';
[x, w, nq, cellNo, pointNo] = makeCellIntegrator(G, cells, 2, 'tri');
W = sparse(cellNo, pointNo, w);

x = x - G.cells.centroids(cellNo,:);
ci = cellfun(@(psi) W*psi(x), basis.psi, 'unif', false);

% ci = W*psi{1}(x-G.cells.centroids(cellNo, :));