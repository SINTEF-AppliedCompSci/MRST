mrstModule add vem vemmech

%%

n = 2;
G = computeGeometry(cartGrid([n,n], [2*n,2]));
G.nodes.coords = G.nodes.coords;
G = computeVEMGeometry(G);

%%

degree = 4;
basis = dgBasis(G, degree, 'legendre');

% basis = DGBasisFunctions(G, degree);

k = basis.k;
p = cell(size(k,1),1);
for pNo = 1:size(k,1)
    p{pNo} = Polynomial(k(pNo,:), 1);
end

cells = (1:G.cells.num)';
[x, w, nq, cellNo, pointNo] = makeCellIntegrator(G, cells, degree);
x = x - G.cells.centroids(cellNo,:);

% xhat = (x - G.cells.centroids(cellNo,:))./(G.cells.diameters(cellNo)/(2*sqrt(G.griddim)));

W = sparse(cellNo, pointNo, w);
ci = cellfun(@(p) W*p(x), p, 'unif', false);
% ci = cellfun(@(psi) W*psi(x, cellNo), basis.psi, 'unif', false);

% ci = W*psi{1}(x-G.cells.centroids(cellNo, :));