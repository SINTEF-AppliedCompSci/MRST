function [x, cellNo, W] = cellBasisIntegrator(disc)

    G      = disc.G;
    nDof   = disc.basis.nDof;
    degree = disc.degree*G.griddim;
%     degree = 2;
    [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, (1:G.cells.num)', degree, 'tri');
    
    x = (x - G.cells.centroids(cellNo,:))./(G.cells.diameters(cellNo)/(2*sqrt(G.griddim)));
    
    [ii, jj] = blockDiagIndex(ones(G.cells.num*nDof, 1), repmat(nq, nDof,1));
    W = sparse(ii, jj, repmat(w, nDof, 1));
    
end