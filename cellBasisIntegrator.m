function [x, cellNo, W] = cellBasisIntegrator(model)

    G = model.G;
    nDof = model.basis.nDof;
    degree = model.degree*G.griddim;
    [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, (1:G.cells.num)', degree, 'tri');
    
    x = (x - G.cells.centroids(cellNo,:))./G.cells.diameters(cellNo);
    
    [ii, jj] = blockDiagIndex(ones(G.cells.num*nDof, 1), repmat(nq, nDof,1));
    W = sparse(ii, jj, repmat(w, nDof, 1));
    
end