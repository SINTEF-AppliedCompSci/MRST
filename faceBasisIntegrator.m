function [x, cellNo, faceNo, W] = faceBasisIntegrator(model)

    G = model.G;
    nDof = model.basis.nDof;
    degree = model.degree;
    [x, w, nq, ii, jj, cellNo, faceNo] = makeFaceIntegrator(G, (1:G.cells.num)', degree);
    
    x = (x - G.cells.centroids(cellNo,:))./G.cells.diameters(cellNo);
    
    [ii, jj] = blockDiagIndex(ones(G.cells.num*nDof, 1), repmat(nq, nDof,1));
    W = sparse(ii, jj, repmat(w, nDof, 1));
    
end