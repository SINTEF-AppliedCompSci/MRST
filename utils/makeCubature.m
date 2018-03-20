function cubature  = makeCubature(disc, type)

    G = disc.G;
    degree = disc.degree + 1;
    issurf = strcmp(type, 'surface');
    dim = G.griddim - issurf;
    basis = dgBasis(dim, degree, 'legendre');
    nDof = basis.nDof;
    psi  = basis.psi;
    
    x = getQuadraturePoints(degree+1, dim);
    
    P = reshape(cell2mat(cellfun(@(p) p(x), psi, 'unif', false)), nDof, nDof)';
    
    [xq, w, ~, ii, jj, cellNo, faceNo] = makeCellIntegrator(G, (1:G.cells.num)', degree, type);
    
    W = sparse(ii, jj, w);
    
    xq = disc.transformCoords(xq, cellNo);
    
    I = cellfun(@(p) W*p(xq), psi, 'unif', false);
    rhs = zeros(nDof, G.cells.num);
    for dofNo = 1:nDof
        rhs(dofNo,  :) = I{dofNo};
    end
    
    w = reshape(P\rhs, [], 1);
    
    nq = nDof;
    cellNo = reshape(repmat((1:G.cells.num), nq, 1), [], 1);
    x = disc.transformCoords(repmat(x, G.cells.num, 1), cellNo, true);
    
%     x = repmat(x, G.cells.num, 1) + rldecode(G.cells.centroids, nq, 1);

    cubature = struct('x', x, 'w', w, 'nq', nq, 'cellNo', cellNo, 'faceNo', faceNo);
    
end