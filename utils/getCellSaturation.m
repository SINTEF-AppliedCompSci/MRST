function state = getCellSaturation(model, state)
    
    disc = model.disc;
    G    = model.G;
    
    
    sfun = @(x,c, phNo) getSatFromDof(x, c, state.sdof(:,phNo), disc);
    [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, (1:G.cells.num)', disc.degree);

    W = sparse(ii, jj, w);

    x = disc.transformCoords(x, cellNo);
    
    nPh = nnz(model.getActivePhases);
    s = zeros(G.cells.num, nPh);
    for phNo = 1:nPh
        s(:,phNo) = (W*sfun(x, cellNo, phNo))./G.cells.volumes;
    end
    
    s0 = s;
    
    nDof = model.disc.basis.nDof;
    ix = 1:nDof:G.cells.num*nDof;
    
    
    s(:,1) = state.sdof(ix);
    s(:,2) = 1 - s(:,1);
    
    state.s = s;
    
end