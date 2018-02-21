function state = getCellSaturation(model, state)
    
    disc = model.disc;
    G    = model.G;
    
    
    sfun = @(x,c, phNo) getSatFromDof(x, c, state.sdof(:,phNo), disc);
    [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, (1:G.cells.num)', disc.degree, 'tri');
    W = sparse(ii, jj, w);
    
    nPh = nnz(model.getActivePhases);
    s = zeros(G.cells.num, nPh);
    for phNo = 1:nPh
        s(:,phNo) = (W*sfun(x, cellNo, phNo))./G.cells.volumes;
    end
    
    state.s = s;
    
end