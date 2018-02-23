function state = getCellSaturation(model, state)
    
    disc = model.disc;
    G    = model.G;
    
    
    sfun = @(x,c, phNo) getSatFromDof(x, c, state.sdof(:,phNo), disc);
    [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, (1:G.cells.num)', disc.degree);
%     if disc.degree > 1
%         x = (x - G.cells.centroids(cellNo))./(G.cells.diameters(cellNo)/(2*sqrt(G.griddim)));
%     end
    W = sparse(ii, jj, w);

    x = disc.transformCoords(x, cellNo);
    
    nPh = nnz(model.getActivePhases);
    s = zeros(G.cells.num, nPh);
    for phNo = 1:nPh
        s(:,phNo) = (W*sfun(x, cellNo, phNo))./G.cells.volumes;
%         s(:,phNo) = disc.cellInt(@(x, c, psi) sfun(x, c, phNo), (1:G.cells.num)')./G.cells.volumes;
%         (W*sfun(x, cellNo, phNo))./G.cells.volumes;
    end
    
    state.s = s;
    
end