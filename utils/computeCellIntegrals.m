function ci = computeCellIntegrals(G, degree)

    [k, nDof] = dgBasis(degree, G.griddim);
%     C = zeros(G.cells.num, nDof*nDof);
    intFun = makeCellIntegrator(G, (1:G.cells.num)', degree*2, 'tri');
    ci = zeros(nDof^2*G.cells.num,1);
    for i = 1:nDof
        psii = Polynomial(k(i,:));
        for j = 1:nDof
            psij = Polynomial(k(j,:));
            psii_psij = psii*psij;
            
            ix = (1:nDof^2:G.cells.num*nDof^2) + (i-1)*nDof + j-1;
            
            ci(ix) = intFun(psii_psij);
%             ci = [ci; intFun(pi_pj)];
            
        end
    end
    
    [ii, jj] = blockDiagIndex(repmat(nDof, G.cells.num, 1));
    ci = sparse(jj, ii, ci);
%     ci = 

end