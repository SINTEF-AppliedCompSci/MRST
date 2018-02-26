function s = getSatFromDof(x, cells, dof, disc)

    G       = disc.G;
    psi     = disc.basis.psi;
    nDof    = disc.basis.nDof;
    
    ii = sum((1:G.cells.num)' == cells',2);
    s = 0;
    for dofNo = 1:nDof
        ind = (1:nDof:G.cells.num*nDof)' + dofNo - 1;
        ix = rldecode(ind, ii, 1);
        s = s + dof(ix).*psi{dofNo}(x);
    end
    
end