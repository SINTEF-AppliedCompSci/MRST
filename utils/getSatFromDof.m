function s = getSatFromDof(x, cell, dof, model)

    G    = model.G;
    psi  = model.basis.psi;
    nDof = model.basis.nDof;
    
    ii = sum((1:G.cells.num)' == cell',2);
    s = 0;
    for dofNo = 1:nDof
        ind = (1:nDof:G.cells.num*nDof) + dofNo - 1;
        dofloc = rldecode(dof(ind), ii, 1);
        s = s + dofloc.*psi{dofNo}(x);
    end
    
end