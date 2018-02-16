function s = getSatFromDof(x, cell, dof, limflag, model)

    G       = model.G;
    psi     = model.basis.psi;
    psi_lim = model.basis.psi_lim;
    nDof    = model.basis.nDof;
    
    ii = sum((1:G.cells.num)' == cell',2);
    s = 0;
    for dofNo = 1:nDof
        
        ind = (1:nDof:G.cells.num*nDof)' + dofNo - 1;
        ix = rldecode(ind, ii, 1);
        s = s + dof(ix).*psi{dofNo}(x).*(~limflag(cell));
        
        if model.degree > 0 && dofNo <= polyDim(1, G.griddim)
            s = s + dof(ix).*psi_lim{dofNo}(x).*limflag(cell);
        end
        
    end
    
%     if model.degree > 0
%         psi_lim = model.basis.psi_lim;
%         for dofNo = 1:3
            
            

end