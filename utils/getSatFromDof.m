function s = getSatFromDof(x, cell, dof, disc)

    G       = disc.G;
    psi     = disc.basis.psi;
%     psi_lim = disc.basis.psi_lim;
    nDof    = disc.basis.nDof;
    
    ii = sum((1:G.cells.num)' == cell',2);
    s = 0;
    for dofNo = 1:nDof
        
        ind = (1:nDof:G.cells.num*nDof)' + dofNo - 1;
        ix = rldecode(ind, ii, 1);
        s = s + dof(ix).*psi{dofNo}(x);%.*(~limflag(cell));
        
%         if disc.degree > 0 && dofNo <= polyDim(1, G.griddim)
%             s = s + dof(ix).*psi_lim{dofNo}(x).*limflag(cell);
%         end
        
    end
    
%     if model.degree > 0
%         psi_lim = model.basis.psi_lim;
%         for dofNo = 1:3
            
            

end