function s = getSatFromDof(x, cell, dof, psi, G)


    ii = sum((1:G.cells.num)' == cell',2);
    nDof = 3;
    w = 1;
%     dof = rldecode(dof, cell, 1);
%     s = @(x) 0*x(:,1);
    s = 0;
%     psi = {@(x) w*prod(x.^[0,0],2), @(x) w*prod(x.^[1,0],2), @(x) w*prod(x.^[0,1],2)};
    for dofNo = 1:nDof
        ind = (1:nDof:G.cells.num*nDof) + dofNo - 1;
        dofloc = rldecode(dof(ind), ii, 1);
        s = s + dofloc.*psi{dofNo}(x);
    end
    
end