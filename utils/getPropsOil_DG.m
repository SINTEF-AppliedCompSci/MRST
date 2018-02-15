function [vO, bO, mobO, rhoO, pO, upcO, dpO, muO] = getPropsOil_DG(model, pO, sO, T, gdz)

    fluid = model.fluid;
    op = model.operators;
    G = model.G;
    
    w = 1;
    psi = {@(x) w*prod(x.^[0,0],2), @(x) w*prod(x.^[1,0],2), @(x) w*prod(x.^[0,1],2)};
%     sO = @(x,c) getSatFromDof(x, c, sOdof, psi, G);

    bO     = fluid.bO(pO);
    rhoO   = bO.*fluid.rhoOS;
    rhoOf  = op.faceAvg(rhoO);
    dpO    = op.Grad(pO) - rhoOf.*gdz;
    
    muO = fluid.muO(pO);
    mobO = @(x,c) fluid.krO(sO(x,c))./muO(c);
    
    vO = @(x,c) -mobO(x,c).*T.*dpO;

    upcO  = (double(dpO)<=0);
    
%     % rhoW on face, average of neighboring cells
%     rhoWf  = op.faceAvg(rhoW);
%     dpW    = op.Grad(pW) - rhoWf.*gdz;
%     % water upstream-index
%     upcw  = (double(dpW)<=0);
%     [krWf, krW] = s.splitFaceCellValue(s, upcw, krW);
%     [muWf, muW] = s.splitFaceCellValue(s, upcw, muW);
%     mobW   = krW./muW;
%     
%     vW = -(krWf./muWf).*T.*dpW;
    if any(bO < 0)
        warning('Negative water compressibility present!')
    end
    
    
end

function s = getSatFromDof(x, cell, dof, psi, G)

    ii = sum((1:G.cells.num)' == cell',2);
    nDof = 3;
    s = 0;
    for dofNo = 1:nDof
        ind = (1:nDof:G.cells.num*nDof) + dofNo - 1;
        dofloc = rldecode(dof(ind), ii, 1);
        s = s + dofloc.*psi{dofNo}(x);
    end
    
end