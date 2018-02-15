function [vW, bW, mobW, rhoW, pW, upcW, dpW, muW] = getPropsWater_DG(model, pO, sW, T, gdz)


    fluid = model.fluid;
    op = model.operators;
    G = model.G;
    
     % Check for capillary pressure (p_cow)
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW);
    end
    pW = pO - pcOW;
    muW = fluid.muW(pW);
    
   

    bW     = fluid.bW(pW);
    rhoW   = bW.*fluid.rhoWS;
    rhoWf  = op.faceAvg(rhoW);
    dpW    = op.Grad(pW) - rhoWf.*gdz;
    
    muW = fluid.muW(pW);
    mobW = @(x,c) fluid.krW(sW(x,c))./muW(c);
    
    vW = @(x,c) -mobW(x,c).*T.*dpW;

    upcW  = (double(dpW)<=0);
    
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
    if any(bW < 0)
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