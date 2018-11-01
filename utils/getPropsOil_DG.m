function [vO, bO, mobO, rhoO, pO, upcO, dpO, muO] = getPropsOil_DG(model, pO, T, gdz, mobMult)

    fluid = model.fluid;
    op = model.operators;

    bO     = fluid.bO(pO);
    rhoO   = bO.*fluid.rhoOS;
    rhoOf  = op.faceAvg(rhoO);
    dpO    = op.Grad(pO) - rhoOf.*gdz;
    
    muO = fluid.muO(pO);
    mobO = @(sO,sT,c) mobMult(c).*fluid.krO(sO./sT)./muO(c);
    
    vO = @(sO,sT,c) -mobO(sO,sT,c).*T.*dpO;

    upcO  = (double(dpO)<=0);

    if any(bO < 0)
        warning('Negative water compressibility present!')
    end
     
end