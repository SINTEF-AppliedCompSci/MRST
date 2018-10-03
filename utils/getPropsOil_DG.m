function [vO, bO, mobO, rhoO, pO, upcO, dpO, muO] = getPropsOil_DG(model, pO, T, gdz, mobMult)

    fluid = model.fluid;
    op = model.operators;

    bO     = fluid.bO(pO);
    rhoO   = bO.*fluid.rhoOS;
    rhoOf  = op.faceAvg(rhoO);
    dpO    = op.Grad(pO) - rhoOf.*gdz;
    
    muO = fluid.muO(pO);
    mobO = @(sO,c) mobMult(c).*fluid.krO(sO)./muO(c);
    
    vO = @(sO,c) -mobO(sO,c).*T.*dpO;

    upcO  = (double(dpO)<=0);

    if any(bO < 0)
        warning('Negative water compressibility present!')
    end
     
end