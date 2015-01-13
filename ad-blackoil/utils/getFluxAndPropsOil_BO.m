function [vO, bO, mobO, rhoO, p, upco] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, isSat)

    if nargin < 7
        assert(~model.disgas, 'RS and saturated flag must be supplied for disgas model');
        rs = 0;
    end
    
    fluid = model.fluid;
    s = model.operators;
    % Oil props
    if model.disgas
        bO  = fluid.bO(p,  rs, isSat);
        muO = fluid.muO(p, rs, isSat);
    else
        bO  = fluid.bO(p);
        muO = fluid.muO(p);
    end
        
    if any(bO < 0)
        warning('Negative oil compressibility present!')
    end
    rhoO   = bO.*(rs*fluid.rhoGS + fluid.rhoOS);
    rhoOf  = s.faceAvg(rhoO);
    mobO   = krO./muO;
    dpO    = s.Grad(p) - rhoOf.*gdz;
    % oil upstream-index
    upco = (double(dpO)<=0);
    vO   = - s.faceUpstr(upco, mobO).*T.*dpO;
end