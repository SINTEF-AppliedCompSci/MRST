function [vO, bO, mobO, rhoO, p, upco, dpO] = getFluxAndPropsOil_BO(model, p, p_prop, sO, krO, T, gdz, rs, isSat)
    disgas = isprop(model, 'disgas') && model.disgas;
    
    if nargin < 7
        assert(~disgas, 'RS and saturated flag must be supplied for disgas model');
        rs = 0;
    end
    
    fluid = model.fluid;
    s = model.operators;
    % Oil props
    if disgas
        bO  = fluid.bO(p_prop,  rs, isSat);
        muO = fluid.muO(p_prop, rs, isSat);
        rhoO   = bO.*(rs*fluid.rhoGS + fluid.rhoOS);
    else
        bO  = fluid.bO(p_prop);
        if isfield(fluid, 'BOxmuO')
            muO = fluid.BOxmuO(p_prop).*bO;
        else
            muO = fluid.muO(p_prop);
        end
        rhoO   = bO.*fluid.rhoOS;
    end
        
    if any(bO < 0)
        warning('Negative oil compressibility present!')
    end
    
    rhoOf  = s.faceAvg(rhoO);
    mobO   = krO./muO;
    dpO    = s.Grad(p) - rhoOf.*gdz;
    % oil upstream-index
    upco = (double(dpO)<=0);
    vO   = - s.faceUpstr(upco, mobO).*T.*dpO;
end
