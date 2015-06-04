function [vG, bG, mobG, rhoG, pG, upcg, dpG] = getFluxAndPropsGas_BO(model, pO, p_prop, sG, krG, T, gdz, rv, isSat)
    fluid = model.fluid;
    s = model.operators;
    % Check for capillary pressure (p_cow)
    pcOG = 0;
    if isfield(fluid, 'pcOG') && ~isempty(sG)
        pcOG  = fluid.pcOG(sG);
    end
    pG = p_prop + pcOG;
    
    
    if nargin < 7
        assert(~model.vapoil, 'RS and saturated flag must be supplied for vapoil model');
        rv = 0;
    end
    
    % Gas props (calculated at oil pressure)
    if model.vapoil
        bG  = fluid.bG(p_prop, rv, isSat);
        muG = fluid.muG(p_prop, rv, isSat);
    else
        bG  = fluid.bG(p_prop);
        muG = fluid.muG(p_prop);
    end
    
    if any(bG < 0)
        warning('Negative gas compressibility present!')
    end
    rhoG   = bG.*(rv*fluid.rhoOS + fluid.rhoGS);
    rhoGf  = s.faceAvg(rhoG);
    mobG   = krG./muG;
    dpG    = s.Grad(pO+pcOG) - rhoGf.*gdz;
    % gas upstream-index
    upcg    = (double(dpG)<=0);
    vG = - s.faceUpstr(upcg, mobG).*T.*dpG;
end
