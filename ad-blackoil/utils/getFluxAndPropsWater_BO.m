function [vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, pO, sW, krW, T, gdz)
    fluid = model.fluid;
    s = model.operators;
    % Check for capillary pressure (p_cow)
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW);
    end
    pW = pO - pcOW;
    
    bW     = fluid.bW(pO);
    rhoW   = bW.*fluid.rhoWS;
    % rhoW on face, average of neighboring cells
    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./fluid.muW(pO);
    dpW    = s.Grad(pO-pcOW) - rhoWf.*gdz;
    % water upstream-index
    upcw  = (double(dpW)<=0);
    vW = -s.faceUpstr(upcw, mobW).*T.*dpW;
    if any(bW < 0)
        warning('Negative water compressibility present!')
    end
end
