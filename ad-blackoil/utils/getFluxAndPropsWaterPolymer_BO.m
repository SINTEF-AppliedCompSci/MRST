function [vW, vP, bW, mobW, mobP, rhoW, pW, upcw, a] = ...
        getFluxAndPropsWaterPolymer_BO(model, pO, sW, c, ads, ...
        krW, T, gdz)
    fluid = model.fluid;
    s = model.operators;
    
    % Check for capillary pressure (p_cow)
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW);
    end
    pW = pO - pcOW;
    
    % Multipliers due to polymer
    mixpar = fluid.mixPar;
    cbar   = c/fluid.cmax;
    a = fluid.muWMult(fluid.cmax).^(1-mixpar);
    b = 1./(1-cbar+cbar./a);
    muWMult = b.*fluid.muWMult(c).^mixpar;
    permRed = 1 + ((fluid.rrf-1)./fluid.adsMax).*ads;
    muWMult  = muWMult.*permRed;
    
    % Water props
    bW     = fluid.bW(pO);
    rhoW   = bW.*fluid.rhoWS;
    % rhoW on face, average of neighboring cells
    rhoWf  = s.faceAvg(rhoW);
    muW    = fluid.muW(pO);
    muWeff = muWMult.*muW;
    mobW   = krW./muWeff;
    dpW    = s.Grad(pO-pcOW) - rhoWf.*gdz;
    % water upstream-index
    upcw = (double(dpW)<=0);
    vW   = -s.faceUpstr(upcw, mobW).*T.*dpW;
    if any(bW < 0)
        warning('Negative water compressibility present!')
    end
    
    % Polymer
    mobP = (mobW.*c)./(a + (1-a)*cbar);
    vP   = - s.faceUpstr(upcw, mobP).*s.T.*dpW;
    
end


