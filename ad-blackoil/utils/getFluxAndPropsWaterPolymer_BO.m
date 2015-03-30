function [vW, vP, bW, muWMult, mobW, mobP, rhoW, pW, upcw, a, dpW] = ...
        getFluxAndPropsWaterPolymer_BO(model, pO, sW, c, ads, ...
        krW, T, gdz)
    f = model.fluid;
    s = model.operators;

    % Check for capillary pressure (p_cow)
    pcOW = 0;
    if isfield(f, 'pcOW') && ~isempty(sW)
        pcOW  = f.pcOW(sW);
    end
    pW = pO - pcOW;

    % Multipliers due to polymer
    mixpar = f.mixPar;
    cbar   = c./f.cmax;
    a      = f.muWMult(f.cmax).^(1-mixpar);
    b      = 1./(1-cbar+cbar./a);
    
    % The viscosity multiplier only result from the polymer mixing.
    muWMult  = f.muWMult(c);
    muWMultT = b.*muWMult.^mixpar;
    permRed  = 1 + ((f.rrf-1)./f.adsMax).*ads;
    muWMultT = muWMultT.*permRed;
    
    % Water props
    bW     = f.bW(pO);
    rhoW   = bW.*f.rhoWS;
    % rhoW on face, average of neighboring cells
    rhoWf  = s.faceAvg(rhoW);
    muW    = f.muW(pO);
    muWeff = muWMultT.*muW;
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
    
    %muPeff = muWeff.*(a + (1-a)*cbar);
    
    
end


