function [bW, rhoW, mobW, dpW, Gw] = propsOW_water(sW, krW, grav, dZ, f, p, s)

    if isfield(f, 'tranMultR');
        trMult = f.tranMultR(p);
    else
        trMult = 1;
    end
    
    hasCap = isfield(f,  'pcOW');
    
    pcOW = 0;
    if hasCap
        pcOW  = f.pcOW(sW);
    end

    bW     = f.bW(p);
    rhoW   = bW.*f.rhoWS;

    rhoWf  = s.faceAvg(rhoW);
    
    mobW   = trMult.*krW./f.muW(p);
    
    Gw = grav*(rhoWf.*dZ);
    if hasCap
        Gw = Gw - s.grad(pcOW);
    end
    dpW = s.grad(p) - Gw;
    
%     bW = ensurePositiveADValue(bW);
%     mobW = ensurePositiveADValue(mobW);

end