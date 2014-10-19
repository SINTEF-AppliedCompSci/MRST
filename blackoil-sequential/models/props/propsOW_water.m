function [bW, rhoW, mobW, dpW, Gw] = propsOW_water(sW, krW, gdz, f, p, s)

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
    
    Gw = rhoWf.*gdz;
    if hasCap
        Gw = Gw - s.Grad(pcOW);
    end
    dpW = s.Grad(p) - Gw;
    
%     bW = ensurePositiveADValue(bW);
%     mobW = ensurePositiveADValue(mobW);

end