function [bW, rhoW, mobW, Gw, muW] = propsOW_water(sW, krW, gdz, f, p, s)

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
    
    muW = f.muW(p);
    mobW   = trMult.*krW./muW;
    
    Gw = rhoWf.*gdz;
    if hasCap
        Gw = Gw - s.Grad(pcOW);
    end

end
