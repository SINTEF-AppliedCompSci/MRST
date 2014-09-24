function [bG, rhoG, mobG, dpG, Gg] = propsBO_gas(sG, krG, grav, dZ, f, p, s)

    if isfield(f, 'tranMultR');
        trMult = f.tranMultR(p);
    else
        trMult = 1;
    end
    
    hasCap = isfield(f, 'pcOG');
    pcOG = 0;
    if hasCap
        pcOG  = f.pcOG(sG);
    end
    
    % oil props
    bG     = f.bG(p);
    rhoG   = bG.*f.rhoGS;
    rhoGf  = s.faceAvg(rhoG);
    mobG = trMult.*krG./f.muG(p);

    
    Gg = grav*(rhoGf.*dZ);
    if hasCap
        Gg = Gg - s.grad(pcOG);
    end
    dpG    = s.grad(p) - Gg;
    
%     bG = ensurePositiveADValue(bG);
%     mobG = ensurePositiveADValue(mobG);

end