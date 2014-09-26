function [bO, rhoO, mobO, dpO, Go] = propsOW_oil(sO, krO, grav, dZ, f, p, s)

    if isfield(f, 'tranMultR');
        trMult = f.tranMultR(p);
    else
        trMult = 1;
    end

    % oil props
    bO     = f.bO(p);
    rhoO   = bO.*f.rhoOS;
    rhoOf  = s.faceAvg(rhoO);
    
    Go = grav*(rhoOf.*dZ);
    dpO    = s.grad(p) - Go;
    
    if isfield(f, 'BOxmuO')
        % mob0 is already multplied with b0
        mobO   = (trMult.*krO./f.BOxmuO(p))./bO;
    else
        mobO   = trMult.*krO./f.muO(p);
    end

end