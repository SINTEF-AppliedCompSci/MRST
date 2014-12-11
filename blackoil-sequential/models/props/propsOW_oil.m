function [bO, rhoO, mobO, Go] = propsOW_oil(sO, krO, gdz, f, p, s)

    if isfield(f, 'tranMultR');
        trMult = f.tranMultR(p);
    else
        trMult = 1;
    end

    % oil props
    bO     = f.bO(p);
    rhoO   = bO.*f.rhoOS;
    rhoOf  = s.faceAvg(rhoO);
    
    Go = rhoOf.*gdz;
    
    if isfield(f, 'BOxmuO')
        % mob0 is already multplied with b0
        mobO   = (trMult.*krO./f.BOxmuO(p))./bO;
    else
        mobO   = trMult.*krO./f.muO(p);
    end

end
