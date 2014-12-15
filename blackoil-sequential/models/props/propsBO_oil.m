function [bO, rhoO, mobO, dpO, Go] = propsBO_oil(state, krO, rs, isSat, grav, dZ, f, p, s)

    if isfield(f, 'tranMultR');
        trMult = f.tranMultR(p);
    else
        trMult = 1;
    end

    % oil props
    bO     = f.bO(p, rs, isSat);
    
    rhoO   = bO.*(rs*f.rhoGS + f.rhoOS);
    rhoOf  = s.faceAvg(rhoO);
    
    Go = grav*(rhoOf.*dZ);
    dpO    = s.grad(p) - Go;
    
    muO = f.muO(p,rs,isSat);
%     muO = double(muO);
%     bad = double(muO) < 0;
%     assert(~any(bad))
%     if any(bad)
%         muO(bad) = f.muO(p(bad), rs(bad), ~isSat(bad));
%     end
    
    mobO   = trMult.*krO./muO;

%     bad = double(muO)<0;
%     if any(bad)
%         mobO(bad) = mobO(bad)./mobO(bad) - 1;
%     end
    %negative viscosity can happen - what to do
%     bO = ensurePositiveADValue(bO);
%     mobO = ensurePositiveADValue(mobO);
end