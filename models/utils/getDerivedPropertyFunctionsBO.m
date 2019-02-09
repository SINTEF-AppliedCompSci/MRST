function [b, mu, rho, mob, f] = getDerivedPropertyFunctionsBO(model, pO, mobMult, status)

    fluid  = model.fluid;
    disgas = model.disgas;
    vapoil = model.vapoil;
    
    if disgas
        isSatO = ~status{1};
    end
    if vapoil
        isSatG = ~status{2};
    end
    
    % Water----------------------------------------------------------------
    if isfield(fluid, 'pcOW')
        pW = @(sW, c) pO(c) - fluid.pcOW(sW);
    else
        pW = @(sW, c) pO(c);
    end
    bW   = @(sW, c) fluid.bW(pW(sW, c));
    muW  = @(sW, c) fluid.muW(pW(sW, c));
    rhoW = @(sW, c) bW(sW, c).*fluid.rhoWS;
    mobW = @(sW, sT, c) mobMult(c).*fluid.krW(sW./sT)./muW(sW, c);
    %----------------------------------------------------------------------
    
    % Oil------------------------------------------------------------------
    if disgas
        bO   = @(rs, c) fluid.bO(pO(c), rs, isSatO(c));
        muO  = @(rs, c) fluid.muO(pO(c), rs, isSatO(c));
        rhoO = @(rs, c) bO(rs, c).*(rs.*fluid.rhoGS + fluid.rhoOS);
    else
        bO   = @(rs, c) fluid.bO(pO(c));
        muO  = @(rs, c) fluid.muO(pO(c));
        rhoO = @(rs, c) fluid.bO(pO(c)).*fluid.rhoOS;
    end
    mobO = @(sO, sT, rs, c) mobMult(c).*fluid.krO(sO./sT)./muO(rs, c);
    %----------------------------------------------------------------------
    
    % Gas------------------------------------------------------------------
    if isfield(fluid, 'pcOG')
        pG = @(sG, c) pO(c) + fluid.pcOG(sG);
    else
        pG = @(sG, c) pO(c);
    end
    if disgas
        bG   = @(sG, rv, c) fluid.bG(pG(sG, c), rv, isSatG(c));
        muG  = @(sG, rv, c) fluid.muG(pG(sG, c), rv, isSatG(c));
        rhoG = @(sG, rv, c) bG(rs, c).*(rs.*fluid.rhoGS + fluid.rhoOS);
    else
        bG   = @(sG, rv, c) fluid.bG(pG(sG, c));
        muG  = @(sG, rv, c) fluid.muG(pG(sG, c));
        rhoG = @(sG, rv, c) fluid.bG(pG(sG, c)).*fluid.rhoGS;
    end
    mobG = @(sG, sT, rv, c) mobMult(c).*fluid.krG(sG./sT)./muG(sG, rv, c);
    %----------------------------------------------------------------------
    
    % Fractional flow functions--------------------------------------------
    mobT = @(sW, sO, sG, sT, cW, cO, cG) mobW(sW, sT, cW) + mobO(sO, sT, cO) + mobG(sG, sT, cG);
    fW   = @(sW, sO, sG, sT, cW, cO, cG) mobW(sW, sT, cW)./mobT(sW, sO, sG, sT, cW, cO, cG);
    fO   = @(sW, sO, sG, sT, cW, cO, cG) mobO(sO, sT, cO)./mobT(sW, sO, sG, sT, cW, cO, cG);
    fG   = @(sW, sO, sG, sT, cW, cO, cG) mobG(sG, sT, cG)./mobT(sW, sO, sG, sT, cW, cO, cG);
    %----------------------------------------------------------------------
    
    % Gather output--------------------------------------------------------
    b   = {bW  , bO  , bG  };
    mu  = {muW , muO , muG };
    rho = {rhoW, rhoO, rhoG};
    mob = {mobW, mobO, mobG};
    f   = {fW  , fO  , fG  };
    %----------------------------------------------------------------------
    
end