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
        pW = @(c, sW) pO(c) - fluid.pcOW(sW);
    else
        pW = @(c, varargin) pO(c);
    end
    bW   = @(c, sW, varargin) fluid.bW(pW(c, sW));
    muW  = @(c, sW, varargin) fluid.muW(pW(c, sW));
    rhoW = @(c, sW, varargin) bW(c, sW).*fluid.rhoWS;
    mobW = @(c, sW, sT, varargin) mobMult(c).*fluid.krW(sW./sT)./muW(c, sW);
    %----------------------------------------------------------------------
    
    % Oil------------------------------------------------------------------
    if disgas
        bO   = @(c, rv, varargin) fluid.bO(pO(c), rs, isSatO(c));
        muO  = @(c, rv, varargin) fluid.muO(pO(c), rs, isSatO(c));
        rhoO = @(c, rv, varargin) bO(c, rs).*(rs.*fluid.rhoGS + fluid.rhoOS);
        mobO = @(c, sO, sT, rv, varargin) mobMult(c).*fluid.krO(sO./sT)./muO(c, rv);
    else
        bO   = @(c, varargin) fluid.bO(pO(c));
        muO  = @(c, varargin) fluid.muO(pO(c));
        rhoO = @(c, varargin) fluid.bO(pO(c)).*fluid.rhoOS;
        mobO = @(c, sO, sT, varargin) mobMult(c).*fluid.krO(sO./sT)./muO(c, varargin);
    end
    %----------------------------------------------------------------------
    
    % Gas------------------------------------------------------------------
    if isfield(fluid, 'pcOG')
        pG = @(c, sG) pO(c) + fluid.pcOG(sG);
    else
        pG = @(c, varargin) pO(c);
    end
    if disgas
        bG   = @(c, sG, rs) fluid.bG(pG(sG, c), rs, isSatG(c));
        muG  = @(c, sG, rs) fluid.muG(pG(sG, c), rs, isSatG(c));
        rhoG = @(c, sG, rs) bG(rs, c).*(rs.*fluid.rhoGS + fluid.rhoOS);
        mobG = @(c, sG, sT, rs) mobMult(c).*fluid.krG(sG./sT)./muG(c, sG, rs);
    else
        bG   = @(c, sG, varargin) fluid.bG(pG(c, sG, varargin));
        muG  = @(c, sG, varargin) fluid.muG(pG(c, sG, varargin));
        rhoG = @(c, sG, varargin) fluid.bG(pG(c, sG, varargin)).*fluid.rhoGS;
        mobG = @(c, sG, sT, varargin) mobMult(c).*fluid.krG(sG./sT)./muG(c, sG, varargin);
    end
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