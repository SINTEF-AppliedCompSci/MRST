function [b, mu, rho, mob] = getDerivedPropertyFunctionsBO(model, pO, mobMult, status)

    fluid  = model.fluid;
    disgas = model.disgas;
    vapoil = model.vapoil;
    
    if disgas
        isSatO = ~status{1};
    end
    if vapoil
        isSatG = ~status{2};
    end
    
    [b, mu, rho, mob] = deal(cell(model.water + model.oil + model.gas,1));
    eqNo = 1;
    % Water----------------------------------------------------------------
    if model.water
        if isfield(fluid, 'pcOW')
            pW = @(c, sW) pO(c) - fluid.pcOW(sW);
        else
            pW = @(c, varargin) pO(c);
        end
        b{eqNo}   = @(c, sW, varargin) fluid.bW(pW(c, sW));
        mu{eqNo}  = @(c, sW, varargin) fluid.muW(pW(c, sW));
        rho{eqNo} = @(c, sW, varargin) b{eqNo}(c, sW).*fluid.rhoWS;
        mob{eqNo} = @(c, sW, sT, varargin) mobMult(c).*fluid.krW(sW./sT)./mu{eqNo}(c, sW);
        eqNo = eqNo + 1;
    end
    %----------------------------------------------------------------------
    
    % Oil------------------------------------------------------------------
    if model.oil
        if disgas
            b{eqNo}   = @(c, sO, rS, varargin) fluid.bO(pO(c), rS, isSatO(c));
            mu{eqNo}  = @(c, sO, rS, varargin) fluid.muO(pO(c), rS, isSatO(c));
            rho{eqNo} = @(c, sO, rS, varargin) bO(c, rS).*(rS.*fluid.rhoGS + fluid.rhoOS);
            mob{eqNo} = @(c, sO, sT, rS, varargin) mobMult(c).*fluid.krO(sO./sT)./mu{eqNo}(c, rS);
        else
            b{eqNo}   = @(c, varargin) fluid.bO(pO(c));
            mu{eqNo}  = @(c, varargin) fluid.muO(pO(c));
            rho{eqNo} = @(c, varargin) fluid.bO(pO(c)).*fluid.rhoOS;
            mob{eqNo} = @(c, sO, sT, varargin) mobMult(c).*fluid.krO(sO./sT)./mu{eqNo}(c, varargin);
            eqNo = eqNo + 1;
        end
    end
    %----------------------------------------------------------------------
    
    % Gas------------------------------------------------------------------
    if model.gas
        if isfield(fluid, 'pcOG')
            pG = @(c, sG) pO(c) + fluid.pcOG(sG);
        else
            pG = @(c, varargin) pO(c);
        end
        if disgas
            b{eqNo}   = @(c, sG, rV) fluid.bG(pG(sG, c), rV, isSatG(c));
            mu{eqNo}  = @(c, sG, rV) fluid.muG(pG(sG, c), rV, isSatG(c));
            rho{eqNo} = @(c, sG, rV) b{eqNo}(rs, c).*(rv.*fluid.rhoOS + fluid.rhoGS);
            mob{eqNo} = @(c, sG, sT, rV) mobMult(c).*fluid.krG(sG./sT)./mu{eqNo}(c, sG, rV);
        else
            b{eqNo}   = @(c, sG, varargin) fluid.bG(pG(c, sG, varargin));
            mu{eqNo}  = @(c, sG, varargin) fluid.muG(pG(c, sG, varargin));
            rho{eqNo} = @(c, sG, varargin) fluid.bG(pG(c, sG, varargin)).*fluid.rhoGS;
            mob{eqNo} = @(c, sG, sT, varargin) mobMult(c).*fluid.krG(sG./sT)./mu{eqNo}(c, sG, varargin);
        end
    end
    %----------------------------------------------------------------------
    
%     % Fractional flow functions--------------------------------------------
%     mobT = @(sW, sO, sG, sT, cW, cO, cG) mobW(sW, sT, cW) + mobO(sO, sT, cO) + mobG(sG, sT, cG);
%     fW   = @(sW, sO, sG, sT, cW, cO, cG) mobW(sW, sT, cW)./mobT(sW, sO, sG, sT, cW, cO, cG);
%     fO   = @(sW, sO, sG, sT, cW, cO, cG) mobO(sO, sT, cO)./mobT(sW, sO, sG, sT, cW, cO, cG);
%     fG   = @(sW, sO, sG, sT, cW, cO, cG) mobG(sG, sT, cG)./mobT(sW, sO, sG, sT, cW, cO, cG);
%     %----------------------------------------------------------------------
%     
%     % Gather output--------------------------------------------------------
%     b   = {bW  , bO  , bG  };
%     mu  = {muW , muO , muG };
%     rho = {rhoW, rhoO, rhoG};
%     mob = {mobW, mobO, mobG};
%     f   = {fW  , fO  , fG  };
%     %----------------------------------------------------------------------
    
end