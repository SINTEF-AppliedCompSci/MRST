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
    
    [b, mu, rho, mob] = deal(cell(model.water + model.oil + model.gas,1));
    phNo = 1;
    % Water----------------------------------------------------------------
    if model.water
        if isfield(fluid, 'pcOW')
            pW = @(c, sW) pO(c) - fluid.pcOW(sW);
        else
            pW = @(c, varargin) pO(c);
        end
        b{phNo}   = @(c, sW, varargin) fluid.bW(pW(c, sW));
        mu{phNo}  = @(c, sW, varargin) fluid.muW(pW(c, sW));
        rho{phNo} = @(c, sW, varargin) b{phNo}(c, sW).*fluid.rhoWS;
        mob{phNo} = @(c, sW, sT, varargin) mobMult(c).*fluid.krW(sW./sT)./mu{phNo}(c, sW);
        phNo = phNo + 1;
    end
    %----------------------------------------------------------------------
    
    % Oil------------------------------------------------------------------
    if model.oil
        if disgas
            b{phNo}   = @(c, rs, varargin) fluid.bO(pO(c), rs, isSatO(c));
            mu{phNo}  = @(c, rs, varargin) fluid.muO(pO(c), rs, isSatO(c));
            rho{phNo} = @(c, rs, varargin) b{phNo}(c, rs, isSatO(c)).*(rs.*fluid.rhoGS + fluid.rhoOS);
            mob{phNo} = @(c, sO, sT, rs, varargin) mobMult(c).*fluid.krO(sO./sT)./mu{phNo}(c, rv);
        else
            b{phNo}   = @(c, varargin) fluid.bO(pO(c));
            mu{phNo}  = @(c, varargin) fluid.muO(pO(c));
            rho{phNo} = @(c, varargin) fluid.bO(pO(c)).*fluid.rhoOS;
            mob{phNo} = @(c, sO, sT, varargin) mobMult(c).*fluid.krO(sO./sT)./mu{phNo}(c, varargin);
        end
        phNo = phNo + 1;
    end
    %----------------------------------------------------------------------
    
    % Gas------------------------------------------------------------------
    if model.gas
        if isfield(fluid, 'pcOG')
            pG = @(c, sG) pO(c) + fluid.pcOG(sG);
        else
            pG = @(c, varargin) pO(c);
        end
        if vapoil
            b{phNo}   = @(c, sG, rs) fluid.bG(pG(sG, c), rs, isSatG(c));
            mu{phNo}  = @(c, sG, rs) fluid.muG(pG(sG, c), rs, isSatG(c));
            rho{phNo} = @(c, sG, rs) b{phNo}(rs, c).*(rs.*fluid.rhoGS + fluid.rhoOS);
            mob{phNo} = @(c, sG, sT, rs) mobMult(c).*fluid.krG(sG./sT)./mu{phNo}(c, sG, rs);
        else
            b{phNo}   = @(c, sG, varargin) fluid.bG(pG(c, sG));
            mu{phNo}  = @(c, sG, varargin) fluid.muG(pG(c, sG));
            rho{phNo} = @(c, sG, varargin) fluid.bG(pG(c, sG)).*fluid.rhoGS;
            mob{phNo} = @(c, sG, sT, varargin) mobMult(c).*fluid.krG(sG./sT)./mu{phNo}(c, sG, varargin);
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