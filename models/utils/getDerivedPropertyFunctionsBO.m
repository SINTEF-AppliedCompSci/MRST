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
    
    phNo = 1;
    % Water----------------------------------------------------------------
    if model.water
        ix = model.getPhaseIndex('W');
        krW = @(s) model.evaluateRelPerm('W', s);
        if isfield(fluid, 'pcOW')
            pW = @(c, sW) pO(c) - fluid.pcOW(sW);
        else
            pW = @(c, varargin) pO(c);
        end
        b{phNo}   = @(c, sW, varargin) fluid.bW(pW(c, sW));
        mu{phNo}  = @(c, sW, varargin) fluid.muW(pW(c, sW));
        rho{phNo} = @(c, sW, varargin) b{phNo}(c, sW).*fluid.rhoWS;
        mob{phNo} = @(c, s , varargin) mobMult(c).*krW(s)./mu{phNo}(c, s{ix});
        phNo = phNo + 1;
    end
    %----------------------------------------------------------------------
    
    % Oil------------------------------------------------------------------
    if model.oil
        ix = model.getPhaseIndex('O');
        krO = @(s) model.evaluateRelPerm('O', s);
        if disgas
            b{phNo}   = @(c, sO, rS, varargin) fluid.bO(pO(c), rS, isSatO(c));
            mu{phNo}  = @(c, sO, rS, varargin) fluid.muO(pO(c), rS, isSatO(c));
            rho{phNo} = @(c, sO, rS, varargin) b{phNo}(c, sO, rS).*(rS.*fluid.rhoGS + fluid.rhoOS);
            mob{phNo} = @(c, s , rS, varargin) mobMult(c).*krO(s)./mu{phNo}(c, s{ix}, rS);
        else
            b{phNo}   = @(c, varargin) fluid.bO(pO(c));
            mu{phNo}  = @(c, varargin) fluid.muO(pO(c));
            rho{phNo} = @(c, varargin) fluid.bO(pO(c)).*fluid.rhoOS;
            mob{phNo} = @(c, s, varargin) mobMult(c).*krO(s)./mu{phNo}(c);
        end
        phNo = phNo + 1;
    end
    %----------------------------------------------------------------------
    
    % Gas------------------------------------------------------------------
    if model.gas
        ix = model.getPhaseIndex('G');
        krG = @(s) model.evaluateRelPerm('G', s);
        if isfield(fluid, 'pcOG')
            pG = @(c, sG) pO(c) + fluid.pcOG(sG);
        else
            pG = @(c, varargin) pO(c);
        end
        if vapoil
            b{phNo}   = @(c, sG, rV) fluid.bG(pG(c, sG), rV, isSatG(c));
            mu{phNo}  = @(c, sG, rV) fluid.muG(pG(c, sG), rV, isSatG(c));
            rho{phNo} = @(c, sG, rV) b{phNo}(c, sG, rV).*(rV.*fluid.rhoOS + fluid.rhoGS);
            mob{phNo} = @(c, s , rV) mobMult(c).*krG(s)./mu{phNo}(c, s{ix}, rV);
        else
            b{phNo}   = @(c, sG, varargin) fluid.bG(pG(c, sG));
            mu{phNo}  = @(c, sG, varargin) fluid.muG(pG(c, sG));
            rho{phNo} = @(c, sG, varargin) fluid.bG(pG(c, sG)).*fluid.rhoGS;
            mob{phNo} = @(c, s , varargin) mobMult(c).*krG(s)./mu{phNo}(c, s{ix});
        end
    end
    %----------------------------------------------------------------------
    
end