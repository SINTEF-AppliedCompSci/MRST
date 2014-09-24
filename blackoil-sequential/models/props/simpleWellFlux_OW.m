function [bWqW, bOqO, bWw, bOw, Rw, wc, bWmobWw, bOmobOw, perf2well] = simpleWellFlux_OW(W, b, mob, rho, pressure, pBHP, grav, varargin)
    useSimple = false;
    if numel(varargin)
        useSimple = varargin{1};
    end
    
    rhoW = rho{1}; rhoO = rho{2};
    bW = b{1}; bO = b{2};
    mobW = mob{1}; mobO = mob{2};
    
    
    %WELLS ----------------------------------------------------------------
    [Tw, dzw, Rw, wc, perf2well, pInx, iInxW, iInxO] = getWellStuff(W);
    bWw     = bW(wc);
    bOw     = bO(wc);
    mobWw  = mobW(wc);
    mobOw  = mobO(wc);

    if ~useSimple
%         scale = 1;
%         mobWw(iInxW) = (mobWw(iInxW) + scale*mobOw(iInxW));
%         mobOw(iInxW) = (1 - scale).*mobOw(iInxW);
        mobWw(iInxW) = (mobWw(iInxW) + mobOw(iInxW));
        mobOw(iInxW) = 0;
%         tmp = mobOw(iInxW);
%         if isa(tmp, 'ADI')
%             tmp.val
%         else
%             tmp
%         end
        mobWw(iInxO) = 0; 
        mobOw(iInxO) = (mobWw(iInxO) + mobOw(iInxO));
    end
    bWmobWw  = bWw.*mobWw;
    bOmobOw  = bOw.*mobOw;
    pw  = pressure(wc);

    bWqW  = -bWmobWw.*Tw.*(pBHP(perf2well) - pw + grav*dzw.*rhoW(wc));
    bOqO  = -bOmobOw.*Tw.*(pBHP(perf2well) - pw + grav*dzw.*rhoO(wc));
end