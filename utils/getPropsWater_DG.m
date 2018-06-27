function [vW, bW, mobW, rhoW, pW, upcW, dpW, muW] = getPropsWater_DG(model, pO, T, gdz, mobMult)


    fluid = model.fluid;
    op = model.operators;
    
    % Check for capillary pressure (p_cow) (currently not supported)
    assert(~isfield(fluid, 'pcOW'));
    pcOW  = 0;
    if 0
    pcOW = @(sW) 0.*sW;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW);
    end
    end
    pW = pO - pcOW;

    bW    = fluid.bW(pW);
    rhoW  = bW.*fluid.rhoWS;
    rhoWf = op.faceAvg(rhoW);
    dpW   = op.Grad(pW) - rhoWf.*gdz;
    
    muW = fluid.muW(pW);
    mobW = @(sW,c) mobMult(c).*fluid.krW(sW)./muW(c);
    
    vW = @(x,c) -mobW(x,c).*T.*dpW;

    upcW  = (double(dpW)<=0);
    
    if any(bW < 0)
        warning('Negative water compressibility present!')
    end
    
end