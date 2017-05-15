function [vW, vO, vG, vS, mobW, mobO, mobG, mobS, upcW, upcO, upcG, upcS] = getFluxAndPropsSolvent(fluid, pO, krW, krO, krG, krS, muW, muO, muG, muS, rhoW, rhoO, rhoG, rhoS, T, gdz, op)
    
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
       pcOW = fluid.pcOW(sW); 
    end
    pW = pO - pcOW;
    
    % rhoW on face, average of neighboring cells
    rhoWf  = op.faceAvg(rhoW);
    mobW   = krW./muW;
    dpW    = op.Grad(pW) - rhoWf.*gdz;
    % water upstream-index
    upcW  = (double(dpW)<=0);
    vW = -op.faceUpstr(upcW, mobW).*T.*dpW;
    
    rhoOf  = op.faceAvg(rhoO);
    mobO   = krO./muO;
    dpO    = op.Grad(pO) - rhoOf.*gdz;
    % oil upstream-index
    upcO = (double(dpO)<=0);
    vO   = - op.faceUpstr(upcO, mobO).*T.*dpO;
    
    pcOG = 0;
    if isfield(fluid, 'pcOG') && ~isempty(sG)
        M = fluid.Mpres(p);
        pcOG_i = fluid.pcOW(sG);
        pcOG_m = fluid.pcOW_m(sG + sS);
        pcOG = M.*pcOG_i + M.*pcOG_m;
    end
    
    [pG, pS] = deal(pO + pcOG);
    
    rhoGf  = op.faceAvg(rhoG);
    mobG   = krG./muG;
    dpG    = op.Grad(pG) - rhoGf.*gdz;
    % gas upstream-index
    upcG    = (double(dpG)<=0);
    vG = - op.faceUpstr(upcG, mobG).*T.*dpG;
    
    rhoSf  = op.faceAvg(rhoS);
    mobS   = krS./muS;
    dpS    = op.Grad(pS) - rhoSf.*gdz;
    % solvent upstream-index
    upcS    = (double(dpS)<=0);
    vS = - op.faceUpstr(upcS, mobS).*T.*dpS;
    
    
end