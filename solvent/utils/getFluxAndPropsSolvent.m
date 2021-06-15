function [vW, vO, vG, vS, mobW, mobO, mobG, mobS, upcW, upcO, upcG, upcS] = getFluxAndPropsSolvent(fluid, pO, krW, krO, krG, krS, muW, muO, muG, muS, rhoW, rhoO, rhoG, rhoS, sW, sG, sS, T, gdz, op)
% Flux and other properties for the black-oil solvent model equations.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
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
    upcW  = (value(dpW)<=0);
    vW = -op.faceUpstr(upcW, mobW).*T.*dpW;
    
    rhoOf  = op.faceAvg(rhoO);
    mobO   = krO./muO;
    dpO    = op.Grad(pO) - rhoOf.*gdz;
    % oil upstream-index
    upcO = (value(dpO)<=0);
    vO   = - op.faceUpstr(upcO, mobO).*T.*dpO;
    
    pcOG = 0;
    if isfield(fluid, 'pcOG') && ~isempty(sG)
        Mp = fluid.Mp(pO);
        pcOG_m = fluid.pcOG(sG);
        pcOG_i = fluid.pcOG(sG + sS);
        pcOG = Mp.*pcOG_m + (1-Mp).*pcOG_i;
    end
    
    [pG, pS] = deal(pO + pcOG);
    
    rhoGf = op.faceAvg(rhoG);
    mobG  = krG./muG;
    dpG   = op.Grad(pG) - rhoGf.*gdz;
    % gas upstream-index
    upcG  = (value(dpG)<=0);
    vG    = - op.faceUpstr(upcG, mobG).*T.*dpG;
    
    rhoSf  = op.faceAvg(rhoS);
    mobS   = krS./muS;
    dpS    = op.Grad(pS) - rhoSf.*gdz;
    % solvent upstream-index
    upcS    = (value(dpS)<=0);
    vS = - op.faceUpstr(upcS, mobS).*T.*dpS;
    
    
end