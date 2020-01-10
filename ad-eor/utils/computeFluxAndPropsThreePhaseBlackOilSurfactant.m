function [p1, dp, mob, upc, b, rho, pvMult, b0, pvMult0, T] =  ...
    computeFluxAndPropsThreePhaseBlackOilSurfactant(model, p0, p, sW, sG, c, ...
                                                      pBH, W, rs, rs0, rv, rv0, st, st0, varargin)
%
%
% SYNOPSIS:
%   function [p, dpO, dpW, mobO, mobW, upco, upcw, bO, bW, pvMult, bO0, bW0, pvMult0, T] =  
%   computeFluxAndPropsThreePhaseBlackOilSurfactant(model, p0, p, sW, c, pBH, W, rs, rv, st, st0, varargin)
%
% DESCRIPTION: Given the state variable (pressure, saturation and
% concentration), compute fluxes and other properties, as listed below.
%
% PARAMETERS:
%   model    - Model instance
%   p0       - Pressure at previous time step
%   p        - Pressure at current time step
%   sW       - Water Saturation
%   sG       - Gas Saturation
%   c        - Surfactant concentration
%   pBH      - bottom hole pressure
%   W        - Well structure
%   varargin -
%
% RETURNS:
%   dpO     - Oil pressure gradient (face-valued)
%   dpW     - Water pressure gradient (face-valued)
%   mobO    - Oil mobility
%   mobW    - Water mobility
%   upco    - Upstream direction for oil (face-values)
%   upcw    - Upstream direction for water (face-values)
%   bO      - Oil formation volume factor
%   bW      - water formation volume factor
%   pvMult  - Pore volume multiplier
%   bO0     - Oil formation volume factor for previous time step
%   bW0     - Water formation volume factor for previous time step
%   pvMult0 - Pore volume multiplier for previous time step
%   T       - Transmissibilities
%
% EXAMPLE:
%
% SEE ALSO:
%

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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


    opt = struct('velocCompMethod', 'square');
    opt = merge_options(opt, varargin{:});

    G     = model.G;
    fluid = model.fluid;
    op    = model.operators;

    Nc = computeCapillaryNumber(p, c, pBH, W, fluid, G, op, 'velocCompMethod', opt.velocCompMethod);
    [krW, krO, krG] = computeThreePhaseRelPermSft(sW, sG, c, Nc, fluid);

    if any(sW > 0.95)
        1;
    end
    
    
    % Multipliers for properties
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

    % Compute transmissibility
    T = op.T.*transMult;

    % Gravity contribution
    gdz = model.getGravityGradient();

    % Modifiy relperm by mobility multiplier (if any)
    krW = mobMult.*krW;
    krO = mobMult.*krO;
    krG = mobMult.*krG;

    fluid = model.fluid;

    % Capillary pressure
    pcOW = 0;
    pcOG = 0;
    if isfield(fluid, 'pcOW') 
        pcOW  = fluid.pcOW(sW);
    end
    if isfield(fluid, 'pcOG') && ~isempty(sG)
        pcOG  = fluid.pcOG(sG);
    end
    pcOW = pcOW.*fluid.ift(c)/fluid.ift(0);
    pO = p;
    pW = pO - pcOW;
    pG = pO + pcOG;

    bW0 = fluid.bW(p0);
%     bO0 = fluid.bO(p0);
    if isprop(model, 'disgas') && model.disgas
        bO0 = fluid.bO(p0, rs0, ~st0{1});
    else
        bO0 = fluid.bO(p0);
    end
    
%     bG0 = fluid.bG(p0);
    if model.vapoil
        bG0 = fluid.bG(p0, rv0, ~st0{2});
    else
        bG0 = fluid.bG(p0);
    end

    % Water flux and properties
    bW      = fluid.bW(pW);
    rhoW    = bW.*fluid.rhoWS;
    rhoWf   = op.faceAvg(rhoW);
    muW     = fluid.muWSft(c);
    multmuW = fluid.muW(p)/fluid.muWr;
    mobW    = krW./(muW.*multmuW);
    dpW     = op.Grad(pW) - rhoWf.*gdz;
    upcW    = (double(dpW)<=0);

    % Oil flux and properties
    isSat = ~st{1};
    disgas = isprop(model, 'disgas') && model.disgas;
    
    if nargin < 7
        assert(~disgas, 'RS and saturated flag must be supplied for disgas model');
        rs = 0;
    end
    if disgas
        bO  = fluid.bO(pO,  rs, isSat);
        muO = fluid.muO(pO, rs, isSat);
        rhoO   = bO.*(rs*fluid.rhoGS + fluid.rhoOS);
    else
        bO  = fluid.bO(pO);
        if isfield(fluid, 'BOxmuO')
            muO = fluid.BOxmuO(pO).*bO;
        else
            muO = fluid.muO(pO);
        end
        rhoO   = bO.*fluid.rhoOS;
    end
    
    if any(bO < 0)
        warning('Negative oil compressibility present!')
    end
    
    rhoOf = op.faceAvg(rhoO);
    mobO = krO./muO;
    dpO  = op.Grad(pO) - rhoOf.*gdz;
    upcO = (double(dpO)<=0);
   
    % Gas flux and properties
    isSat = ~st{2};
    if model.vapoil
        bG  = fluid.bG(pG, rv, isSat);
        muG = fluid.muG(pG, rv, isSat);
    else
        bG  = fluid.bG(pG);
        muG = fluid.muG(pG);
    end
    if any(bG < 0)
        warning('Negative gas compressibility present!')
    end
    rhoG   = bG.*(rv*fluid.rhoOS + fluid.rhoGS);
    rhoGf  = op.faceAvg(rhoG);
    mobG = krG./muG;
    dpG  = op.Grad(pG) - rhoGf.*gdz;
    upcG = (double(dpG)<=0);
    
    p1  = {pW, pO, pG};
    dp  = {dpW, dpO, dpG};
    mob = {mobW, mobO, mobG};
    rho = {rhoW, rhoO, rhoG};
    upc = {upcW, upcO, upcG};
    b   = {bW, bO, bG};
    b0  = {bW0, bO0, bG0};

end
