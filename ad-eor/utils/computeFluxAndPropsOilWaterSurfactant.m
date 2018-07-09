function [dp, mob, upc, b, rho, pvMult, b0, pvMult0, T] =  computeFluxAndPropsOilWaterSurfactant(model, p0, p, sW, c, ...
                                                      pBH, W, varargin)
%
%
% SYNOPSIS:
%   function [dpO, dpW, mobO, mobW, upco, upcw, bO, bW, pvMult, bO0, bW0, pvMult0, T] =  computeFluxAndPropsOilWaterSurfactant(model, p0, p, sW, c, pBH, W, varargin)
%
% DESCRIPTION: Given the state variable (pressure, saturation and
% concentration), compute fluxes and other properties, as listed below.
%
% PARAMETERS:
%   model    - Model instance
%   p0       - Pressure at previous time step
%   p        - Pressure at current time step
%   sW       - Saturation
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
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
    [krW, krO] = computeRelPermSft(sW, c, Nc, fluid);

    % Multipliers for properties
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

    % Compute transmissibility
    T = op.T.*transMult;

    % Gravity contribution
    gdz = model.getGravityGradient();

    % Modifiy relperm by mobility multiplier (if any)
    krW = mobMult.*krW;
    krO = mobMult.*krO;

    fluid = model.fluid;

    % Capillary pressure
    pcOW = 0;
    if isfield(fluid, 'pcOW')
        pcOW  = fluid.pcOW(sW);
    end
    pcOW = pcOW.*fluid.ift(c)/fluid.ift(0);
    pO = p;
    pW = pO - pcOW;

    bW0 = fluid.bW(p0);
    bO0 = fluid.bO(p0);

    % Water flux and properties
    bW      = fluid.bW(p);
    rhoW    = bW.*fluid.rhoWS;
    rhoWf   = op.faceAvg(rhoW);
    muW     = fluid.muWSft(c);
    multmuW = fluid.muW(p)/fluid.muWr;
    mobW    = krW./(muW.*multmuW);
    dpW     = op.Grad(pW) - rhoWf.*gdz;
    upcW    = (double(dpW)<=0);

    % Oil flux and properties
    bO    = fluid.bO(pO);
    rhoO  = bO.*fluid.rhoOS;
    rhoOf = op.faceAvg(rhoO);
    if isfield(fluid, 'BOxmuO')
        muO = fluid.BOxmuO(pO).*bO;
    else
        muO = fluid.muO(pO);
    end
    mobO = krO./muO;
    dpO  = op.Grad(pO) - rhoOf.*gdz;
    upcO = (double(dpO)<=0);

    dp  = {dpW, dpO};
    mob = {mobW, mobO};
    rho = {rhoW, rhoO};
    upc = {upcW, upcO};
    b   = {bW, bO};
    b0  = {bW0, bO0};

end
