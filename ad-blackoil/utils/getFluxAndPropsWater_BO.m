function [vW, bW, mobW, rhoW, pW, upcw, dpW, muW] = getFluxAndPropsWater_BO(model, pO, sW, krW, T, gdz)
%Get flux and properties for the water phase for a black-oil problem
%
% SYNOPSIS:
%   [vW, bW, mobW, rhoW, pW, upcw, dpW] = ...
%      getFluxAndPropsWater_BO(model, pO, sW, krW, T, gdz)
% DESCRIPTION:
%   Utility function for evaluating mobilities, densities and fluxes for
%   the water phase in the black-oil (BO)-model
%
% REQUIRED PARAMETERS:
%   model - ThreePhaseBlackOil or derived subclass instance
%
%   pO    - Oil phase pressure. One value per cell in simulation grid.
%
%   sW    - Water phase saturation. One value per cell in simulation grid.
%
%   krW   - Water phase relative permeability, possibly evaluated some
%           three-phase relperm function. One value per cell in simulation
%           grid.
%
%   T     - Transmissibility for the internal interfaces.
%
%   gdz   - Let dz be the gradient of the cell depths, on each interface.
%   Then gdz is the dot product of the dz with the gravity vector of the
%   model. See getGravityGradient for more information.
%
% RETURNS:
%   vW    - Water flux, at reservoir conditions, at each interface.
%
%   bW    - Water reciprocal formation volume factors. One value per cell in
%           simulation grid.
%
%   rhoW  - Water density in each cell.
%
%   pW    - Water phase pressure. One value per cell in simulation grid.
%           Pressure will be different from the oil pressure if capillary
%           pressure is present in the model.
%
%   upcw  - Upwind indicators for the water phase. One value per interface.
%
%   dpW   - Phase pressure differential for each interface. Was used to
%           evaluate the phase flux.
%
% SEE ALSO:
%   getFluxAndPropsGas_BO, getFluxAndPropsOil_BO

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
    fluid = model.fluid;
    s = model.operators;
    fp = model.FlowPropertyFunctions;
    pp = model.PVTPropertyFunctions;
    if isempty(fp)
        fp_sat = StateFunction(model);
        fp_pvt = fp_sat;
    else
        fp_pvt = pp.Density;
        fp_sat = fp.RelativePermeability;
    end
    pvtfun = @(varargin) fp_pvt.evaluateFunctionOnDomainWithArguments(varargin{:});
    satfun = @(varargin) fp_sat.evaluateFunctionOnDomainWithArguments(varargin{:});
    % Check for capillary pressure (p_cow)
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = satfun(fluid.pcOW, sW);
    end
    pW = pO - pcOW;
    muW = pvtfun(fluid.muW, pW);
    
    bW     = pvtfun(fluid.bW, pW);
    rhoW   = bW.*fluid.rhoWS;
    % rhoW on face, average of neighboring cells
    rhoWf  = s.faceAvg(rhoW);
    dpW    = s.Grad(pW) - rhoWf.*gdz;
    % water upstream-index
    upcw  = (value(dpW)<=0);
    [krWf, krW] = s.splitFaceCellValue(s, upcw, krW);
    [muWf, muW] = s.splitFaceCellValue(s, upcw, muW);
    mobW   = krW./muW;
    
    vW = -(krWf./muWf).*T.*dpW;
    if any(bW < 0)
        warning('Negative water compressibility present!')
    end
end


