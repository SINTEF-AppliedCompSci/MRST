function [vG, bG, mobG, rhoG, pG, upcg, dpG, muG] = getFluxAndPropsGas_BO(model, pO, sG, krG, T, gdz, rv, isSat)
%Get flux and properties for the gas phase for a black-oil problem
%
% SYNOPSIS:
%   [vG, bG, mobG, rhoG, pG, upcg, dpG] = ...
%      getFluxAndPropsOil_BO(model, pO, sG, krG, T, gdz, rv, isSat)
% DESCRIPTION:
%   Utility function for evaluating mobilities, densities and fluxes for
%   the gas phase in the black-oil (BO)-model
%
% REQUIRED PARAMETERS:
%   model - ThreePhaseBlackOil or derived subclass instance
%
%   p     - Oil phase pressure. One value per cell in simulation grid.
%
%   sG    - Gas phase saturation. One value per cell in simulation grid.
%
%   krG   - Gas phase relative permeability, possibly evaluated some
%           three-phase relperm function. One value per cell in simulation
%           grid.
%
%   T     - Transmissibility for the internal interfaces.
%
%   gdz   - Let dz be the gradient of the cell depths, on each interface.
%   Then gdz is the dot product of the dz with the gravity vector of the
%   model. See getGravityGradient for more information.
%
%   rv    - Amount of oil vaporized into the gas phase. Use to evaluate the
%           viscosity and density. One value per cell in simulation grid.
%
%   isSat - Boolean indicator if the gas in the cell is saturated with oil.
%   True indicates that the maximum vaporized gas is present inside the gas
%   phase in the cell, and that we need to look at the saturated tables
%   when evaluating rv-dependent properties. One value per cell in
%   simulation grid.
%
% RETURNS:
%   vG    - Gas flux, at reservoir conditions, at each interface.
%
%   bG    - Gas reciprocal formation volume factors. One value per cell in
%           simulation grid.
%
%   rhoG  - Gas density in each cell.
%
%   pG     - Gas phase pressure. One value per cell in simulation grid.
%           Pressure will be different from the oil pressure if capillary
%           pressure is present in the model.
%
%   upcg  - Upwind indicators for the gas phase. One value per interface.
%
%   dpG   - Phase pressure differential for each interface. Was used to
%           evaluate the phase flux.
%
% SEE ALSO:
%   getFluxAndPropsOil_BO, getFluxAndPropsWater_BO

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
    % Check for capillary pressure (p_cow)
    pcOG = 0;
    if isfield(fluid, 'pcOG') && ~isempty(sG)
        pcOG  = fluid.pcOG(sG);
    end
    pG = pO + pcOG;
    
    
    if nargin < 7
        assert(~model.vapoil, 'RS and saturated flag must be supplied for vapoil model');
        rv = 0;
    end
    
    % Gas props (calculated at oil pressure)
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
    rhoGf  = s.faceAvg(rhoG);
    dpG    = s.Grad(pG) - rhoGf.*gdz;
    % gas upstream-index
    upcg    = (value(dpG)<=0);
    
    [krGf, krG] = s.splitFaceCellValue(s, upcg, krG);
    [muGf, muG] = s.splitFaceCellValue(s, upcg, muG);
    mobG   = krG./muG;
    
    vG = - (krGf./muGf).*T.*dpG;
end

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
