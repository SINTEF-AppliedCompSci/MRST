function [vO, bO, mobO, rhoO, p, upco, dpO, muO] = getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, isSat)
%Get flux and properties for the oil phase for a black-oil problem
%
% SYNOPSIS:
%   [vO, bO, mobO, rhoO, p, upco, dpO] = ...
%      getFluxAndPropsOil_BO(model, p, sO, krO, T, gdz, rs, isSat)
% DESCRIPTION:
%   Utility function for evaluating mobilities, densities and fluxes for
%   the oil phase in the black-oil (BO)-model
%
% REQUIRED PARAMETERS:
%   model - ThreePhaseBlackOil or derived subclass instance
%
%   p     - Oil phase pressure. One value per cell in simulation grid.
%
%   sO    - Oil phase saturation. One value per cell in simulation grid.
%
%   krO   - Oil phase relative permeability, possibly evaluated some
%           three-phase relperm function. One value per cell in simulation
%           grid.
%
%   T     - Transmissibility for the internal interfaces.
%
%   gdz   - Let dz be the gradient of the cell depths, on each interface.
%   Then gdz is the dot product of the dz with the gravity vector of the
%   model. See getGravityGradient for more information.
%
%   rs    - Amount of gas dissolved into the oil phase. Use to evaluate the
%           viscosity and density. One value per cell in simulation grid.
%
%   isSat - Boolean indicator if the oil in the cell is saturated with gas.
%   True indicates that there is free gas in the cell, and that we need to
%   look at the saturated tables when evaluating rs-dependent properties.
%   One value per cell in simulation grid.
%
% RETURNS:
%   vO    - Oil flux, at reservoir conditions, at each interface.
%
%   bO    - Oil reciprocal formation volume factors. One value per cell in
%           simulation grid.
%
%   rhoO  - Oil density in each cell.
%
%   p     - Oil phase pressure. One value per cell in simulation grid.
%           Should be unmodified. Included simply to have parity with
%           analogous gas and water functions.
%
%   upco  - Upwind indicators for the oil phase. One value per interface.
%
%   dpO   - Phase pressure differential for each interface. Was used to
%           evaluate the phase flux.
%
% SEE ALSO:
%   getFluxAndPropsGas_BO, getFluxAndPropsWater_BO

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
    disgas = isprop(model, 'disgas') && model.disgas;
    
    if nargin < 7
        assert(~disgas, 'RS and saturated flag must be supplied for disgas model');
        rs = 0;
    end
    
    fluid = model.fluid;
    s = model.operators;
    % Oil props
    if disgas
        bO  = fluid.bO(p,  rs, isSat);
        muO = fluid.muO(p, rs, isSat);
        rhoO   = bO.*(rs*fluid.rhoGS + fluid.rhoOS);
    else
        bO  = fluid.bO(p);
        if isfield(fluid, 'BOxmuO')
            muO = fluid.BOxmuO(p).*bO;
        else
            muO = fluid.muO(p);
        end
        rhoO   = bO.*fluid.rhoOS;
    end
        
    if any(bO < 0)
        warning('Negative oil compressibility present!')
    end
    
    rhoOf  = s.faceAvg(rhoO);
    dpO    = s.Grad(p) - rhoOf.*gdz;
    % oil upstream-index
    upco = (value(dpO)<=0);
    
    [krOf, krO] = s.splitFaceCellValue(s, upco, krO);
    [muOf, muO] = s.splitFaceCellValue(s, upco, muO);
    mobO   = krO./muO;
    
    vO   = -(krOf./muOf).*T.*dpO;
end

