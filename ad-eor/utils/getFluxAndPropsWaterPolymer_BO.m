function [vW, vP, bW, muWeffMult, mobW, mobP, rhoW, pW, upcw, a] = getFluxAndPropsWaterPolymer_BO(model, pO, sW, c, ads, krW, T, gdz)
%
%
% SYNOPSIS:
%   function [vW, vP, bW, muWeffMult, mobW, mobP, rhoW, pW, upcw, a] = getFluxAndPropsWaterPolymer_BO(model, pO, sW, c, ads, krW, T, gdz)
%
% DESCRIPTION: Given pressure, saturation and polymer concentration and some
% other input variables, compute the fluxes and other properties, as listed
% below. Used in the assembly of the blackoil equations with polymer.

%
% PARAMETERS:
%   model - Instance of the model
%   pO    - Pressure in oil phase
%   sW    - Saturation
%   c     - Polymer concentration
%   ads   - Adsorption value
%   krW   - Water relative permeability
%   T     - Transmissibilities
%   gdz   - z-gradient (to be used in computation of gravitational flux)
%
% RETURNS:
%   vW         - Water flux
%   vP         - Polymer flux
%   bW         - Water surface volume factor
%   muWeffMult - Effective viscosity multiplier 
%   mobW       - Water mobility
%   mobP       - Polymer mobility
%   rhoW       - Water density
%   pW         - Pressure in water phase
%   upcw       - Water upstream direction
%   a          - Coefficient appearing in the polymer model 
%
% EXAMPLE:
%
% SEE ALSO: equationsThreePhaseBlackOilPolymer
%

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW);
    end
    pW = pO - pcOW;

    % Multipliers due to polymer
    mixpar = fluid.mixPar;
    cbar   = c/fluid.cmax;
    a = fluid.muWMult(fluid.cmax).^(1-mixpar);
    b = 1./(1-cbar+cbar./a);
    % The viscosity multiplier only result from the polymer mixing.
    muWeffMult = b.*fluid.muWMult(c).^mixpar;
    permRed = 1 + ((fluid.rrf-1)./fluid.adsMax).*ads;
    muWMult  = muWeffMult.*permRed;

    % Water props
    bW     = fluid.bW(pO);
    rhoW   = bW.*fluid.rhoWS;
    % rhoW on face, average of neighboring cells
    rhoWf  = s.faceAvg(rhoW);
    muW    = fluid.muW(pO);
    muWeff = muWMult.*muW;
    mobW   = krW./muWeff;
    dpW    = s.Grad(pW) - rhoWf.*gdz;
    % water upstream-index
    upcw = (double(dpW)<=0);
    vW   = -s.faceUpstr(upcw, mobW).*T.*dpW;
    if any(bW < 0)
        warning('Negative water compressibility present!')
    end

    % Polymer
    mobP = mobW./(a + (1-a)*cbar);
    vP   = - s.faceUpstr(upcw, mobP.*c).*s.T.*dpW;

end


