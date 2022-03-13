function [vO, bO, mobO, rhoO, p, upco, dpO] = getFluxAndPropsOil_BO_BCP(...
        model, p, p_prop, sO, krO, T, gdz, rs, isSat, bcp)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
        bO  = fluid.bO(p_prop,  rs, isSat);
        muO = fluid.muO(p_prop, rs, isSat);
        rhoO   = bO.*(rs*fluid.rhoGS + fluid.rhoOS);
    else
        bO  = fluid.bO(p_prop);
        if isfield(fluid, 'BOxmuO')
            muO = fluid.BOxmuO(p_prop).*bO;
        else
            muO = fluid.muO(p_prop);
        end
        rhoO   = bO.*fluid.rhoOS;
    end
    
    if any(bO < 0)
        warning('Negative oil compressibility present!')
    end
    
    rhoOf  = s.faceAvg(rhoO);
    mobO   = krO./muO;
    dpO    = s.Grad(p) - rhoOf.*gdz;
    
    % Periodic boundary conditions
    if ~isempty(bcp)
        dpO(bcp.face) = dpO(bcp.face) + bcp.value;
    end
    
    % oil upstream-index
    upco = (value(dpO)<=0);
    vO   = - s.faceUpstr(upco, mobO).*T.*dpO;
end
