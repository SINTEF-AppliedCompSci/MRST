function [vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO_BCP(...
        model, pO, p_prop, sW, krW, T, gdz, bcp)
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

    fluid = model.fluid;
    s = model.operators;
    % Check for capillary pressure (p_cow)
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW);
    end
    pW = p_prop - pcOW;
    
    bW     = fluid.bW(p_prop);
    rhoW   = bW.*fluid.rhoWS;
    % rhoW on face, average of neighboring cells
    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./fluid.muW(p_prop);
    dpW    = s.Grad(pO-pcOW) - rhoWf.*gdz;
    
    % Periodic boundary conditions
    if ~isempty(bcp)
        dpcOW = 0;
        if isfield(fluid, 'pcOW') && ~isempty(sW)
            dpcOW = pcOW(s.N(bcp.face,1)) - pcOW(s.N(bcp.face,2));
        end
        dpW(bcp.face) = dpW(bcp.face) + bcp.value + dpcOW;
    end
    
    % water upstream-index
    upcw  = (value(dpW)<=0);
    vW = -s.faceUpstr(upcw, mobW).*T.*dpW;
    if any(bW < 0)
        warning('Negative water compressibility present!')
    end
end
