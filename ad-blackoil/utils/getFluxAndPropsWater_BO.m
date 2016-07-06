function [vW, bW, mobW, rhoW, pW, upcw, dpW] = getFluxAndPropsWater_BO(model, pO, sW, krW, T, gdz)
    fluid = model.fluid;
    s = model.operators;
    % Check for capillary pressure (p_cow)
    pcOW = 0;
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW);
    end
    pW = pO - pcOW;
    
    bW     = fluid.bW(pO);
    rhoW   = bW.*fluid.rhoWS;
    % rhoW on face, average of neighboring cells
    rhoWf  = s.faceAvg(rhoW);
    mobW   = krW./fluid.muW(pO);
    dpW    = s.Grad(pO-pcOW) - rhoWf.*gdz;
    % water upstream-index
    upcw  = (double(dpW)<=0);
    vW = -s.faceUpstr(upcw, mobW).*T.*dpW;
    if any(bW < 0)
        warning('Negative water compressibility present!')
    end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
