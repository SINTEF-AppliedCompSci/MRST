function [pvMult, transMult, mobMult, pvMult0, transMult0, mobMult0] = getMultipliers(fluid, p, p0)
%Get dynamic multiplier values for reservoir quantities
%
% SYNOPSIS:
%   [pvMult, transMult, mobMult] = getMultipliers(fluid, p, p0)
%
% REQUIRED PARAMETERS:
%   fluid    - Fluid model, typically from initDeckADIFluid or some other
%              constructor. 
%
%   p        - Pressure per cell.
%
%   p0       - Pressure per cell (previous timestep).
%
% RETURNS:
%   pvMult   - Pore volume multiplier per cell. Multiplicative modifier for
%              the pore volume based on a simplified model for how the
%              pores grow when fluid pressure is increasing.
%
%  transMult - Transmissibility multiplier, modelling pressure dependent
%              permeability. One value per interface.
%
%  mobMult   - Mobility multiplier per cell.
%
%  pvMult0, transMult0, mobMult0 - Multipliers for previous timestep.
%

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
    [pvMult, transMult, mobMult, pvMult0, transMult0, mobMult0] = deal(1);
    % Pressure dependent pore volume multiplier
    if isfield(fluid, 'pvMultR')
        pvMult =  fluid.pvMultR(p);
        if nargout > 3
            pvMult0 = fluid.pvMultR(p0);
        end
    end
    % Pressure dependent mobility multiplier 
    if isfield(fluid, 'tranMultR')
        mobMult = fluid.tranMultR(p);
        if nargout > 4
            mobMult0 = fluid.tranMultR(p0);
        end
    end
    % Pressure dependent transmissibility multiplier
    if isfield(fluid, 'transMult')
       transMult = fluid.transMult(p); 
       if nargout > 5
           transMult0 = fluid.transMult(p0); 
       end
    end
end
