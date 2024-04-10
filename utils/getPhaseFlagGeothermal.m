function [isL, isV, isLV] = getPhaseFlagGeothermal(p, h, hL, hV, pcrit)
%Get phase flag indicating if the cell contains liqiud, vapor, or a mixture

    isL  = value(h) < hL;% & ~isSC; % Liquid
    isV  = value(h) > hV;% & ~isSC; % Vapor
    isSC = value(p) > pcrit;      % Supercritical
    isV  = isV | (isSC & ~isL);            % Label supercritical as vapor
    isLV = ~isL & ~isV;           % Two-phase
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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