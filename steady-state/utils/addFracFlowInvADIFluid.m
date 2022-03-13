function fluid = addFracFlowInvADIFluid(fluid, deck, varargin)
% Add a fluid function for computing the inverse of the water fractional
% flow curve. This may be called later by using
%   sW = fluid.fracFlowInv(ff)
% The input deck must have the properties SWOF, PVTW and PVCDO.
% 
% Fractional flow = (krW/muW) / ( (krW/muW) + (krO/muO) )

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

T = createFracFlowTablesFromDeck(deck);
T = extendTab(T);

% Add to fluid
reg   = handleRegions(deck, varargin{:});
fluid = assignFracFlowInv(fluid, T, reg);

end

% Remaining code is a slight modification of assignSWOF

function f = assignFracFlowInv(f, T, reg)
f.fracFlowInv = @(ff, varargin)FracFlowInv(ff, T, reg, varargin{:});
end


function v = FracFlowInv(ff, T, reg, varargin)
satinx = getRegMap(ff, reg.SATNUM, reg.SATINX, varargin{:});
v = interpReg(T, ff, satinx);
end
