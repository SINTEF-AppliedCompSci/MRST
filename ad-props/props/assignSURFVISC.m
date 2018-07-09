function f = assignSURFVISC(f, surfvisc, reg)
f.muWSft = @(c, varargin) muWSft(c, surfvisc, reg, varargin{:});
end

function v = muWSft(c, surfvisc, reg, varargin)
satinx = getRegMap(c, reg.PVTNUM, reg.PVTINX, varargin{:});
surfvisc = extendTab(surfvisc);
v = interpReg(surfvisc, c, satinx);
end

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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
