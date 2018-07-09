function f = assignSOF3(f, sof3, reg)
f.krOW = @(so, varargin)krOW(so, sof3, reg, varargin{:});
f.krOG = @(so, varargin)krOG(so, sof3, reg, varargin{:});
end

function v = krOW(so, sof3, reg, varargin)
satinx = getRegMap(so, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), sof3, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, so, satinx);
end

function v = krOG(so, sof3, reg, varargin)
satinx = getRegMap(so, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), sof3, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, so, satinx);
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

