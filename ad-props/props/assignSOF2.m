function f = assignSOF2(f, sof2, reg)
f.krO  = @(so, varargin)krO(so, sof2, reg, varargin{:});

if isfield(reg, 'SURFNUM')
   % Assign miscible relperm for surfactant
   f.krOSft  = @(so, varargin)krOSft(so, sof2, reg, varargin{:});
end

end

function v = krO(so, sof2, reg, varargin)
satinx = getRegMap(so, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), sof2, 'UniformOutput', false);
v = interpReg(T, so, satinx);
end

function v = krOSft(so, sof2, reg, varargin)
surfinx = getRegMap(so, reg.SURFNUM, reg.SURFINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), sof2, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, so, surfinx);
end


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
