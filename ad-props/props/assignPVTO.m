function f = assignPVTO(f, pvto, reg)
f.bO  = @(po, rs, flag, varargin)bO(po, rs, pvto, flag, reg, varargin{:});
f.muO = @(po, rs, flag, varargin)muO(po, rs, pvto, flag, reg, varargin{:});
f.rsSat = @(po, varargin)rsSat(po, pvto, reg, varargin{:});
end

function v = bO(po, rs, pvto, flag, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvto;
for k = 1:numel(T), T{k}.data = [T{k}.data(:,1), 1./T{k}.data(:,2)]; end
v = interpRegPVT(T, po, rs, flag, pvtinx);
end

function v = muO(po, rs, pvto, flag, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvto;
for k = 1:numel(T), T{k}.data = T{k}.data(:,[1 3]); end
v = interpRegPVT(T, po, rs, flag, pvtinx);
end

function v = rsSat(po, pvto, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = cellfun(@(x)[x.data(x.pos(1:end-1),1) x.key], pvto, 'UniformOutput', false);
v = interpReg(T, po, pvtinx);
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
