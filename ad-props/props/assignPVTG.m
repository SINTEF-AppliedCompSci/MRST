function f = assignPVTG(f, pvtg, reg)
f.BG  = @(pg, rv, flag, varargin)BG(pg, rv, pvtg, flag, reg, varargin{:});
f.bG  = @(pg, rv, flag, varargin)bG(pg, rv, pvtg, flag, reg, varargin{:});
f.muG = @(pg, rv, flag, varargin)muG(pg, rv, pvtg, flag, reg, varargin{:});
f.rvSat = @(pg, varargin)rvSat(pg, pvtg, reg, varargin{:});
end

function v = BG(pg, rv, pvtg, flag, reg, varargin)
pvtinx = getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvtg;
for k = 1:numel(T), T{k}.data = T{k}.data(:,1:2); end
v = interpRegPVT(T, rv, pg, flag, pvtinx);
end

function v = bG(pg, rv, pvtg, flag, reg, varargin)
pvtinx = getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvtg;
for k = 1:numel(T), T{k}.data = [T{k}.data(:,1), 1./T{k}.data(:,2)]; end
v = interpRegPVT(T, rv, pg, flag, pvtinx);
end

function v = muG(pg, rv, pvtg, flag, reg, varargin)
pvtinx = getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvtg;
for k = 1:numel(T), T{k}.data = T{k}.data(:,[1 3]); end
v = interpRegPVT(T, rv, pg, flag, pvtinx);
end

function v = rvSat(pg, pvtg, reg, varargin)
pvtinx = getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});
T = cellfun(@(x)[x.key x.data(x.pos(1:end-1),1)], pvtg, 'UniformOutput', false);
v = interpReg(T, pg, pvtinx);
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
