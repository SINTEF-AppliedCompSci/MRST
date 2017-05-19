function f = assignPVCDO(f, pvcdo, reg)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    f.cO  = pvcdo(1, 3);
else
    f.cO  = pvcdo(reg.PVTNUM, 3);
end
f.BO     = @(po, varargin)BO(po, pvcdo, reg, varargin{:});
f.bO     = @(po, varargin)bO(po, pvcdo, reg, varargin{:});
f.BOxmuO = @(po, varargin)BOxmuO(po, pvcdo, reg, varargin{:});

f.muO = @(po, varargin) bO(po, pvcdo, reg, varargin{:}).*...
                    BOxmuO(po, pvcdo, reg, varargin{:});
end

function v = BO(po, pvcdo, reg, varargin)
v = 1./bO(po, pvcdo, reg, varargin{:});
end

function v = bO(po, pvcdo, reg, varargin)
pvtnum = getPVTNUM(po, reg, varargin{:});

por  = pvcdo(pvtnum,1); % ref pres
bor  = pvcdo(pvtnum,2); % ref fvf
co   = pvcdo(pvtnum,3); % compress
X = co.*(po-por);
v = exp(X)./bor;
end

function v = BOxmuO(po, pvcdo, reg, varargin)
pvtnum = getPVTNUM(po, reg, varargin{:});

por  = pvcdo(pvtnum,1); % ref pres
bor  = pvcdo(pvtnum,2); % ref fvf
co   = pvcdo(pvtnum,3); % compress
muor = pvcdo(pvtnum,4); % ref visc
vbo  = pvcdo(pvtnum,5); % viscosibility
Y = (co-vbo).*(po-por);
v = bor.*muor.*exp(-Y);
end


function pvtnum= getPVTNUM(po, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});

if strcmp(pvtinx{1}, ':')
   pvtnum=ones(size(po));
   assert(numel(pvtinx)==1);
else
    pvtnum=nan(size(po));
    for i=1:numel(pvtinx)
       pvtnum(pvtinx{i})=i;
    end
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


