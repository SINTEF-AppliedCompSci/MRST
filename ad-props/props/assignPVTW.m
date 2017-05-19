function f = assignPVTW(f, pvtw, reg)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    f.cW   = pvtw(1, 3);
    f.muWr = pvtw(1, 4);
else
    f.cW   = pvtw(reg.PVTNUM, 3);
    f.muWr = pvtw(reg.PVTNUM, 4);
end
f.BW   = @(pw, varargin)BW(pw, pvtw, reg, varargin{:});
f.bW   = @(pw, varargin)bW(pw, pvtw, reg, varargin{:});
f.muW  = @(pw, varargin)muW(pw, pvtw, reg, varargin{:});
end


function v = BW(pw, pvtw, reg, varargin)
v = 1./bW(pw, pvtw, reg, varargin{:});
end

function v = bW(pw, pvtw, reg, varargin)
pvtnum= getPVTNUM(pw, reg, varargin{:});


pwr  = pvtw(pvtnum,1); % ref pres
bwr  = pvtw(pvtnum,2); % ref fvf
cw   = pvtw(pvtnum,3); % compress
X = cw.*(pw-pwr);
v = (1+X+X.^2/2)./bwr;
end

function v = muW(pw, pvtw, reg, varargin)
pvtnum= getPVTNUM(pw, reg, varargin{:});

pwr  = pvtw(pvtnum,1); % ref pres
muwr = pvtw(pvtnum,4); % ref visc
vbw  = pvtw(pvtnum,5); % viscosibility
Y = -vbw.*(pw-pwr);
v = muwr./(1+Y+Y.^2/2);
end

function pvtnum= getPVTNUM(pw, reg, varargin)
pvtinx = getRegMap(pw, reg.PVTNUM, reg.PVTINX, varargin{:});

if strcmp(pvtinx{1}, ':')
   pvtnum=ones(size(pw));
   assert(numel(pvtinx)==1);
else
    pvtnum=nan(size(pw));
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
