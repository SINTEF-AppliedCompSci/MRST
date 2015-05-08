% Files
%   assignDENSITY.m    - dens of size ntpvtx3
%   assignPVCDO.m      - por  = pvcdo(pvtnum,1);  ref pres
%   assignPVDO.m       - f.muO = @(po, varargin)muO(po, pvdo, reg, varargin{:});
%   assignPVTW.m       - pwr  = pvtw(pvtnum,1);  ref pres
%   assignRelPerm.m    - if ~isfield(f, 'krOG')       two-phase water/oil
%   assignSOF3.m       - f.relperm3ph = @(sw, sg, varargin)relperm3ph(sw, sg, f, varargin);
%   initDeckADIFluid.m - props

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
