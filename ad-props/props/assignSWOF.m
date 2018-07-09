function f = assignSWOF(f, swof, reg)
f.krW  = @(sw, varargin)krW(sw, swof, reg, varargin{:});
f.krOW = @(so, varargin)krOW(so, swof, reg, varargin{:});
f.pcOW = @(sw, varargin)pcOW(sw, swof, reg, varargin{:});
swcon  = cellfun(@(x)x(1,1), swof);
ntsat = numel(reg.SATINX);
if ntsat == 1
    f.sWcon = swcon(1);
else
    f.sWcon = swcon(reg.SATNUM);
end

if isfield(reg, 'SURFNUM')
   % Assign miscible relperm for surfactant
   f.krWSft  = @(sw, varargin)krWSft(sw, swof, reg, varargin{:});
   f.krOWSft  = @(so, varargin)krOWSft(so, swof, reg, varargin{:});
   % Assign residual water saturation for surfactant
   f.sWconSft = swcon(reg.SURFNUM);
   % Assign residual oil saturation
   sOres  = cellfun(@(x)x(end, 1), swof);
   f.sOres = 1 - sOres(reg.SATNUM);
   f.sOresSft = 1 - sOres(reg.SURFNUM);
end

end

function v = krW(sw, swof, reg, varargin)
satinx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, satinx);
end

function v = krOW(so, swof, reg, varargin)
satinx = getRegMap(so, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, 1 - so, satinx);
end

function v = pcOW(sw, swof, reg, varargin)
satinx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,4]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, satinx);
end

function v = krWSft(sw, swof, reg, varargin)
surfinx = getRegMap(sw, reg.SURFNUM, reg.SURFINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, surfinx);
end

function v = krOWSft(so, swof, reg, varargin)
surfinx = getRegMap(so, reg.SURFNUM, reg.SURFINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, 1 - so, surfinx);
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

