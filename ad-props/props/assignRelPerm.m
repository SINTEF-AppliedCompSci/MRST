function f = assignRelPerm(f)

if ~isfield(f, 'krOG')      % two-phase water/oil
   f.relPerm = @(sw, varargin)relPermWO(sw, f, varargin{:});
   if isfield(f, 'krWSft') % Surfactant is used
      f.relPermSft = @(sw, varargin)relPermWOSft(sw, f, varargin{:});
   end
elseif ~isfield(f, 'krOW')  % two-phase oil/gas
    f.relPerm = @(sg, varargin)relPermOG(sg, f, varargin{:});
else                        % three-phase
    f.relPerm = @(sw, sg, varargin)relPermWOG(sw, sg, f, varargin{:});
end

end

function [krW, krO] = relPermWO(sw, f, varargin)
krW = f.krW(sw, varargin{:});
if isfield(f, 'krO')
    krO = f.krO(1-sw, varargin{:});
else
    krO = f.krOW(1-sw, varargin{:});
end
end


function [krWSft, krOSft] = relPermWOSft(sw, f, varargin)
krWSft = f.krWSft(sw, varargin{:});
if isfield(f, 'krO')
    krOSft = f.krOSft(1-sw, varargin{:});
else
    krOSft = f.krOWSft(1-sw, varargin{:});
end
end

function [krO, krG] = relPermOG(sg, f, varargin)
krG = f.krG(sg, varargin{:});
if isfield(f, 'krO')
    krO = f.krO(1-sg, varargin{:});
else
    krO = f.krOG(1-sg, varargin{:});
end
end

function [krW, krO, krG] = relPermWOG(sw, sg, f, varargin)
if isfield(f, 'sWcon')
    swcon = f.sWcon;
else
    swcon = 0;
end
swcon = min(swcon, value(sw)-1e-5);

d  = (sg+sw-swcon);
ww = (sw-swcon)./d;
%krW = ww.*f.krW(sg+sw, varargin{:});
krW = f.krW(sw, varargin{:});

wg = 1-ww;
%krG = wg.*f.krG(sg+sw-swcon, varargin{:});
krG = f.krG(sg, varargin{:});

so = 1-sw-sg;
krow = f.krOW(so, varargin{:});
krog = f.krOG(so,  varargin{:});
krO  = wg.*krog + ww.*krow;
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
