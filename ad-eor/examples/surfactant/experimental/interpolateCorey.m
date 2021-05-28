function [s, krW, krO] = interpolateCorey()
%Undocumented Utility Function

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

    n      = 3;
    swc    = 0.2;
    sor    = 0.2;
    krOres = 0.5;
    krWres = 0.6;
    n      = 1.5;
    swc    = 0.05;
    sor    = 0.05;
    krOres = 1;
    krWres = 1;

    s = (0:0.01:1)';
    [krW, krO] = corey(s, n, swc, sor, krOres, krWres);

end

function [krW, krO] = corey(s, n, swc, sor, krOres, krWres)
    s = min(max(s, swc), 1 - sor);
    krW = krWres*((s - swc)./(1 - swc - sor)).^n;
    sO = 1 - s;
    krO = krOres*((sO - sor)./(1 - swc - sor)).^n;
end
