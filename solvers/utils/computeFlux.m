function flux = computeFlux(G, u, T)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    dispif(mrstVerbose, 'computeFlux\n');
    flux = zeros(numelValue(T{1}), 1);
    s = getSampleAD(u, T{:});
    if isa(s, 'ADI')
        flux = double2ADI(flux, s);
    end
    ind = all(G.faces.neighbors ~= 0, 2);
    c1 = G.faces.neighbors(ind, 1);
    c2 = G.faces.neighbors(ind, 2);
    flux(ind) = T{1}(ind) .* u(c1) - T{2}(ind) .* u(c2);

    c1 = max(G.faces.neighbors(~ind, :), [], 2);
    flux(~ind) = T{1}(~ind) .* u(c1) - T{2}(~ind);

    ind = G.faces.neighbors(:, 1) == 0;
    flux(ind) = -flux(ind);
end
