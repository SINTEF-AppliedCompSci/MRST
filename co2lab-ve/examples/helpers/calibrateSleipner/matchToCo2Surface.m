function val = matchToCo2Surface(sG,surface,G,fluid)
% "surface.h" is permitted to contain nan's (i.e., missing data), in which
% case the following calculation omits nans in order to quantify the match
% between simulated heights and (present) data heights

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    h = G.cells.H(~isnan(surface.h)) .* sG(~isnan(surface.h)) / (1-fluid.res_water) ...
        - surface.h(~isnan(surface.h));
    val = sum( h.^2 .* G.cells.volumes(~isnan(surface.h)) );
    
    %h = G.cells.H.*sG/(1-fluid.res_water)-surface.h;
    %val=sum(h.^2.*G.cells.volumes, 'omitnan'); %@@ doesn't work for ADIs
end
