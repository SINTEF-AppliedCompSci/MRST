function rhoMix = crossFlowMixtureDensity(massFlux, volumeTotalFlux, massFluxFromSurface, map)
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

    massIntoWell = sum_perf(max(-massFlux, 0), map);
    volumeIntoWell = sum_perf(max(-volumeTotalFlux, 0), map);
    volumeFromWell = sum_perf(max(volumeTotalFlux, 0), map);
    % Net volume - this is the volume corresponding to the massFluxTop
    % variable
    volumeExchangeSurface = volumeIntoWell - volumeFromWell;
    
    totalMassIntoWell = (massIntoWell + max(massFluxFromSurface, 0));
    bad = sum(totalMassIntoWell, 2) == 0 | volumeFromWell == 0;
    rhoMix = totalMassIntoWell./(volumeFromWell + max(volumeExchangeSurface, 0));
    rhoMix(bad, :) = 1;
end

function out = sum_perf(v, map)
    nph = size(v, 2);
    nw = numel(map.W);
    out = zeros(nw, nph);
    for i = 1:nph
        out(:, i) = accumarray(map.perf2well, v(:, i));
    end
end
