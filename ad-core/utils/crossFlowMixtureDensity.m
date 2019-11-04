function rhoMix = crossFlowMixtureDensity(massFlux, volumeTotalFlux, massFluxFromSurface, map)
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