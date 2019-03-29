function compi = crossFlowMixture(flux, compi, map)
    % Into wellbore
    flux_in = -min(flux, 0);
    if all(flux_in == 0)
        return
    end

    % Find net flux - this is what is possibily injected from the
    % surface connection with the top composition
    net_flux = sum_perf(sum(flux, 2), map);
    net_injection = max(net_flux, 0);
    % Flux into well-bore plus net injection weighted with topside
    % composition
    comp = sum_perf(flux_in, map) + bsxfun(@times, net_injection, compi);

    % Normalize to get fractions
    compT = sum(comp, 2);
    comp = bsxfun(@rdivide, comp, compT);

    active = compT > 0;
    compi(active, :) = comp(active, :);
end

function out = sum_perf(v, map)
    nph = size(v, 2);
    nw = numel(map.W);
    out = zeros(nw, nph);
    for i = 1:nph
        out(:, i) = accumarray(map.perf2well, v(:, i));
    end
end