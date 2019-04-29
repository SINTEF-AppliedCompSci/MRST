function compi = crossFlowMixture(flux, compi, map, conserveMass)
    if nargin < 4
        conserveMass = false;
    end
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
    sum_in = sum_perf(flux_in, map);
    top_in =  bsxfun(@times, net_injection, compi);
    comp = sum_in + top_in;
    compT = sum(comp, 2);
    % Normalize to get fractions
    comp = bsxfun(@rdivide, comp, compT);
    if conserveMass
        % Ensure exact re-injection with "fractions" which are not in unit
        % range for injectors
        act = sum(top_in, 2) > 0;
        % Top net injected composition + sum of composition into well-bore,
        % divided by total flux out from well-bore
        tmp = (top_in(act, :) + sum_in(act, :))./sum(max(net_flux(act, :), 0), 2);
        compi(act, :) = tmp;
    end
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